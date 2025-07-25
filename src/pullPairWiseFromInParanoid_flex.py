import pandas as pd
import requests
from io import StringIO
from itertools import product

### IMPORTANT ####
# Warning: First run "src/convertOrthUniprotToEnsembl.r" in R
##################

# === USER INPUT ===
data_dir = "data"
# original species
orig_species_input = "human"
# Ortholog species
species_input = "frog"  # Options: "mouse", "zebrafish", # No sheep in inParanoid

# species_to_process = ["mouse", "rat", "human", "zebrafish", "chimpanzee", "chicken", "pig", "cow", "dog", "horse", "marmoset",   "macaque"]

# === INTERNAL MAPPINGS ===

orig_species_info = {
    "mouse":         {"taxid": "10090", "code": "mmusculus"},
    "rat":           {"taxid": "10116", "code": "rnorvegicus"},
    "zebrafish":     {"taxid": "7955",  "code": "drerio"},
    "chimpanzee":    {"taxid": "9598",  "code": "ptroglodytes"},
    "chicken":       {"taxid": "9031",  "code": "ggallus"},
    "pig":           {"taxid": "9823",  "code": "sscrofa"},
    "cow":           {"taxid": "9913",  "code": "btaurus"},
    "dog":           {"taxid": "9615",  "code": "clfamiliaris"},
    "horse":         {"taxid": "9796",  "code": "ecaballus"},
    #"sheep":         {"taxid": "9940",  "code": "oarambouillet"},
    "marmoset":      {"taxid": "9483",  "code": "cjacchus"},
    "macaque":       {"taxid": "9544",  "code": "mmulatta"},
    "frog":          {"taxid": "8364",  "code": "xtropicalis"},
    "rabbit":        {"taxid": "9986",  "code": "ocuniculus"},
    "guineapig":     {"taxid": "10141",  "code": "cporcellus"},
    "pufferfish":    {"taxid": "99883",  "code": "tnigroviridis"},
    "human":         {"taxid": "9606",  "code": "hsapiens"},
}


species_info = {
    "mouse":         {"taxid": "10090", "code": "mmusculus"},
    "rat":           {"taxid": "10116", "code": "rnorvegicus"},
    "zebrafish":     {"taxid": "7955",  "code": "drerio"},
    "chimpanzee":    {"taxid": "9598",  "code": "ptroglodytes"},
    "chicken":       {"taxid": "9031",  "code": "ggallus"},
    "pig":           {"taxid": "9823",  "code": "sscrofa"},
    "cow":           {"taxid": "9913",  "code": "btaurus"},
    "dog":           {"taxid": "9615",  "code": "clfamiliaris"},
    "horse":         {"taxid": "9796",  "code": "ecaballus"},
    #"sheep":         {"taxid": "9940",  "code": "oarambouillet"},
    "marmoset":      {"taxid": "9483",  "code": "cjacchus"},
    "macaque":       {"taxid": "9544",  "code": "mmulatta"},
    "frog":          {"taxid": "8364",  "code": "xtropicalis"},
    "rabbit":        {"taxid": "9986",  "code": "ocuniculus"},
    "guineapig":     {"taxid": "10141",  "code": "cporcellus"},
    "pufferfish":    {"taxid": "99883",  "code": "tnigroviridis"},
    "human":         {"taxid": "9606",  "code": "hsapiens"},
}

if orig_species_input not in orig_species_info:
    raise ValueError(f"Species '{orig_species_input}' not supported. Choose from: {list(orig_species_info)}")

if species_input not in species_info:
    raise ValueError(f"Species '{species_input}' not supported. Choose from: {list(species_info)}")

species = species_info[species_input]
taxid = species["taxid"]
code = species["code"]

orig_species = orig_species_info[orig_species_input]
orig_taxid = orig_species["taxid"]
orig_code = orig_species["code"]

# === Step 1: Download inParanoid file for human vs species ===
url = f"https://inparanoidb.sbc.su.se/download/sqltable/{orig_taxid}&{taxid}&prot"
r = requests.get(url)
r.raise_for_status()

df = pd.read_csv(StringIO(r.text.strip()), sep="\t", header=None)
df.columns = ["cluster_id", "bitscore", "source_file", "inparalog_score", "protein_id", "seed_score"]

# === Step 2: Add species labels ===
def infer_species(src):
    if orig_taxid in src:
        return orig_species_input
    elif taxid in src:
        return species_input
    return "unknown"

df["species"] = df["source_file"].apply(infer_species)

# === Step 3: Build ortholog pairs ===
records = []
for cid, grp in df.groupby("cluster_id"):
    orig_spec = grp[grp["species"] == orig_species_input]
    others = grp[grp["species"] == species_input]
    
    for h, o in product(orig_spec.itertuples(index=False), others.itertuples(index=False)):
        records.append({
            "cluster_id": cid,
            f"{orig_species_input}_protein": h.protein_id,
            f"{orig_species_input}_inparalog_score": h.inparalog_score,
            f"{orig_species_input}_seed_score": h.seed_score,
            f"{species_input}_protein": o.protein_id,
            f"{species_input}_inparalog_score": o.inparalog_score,
            f"{species_input}_seed_score": o.seed_score,
            "bitscore": (h.bitscore + o.bitscore) / 2
        })


df_orthologs = pd.DataFrame(records)
df_orthologs.to_csv(f"data/{orig_species_input}_centric_inParanoid_{species_input}.csv", index=False)

# === Step 4: UniProt → Gene Name ===
print("Starting annotation process...")

# --- Load Human UniProt Mapping once ---
orig_uniprot_file = os.path.join(data_dir, f"uniprotMapping_{orig_species_input}.csv")
if not os.path.exists(orig_uniprot_file):
    print(f"Error: {orig_species_input} UniProt mapping file not found at {orig_uniprot_file}. Please run the UniProt fetching script for human first.")
    exit() # Exit if the essential human file is missing

print(f"Loading {orig_species_input} UniProt mapping from {orig_uniprot_file}...")
# We only need 'Entry' (Accession) and 'Gene Names'
orig_uniprot_df = pd.read_csv(orig_uniprot_file, usecols=["Entry", "Gene Names"])
# Rename 'Gene Names' to distinguish it as human gene name
orig_uniprot_df = orig_uniprot_df.rename(columns={"Gene Names": f"{orig_species_input}_Gene_Name"})
print(f"Loaded {len(orig_uniprot_df)} {orig_species_input} UniProt entries.")

# --- Process each species ---
print(f"\n--- Processing {species_input} ---")

inparanoid_file = os.path.join(data_dir, f"{orig_species_input}_centric_inParanoid_{species_input}.csv")

species_uniprot_file = os.path.join(data_dir, f"uniprotMapping_{species_input}.csv")
output_file = os.path.join(data_dir, f"{orig_species_input}_centric_inParanoid_{species_input}_AnnWithUniProt.csv")

# Check if input files exist
if not os.path.exists(inparanoid_file):
    print(f"Skipping {species_input}: InParanoid file not found at {inparanoid_file}")
    if not os.path.exists(species_uniprot_file):
        print(f"Skipping {species_input}: UniProt mapping file not found at {species_uniprot_file}")

print(f"Loading InParanoid data from {inparanoid_file}...")
df_inparanoid = pd.read_csv(inparanoid_file)
print(f"Loaded {len(df_inparanoid)} InParanoid entries for {species_input}.")

print(f"Loading {species_input} UniProt mapping from {species_uniprot_file}...")
species_uniprot_df = pd.read_csv(species_uniprot_file, usecols=["Entry", "Gene Names"])
# Rename 'Gene Names' to distinguish it as the species' gene name
species_uniprot_df = species_uniprot_df.rename(columns={"Gene Names": f"{species_input}_Gene_Name"})
print(f"Loaded {len(species_uniprot_df)} {species_input} UniProt entries.")

# 1. Annotate based on "orig_species_protein" with orig UniProt data
print(f"Merging {orig_species_input} gene names...")
# Use left merge to keep all rows from df_inparanoid
df_merged = pd.merge(
    df_inparanoid,
    orig_uniprot_df,
    left_on=f"{orig_species_input}_protein",
    right_on="Entry",
    how="left"
)
    # Drop the redundant 'Entry' column from the merge
df_merged = df_merged.drop(columns=["Entry"])

# 2. Annotate based on "{species_input}_protein" with species UniProt data
print(f"Merging {species_input} gene names...")
df_merged = pd.merge(
    df_merged, # Merge into the already merged dataframe
    species_uniprot_df,
    left_on=f"{species_input}_protein",
    right_on="Entry",
    how="left"
)
# Drop the redundant 'Entry' column from the second merge
df_merged = df_merged.drop(columns=["Entry"])

# 3. Save merged data
print(f"Saving merged data to {output_file}...")
df_merged.to_csv(output_file, index=False)
print(f"Successfully saved {len(df_merged)} annotated entries for {species_input}.")

print("\nAnnotation process completed for all specified species.")