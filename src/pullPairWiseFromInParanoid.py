import pandas as pd
import requests
from io import StringIO
from itertools import product

### IMPORTANT ####
# Warning: First run "src/convertOrthUniprotToEnsembl.r" in R
##################

# === USER INPUT ===
species_input = "horse"  # Options: "mouse", "zebrafish", # No sheep in inParanoid

# === INTERNAL MAPPINGS ===
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
    "sheep":         {"taxid": "9940",  "code": "oarambouillet"},
    "marmoset":      {"taxid": "9483",  "code": "cjacchus"},
    "macaque":       {"taxid": "9544",  "code": "mmulatta"},
}

if species_input not in species_info:
    raise ValueError(f"Species '{species_input}' not supported. Choose from: {list(species_info)}")

species = species_info[species_input]
taxid = species["taxid"]
code = species["code"]

# === Step 1: Download inParanoid file for human vs species ===
url = f"https://inparanoidb.sbc.su.se/download/sqltable/9606&{taxid}&prot"
r = requests.get(url)
r.raise_for_status()

df = pd.read_csv(StringIO(r.text.strip()), sep="\t", header=None)
df.columns = ["cluster_id", "bitscore", "source_file", "inparalog_score", "protein_id", "seed_score"]

# === Step 2: Add species labels ===
def infer_species(src):
    if "9606" in src:
        return "human"
    elif taxid in src:
        return species_input
    return "unknown"

df["species"] = df["source_file"].apply(infer_species)

# === Step 3: Build ortholog pairs ===
records = []
for cid, grp in df.groupby("cluster_id"):
    humans = grp[grp["species"] == "human"]
    others = grp[grp["species"] == species_input]
    for h, o in product(humans.itertuples(index=False), others.itertuples(index=False)):
        records.append({
            "cluster_id": cid,
            "human_protein": h.protein_id,
            "human_inparalog_score": h.inparalog_score,
            "human_seed_score": h.seed_score,
            f"{species_input}_protein": o.protein_id,
            f"{species_input}_inparalog_score": o.inparalog_score,
            f"{species_input}_seed_score": o.seed_score,
            "bitscore": (h.bitscore + o.bitscore) / 2
        })

df_orthologs = pd.DataFrame(records)
df_orthologs.to_csv(f"data/inParanoid_{species_input}.csv", index=False)

# === Step 4: Human UniProt → HGNC/Ensembl mapping ===
hgnc_df = pd.read_csv("data/HGNC_gene_info_full.tsv", sep="\t", dtype=str)
hgnc_df = hgnc_df.dropna(subset=["uniprot_ids", "ensembl_gene_id"])
hgnc_exploded = hgnc_df.assign(uniprot_id=hgnc_df["uniprot_ids"].str.split(",")).explode("uniprot_id")
hgnc_exploded["uniprot_id"] = hgnc_exploded["uniprot_id"].str.strip()

df_merged = df_orthologs.merge(
    hgnc_exploded[["uniprot_id", "ensembl_gene_id", "symbol"]],
    left_on="human_protein",
    right_on="uniprot_id",
    how="left"
)

df_merged = df_merged.rename(columns={
    "symbol": "human_gene",
    "ensembl_gene_id": "human_ensembl_gene_id"
}).drop(columns=["uniprot_id"])

df_merged = df_merged.dropna(subset=["human_ensembl_gene_id"])
df_merged.to_csv(f"data/{species_input}_inParanoid_withHGNC.tsv", sep="\t", index=False)

# === Step 5: Optional - Species UniProt → Ensembl mapping ===
map_path = f"data/{species_input}_uniprot_to_ensembl.tsv"
try:
    species_map = pd.read_csv(map_path, sep="\t", dtype=str)
    df_merged = df_merged.merge(
        species_map,
        left_on=f"{species_input}_protein",
        right_on="uniprotswissprot",
        how="left"
    ).rename(columns={"ensembl_gene_id": f"{species_input}_ensembl_gene_id"}) \
     .drop(columns=["uniprotswissprot"])
except FileNotFoundError:
    print(f"⚠️  Mapping file not found: {map_path}")

df_merged.to_csv(f"data/df_merged_with_{species_input}_ensembl.tsv", sep="\t", index=False)
