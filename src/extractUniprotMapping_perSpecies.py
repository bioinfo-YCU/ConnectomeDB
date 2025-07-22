import pandas as pd
import requests
from io import StringIO
import os

# === USER INPUT ===
species_input = "human"  # e.g., "mouse", "zebrafish", "horse"

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
    "human":         {"taxid": "9606",  "code": "hsapiens"},
}

# === VALIDATE INPUT ===
if species_input not in species_info:
    raise ValueError(f"Species '{species_input}' not supported. Choose from: {list(species_info)}")

taxid = species_info[species_input]["taxid"]

# === UNIPROT API SETTINGS ===
base_url = "https://rest.uniprot.org/uniprotkb/stream"
query = f"reviewed:true AND taxonomy_id:{taxid}"
fields = "accession,gene_names,protein_name,organism_name"
format = "tsv"

params = {
    "query": query,
    "fields": fields,
    "format": format,
    "size": 500  # Max allowed is 500 per request; use paging if needed
}

# === API REQUEST ===
print(f"Querying UniProt for Swiss-Prot entries (taxid={taxid})...")
response = requests.get(base_url, params=params)

if response.status_code != 200:
    raise RuntimeError(f"Failed to fetch data: {response.status_code} - {response.text}")

# === PARSE AND SAVE ===
df = pd.read_csv(StringIO(response.text), sep="\t")

# === SPLIT MULTIPLE GENE NAMES INTO ROWS ===
df["Gene Names"] = df["Gene Names"].fillna("")  # Just in case
df_expanded = df.assign(**{
    "Gene Names": df["Gene Names"].str.split()
}).explode("Gene Names").reset_index(drop=True)

# Ensure output directory exists
os.makedirs("data", exist_ok=True)

# Save to file
out_file = f"data/uniprotMapping_{species_input}.csv"
df_expanded.to_csv(out_file, index=False)

print(f"Saved {len(df)} entries to {out_file}")
