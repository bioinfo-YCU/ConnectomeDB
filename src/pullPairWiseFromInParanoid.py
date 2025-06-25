# import inParanoid pairwise Data
# Column Name	Description
# bitscore	BLAST score of the orthologous pair
# score	InParalog score (how well the protein aligns with the cluster seed)
# seed_score	Score of the seed protein (usually 1.0)

import pandas as pd
import requests
from io import StringIO
from itertools import product

# Download InParanoid prot table
url = "https://inparanoidb.sbc.su.se/download/sqltable/9606&10090&prot"
r = requests.get(url)
r.raise_for_status()

df = pd.read_csv(StringIO(r.text.strip()), sep="\t", header=None)
df.columns = ["cluster_id", "bitscore", "source_file", "inparalog_score", "protein_id", "seed_score"]

# Tag each row by species
def infer_species(src):
    if "9606" in src:
        return "human"
    elif "10090" in src:
        return "mouse"
    return "unknown"

df["species"] = df["source_file"].apply(infer_species)

# Expand ortholog pairs within each cluster
records = []
for cid, grp in df.groupby("cluster_id"):
    humans = grp[grp["species"] == "human"]
    mice = grp[grp["species"] == "mouse"]
    for h, m in product(humans.itertuples(index=False), mice.itertuples(index=False)):
        records.append({
            "cluster_id": cid,
            "human_protein": h.protein_id,
            "human_inparalog_score": h.inparalog_score,
            "human_seed_score": h.seed_score,
            "mouse_protein": m.protein_id,
            "mouse_inparalog_score": m.inparalog_score,
            "mouse_seed_score": m.seed_score,
            "bitscore": (h.bitscore + m.bitscore) / 2  # average for now
        })

df_orthologs = pd.DataFrame(records)

df_orthologs.to_csv("data/inParanoid_mmusculus.csv")

hgnc_df = pd.read_csv("data/HGNC_gene_info_full.tsv", sep="\t", dtype=str)
hgnc_df = hgnc_df.dropna(subset=["uniprot_ids", "ensembl_gene_id"])
# Split uniprot_ids on comma and explode
hgnc_exploded = hgnc_df.assign(uniprot_id=hgnc_df["uniprot_ids"].str.split(",")).explode("uniprot_id")
hgnc_exploded["uniprot_id"] = hgnc_exploded["uniprot_id"].str.strip()

uniprot_to_ensembl = hgnc_exploded.set_index("uniprot_id")["ensembl_gene_id"].to_dict()

# Left join on human_protein
df_merged = df_orthologs.merge(
    hgnc_exploded[["uniprot_id", "ensembl_gene_id", "symbol"]],
    left_on="human_protein",
    right_on="uniprot_id",
    how="left"
)

# Optionally rename
df_merged = df_merged.rename(columns={
    "symbol": "human_gene",
    "ensembl_gene_id": "human_ensembl_gene_id"
}).drop(columns=["uniprot_id"])
df_merged = df_merged.dropna(subset=["human_ensembl_gene_id"])
df_merged.to_csv("data/mmusculus_inParanoid_uniProt_withHGNCAnn.tsv", sep="\t", index=False)

# Step 1: Load the mapping file
mouse_map = pd.read_csv("data/mouse_uniprot_to_ensembl.tsv", sep="\t", dtype=str)

# Step 2: Merge with df_merged on UniProt ID
# Assuming your UniProt column in df_merged is named 'mouse_protein'
df_merged = df_merged.merge(
    mouse_map,
    left_on="mouse_protein",
    right_on="uniprotswissprot",
    how="left"
)

# Step 3: Rename the new column for clarity (optional)
df_merged = df_merged.rename(columns={"ensembl_gene_id": "mouse_ensembl_gene_id"})

# Step 4: Drop the helper merge column if not needed
df_merged = df_merged.drop(columns=["uniprotswissprot"])

df_merged.to_csv("data/df_merged_with_mouse_ensembl.tsv", sep="\t", index=False)
