import sys
import os
import pandas as pd
import warnings
import re
from itables import init_notebook_mode, show
from IPython.display import display, Javascript
import itables.options as opt
import createDataTable_perSpecies

warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

# === Species Parameter === #
species = "horse"  # Change to "zebrafish", "sheep", etc.
species_file_prefix = {
    "mouse": "mmusculus",
    "rat": "rnorvegicus",
    "zebrafish": "drerio",
    "chimpanzee": "ptroglodytes",
    "chicken": "ggallus",
    "pig": "sscrofa",
    "cow": "btaurus",
    "dog": "clfamiliaris",
    "horse": "ecaballus",
    "marmoset": "cjacchus",
    "macaque": "mmulatta",
    "sheep": "oarambouillet"
}[species]

# === Load gene pair === #
gene_pair_var = f"{species}_gene_pair1"
gene_pair_df = getattr(createDataTable_perSpecies, gene_pair_var)
ligand_ens_id = [col for col in gene_pair_df.columns if "Ligand Ensembl ID" in col][0]
receptor_ens_id = [col for col in gene_pair_df.columns if "Receptor Ensembl ID" in col][0]


gene_pair_df = gene_pair_df[[
    '<span title="Double-click header of Interaction ID to ensure all values are shown">Interaction ID&nbsp;</span>',
    '<span title="Double-click header of Ligand to ensure all values are shown">Ligand&nbsp;</span>',
    '<span title="Double-click header of Receptor to ensure all values are shown">Receptor&nbsp;</span>',
    'LR Pair Card', f'{species.capitalize()} LR Pair',
    '<span title="HUGO Gene Nomenclature Committee (HGNC) ID. Click on the link for more details">Ligand HGNC ID&nbsp;&nbsp;</span>',
    '<span title="HUGO Gene Nomenclature Committee (HGNC) ID. Click on the link for more details">Receptor HGNC ID&nbsp;&nbsp;</span>',
    'Ligand GOC score', 'Ligand WGA coverage',
    'Ligand % Identity', 'Ligand Target % Identity',
    'Ligand Orthology Confidence', ligand_ens_id,
    'Receptor GOC score', 'Receptor WGA coverage', 'Receptor % Identity',
    'Receptor Target % Identity', 'Receptor Orthology Confidence', receptor_ens_id]]

# Rename columns
rename_dict = dict(zip(gene_pair_df.columns[:7], [
    "Interaction ID", "Ligand", "Receptor", "LR Pair Card",
    f"{species.capitalize()} LR Pair", "Ligand HGNC ID", "Receptor HGNC ID"
]))
gene_pair_df.rename(columns=rename_dict, inplace=True)

gene_pair_df = gene_pair_df.rename(columns={
    ligand_ens_id: "Ligand Ensembl ID",
    receptor_ens_id: "Receptor Ensembl ID"
})
    

# Load ortholog mapping
biomart_df = pd.read_csv(f"data/{species_file_prefix}_ID_biomart.csv", dtype=str)
biomart_df = biomart_df.dropna(subset=[f"{species_file_prefix}_homolog_ensembl_gene", "ensembl_gene_id"])

# Extract ID from anchor tags
def extract_link_text(html_string):
    match = re.search(r'<a[^>]*>(.*?)</a>', html_string)
    return match.group(1).strip() if match else None

def extract_hgnc_id(col):
    match = re.search(r'HGNC:(\d+)', col)
    return 'HGNC:' + str(match.group(1)) if match else None
def extract_paircard_id(col):
    """Use regular expression to extract the HGNC ID after 'cards/'."""
    match = re.search(r'cards/([^/]+)\.html', col)
    if match:
        return str(match.group(1))
    return None
    
# Process columns
gene_pair_df['LR Pair Card'] = gene_pair_df['LR Pair Card'].apply(extract_paircard_id)
gene_pair_df['Ligand HGNC ID'] = gene_pair_df['Ligand HGNC ID'].apply(extract_hgnc_id)
gene_pair_df['Receptor HGNC ID'] = gene_pair_df['Receptor HGNC ID'].apply(extract_hgnc_id)

# Mapping
species_to_human_map = dict(zip(
    biomart_df[f"{species_file_prefix}_homolog_ensembl_gene"],
    biomart_df["ensembl_gene_id"]
))

gene_pair_df["Human Ligand Ensembl ID"] = gene_pair_df["Ligand Ensembl ID"].map(species_to_human_map)
gene_pair_df["Human Receptor Ensembl ID"] = gene_pair_df["Receptor Ensembl ID"].map(species_to_human_map)

# Load df_merged
merged_df = pd.read_csv(f"data/df_merged_with_{species}_ensembl.tsv", sep="\t")

# Index for merge
gene_pair_indexed = gene_pair_df.reset_index(drop=False).rename(columns={"index": "orig_row"})

# LIGAND MERGE
df_ligand = merged_df.add_prefix("Ligand_")
ligand_merge = gene_pair_indexed.merge(
    df_ligand,
    left_on="Human Ligand Ensembl ID",
    right_on="Ligand_human_ensembl_gene_id",
    how="left"
)
if isinstance(ligand_merge.columns, pd.MultiIndex):
    ligand_merge.columns = ligand_merge.columns.get_level_values(0)
if ligand_merge.columns.duplicated().any():
    ligand_merge = ligand_merge.loc[:, ~ligand_merge.columns.duplicated()]

def resolve_ligand_row(group):
    match = group[group[f"Ligand_{species}_ensembl_gene_id"] == group["Ligand Ensembl ID"]]
    return match.iloc[[0]] if len(match) else group.iloc[[0]]

ligand_final = ligand_merge.groupby("orig_row", group_keys=False).apply(resolve_ligand_row).reset_index(drop=True)

# RECEPTOR MERGE
df_receptor = merged_df.add_prefix("Receptor_")
receptor_merge = ligand_final.merge(
    df_receptor,
    left_on="Human Receptor Ensembl ID",
    right_on="Receptor_human_ensembl_gene_id",
    how="left"
)
if isinstance(receptor_merge.columns, pd.MultiIndex):
    receptor_merge.columns = receptor_merge.columns.get_level_values(0)
if receptor_merge.columns.duplicated().any():
    receptor_merge = receptor_merge.loc[:, ~receptor_merge.columns.duplicated()]

def resolve_receptor_row(group):
    match = group[group[f"Receptor_{species}_ensembl_gene_id"] == group["Receptor Ensembl ID"]]
    return match.iloc[[0]] if len(match) else group.iloc[[0]]

final_result = receptor_merge.groupby("orig_row", group_keys=False).apply(resolve_receptor_row).reset_index(drop=True).drop(columns=["orig_row"])

assert len(final_result) == len(gene_pair_df), f"Row mismatch: {len(final_result)} != {len(gene_pair_df)}"

final_result.to_csv(f"data/human_{species}_merged_ensemblBiomaRt_inParanoid.csv")

