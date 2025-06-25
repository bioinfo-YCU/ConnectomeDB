import sys
import os
import pandas as pd
import warnings
import re
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)
from itables import init_notebook_mode, show
from IPython.display import display, Javascript
import itables.options as opt
from createDataTable_perSpecies import mouse_gene_pair1

mouse_gene_pair1= mouse_gene_pair1[['<span title="Double-click header of Interaction ID to ensure all values are shown">Interaction ID&nbsp;</span>','<span title="Double-click header of Ligand to ensure all values are shown">Ligand&nbsp;</span>',
       '<span title="Double-click header of Receptor to ensure all values are shown">Receptor&nbsp;</span>',
                  'LR Pair Card', 'Mouse LR Pair','<span title="HUGO Gene Nomenclature Committee (HGNC) ID. Click on the link for more details">Ligand HGNC ID&nbsp;&nbsp;</span>',
       '<span title="HUGO Gene Nomenclature Committee (HGNC) ID. Click on the link for more details">Receptor HGNC ID&nbsp;&nbsp;</span>','Ligand GOC score', 'Ligand WGA coverage',
       'Ligand % Identity', 'Ligand Target % Identity',
       'Ligand Orthology Confidence', 'Ligand Ensembl ID',
       'Receptor GOC score', 'Receptor WGA coverage', 'Receptor % Identity',
       'Receptor Target % Identity', 'Receptor Orthology Confidence','Receptor Ensembl ID']]

mouse_gene_pair1.columns = [
    "Interaction ID",
    "Ligand",
    "Receptor",
    "LR Pair Card", 
    "Mouse LR Pair",
    "Ligand HGNC ID",
    "Receptor HGNC ID",
    *mouse_gene_pair1.columns[7:]
]
mousebioM_df = pd.read_csv("data/mmusculus_ID_biomart.csv", dtype=str)
mousebioM_df = mousebioM_df.dropna(subset=["mmusculus_homolog_ensembl_gene", "ensembl_gene_id"])


def extract_link_text(html_string):
    """Extract visible text from an anchor tag <a>...</a>."""
    match = re.search(r'<a[^>]*>(.*?)</a>', html_string)
    if match:
        return match.group(1).strip()
    return None

mouse_gene_pair1['LR Pair Card'] = mouse_gene_pair1['LR Pair Card'].apply(extract_paircard_id)
# Create the mapping dictionary from mouse to human Ensembl gene ID
mouse_to_human_map = dict(zip(
    mousebioM_df["mmusculus_homolog_ensembl_gene"],
    mousebioM_df["ensembl_gene_id"]
))

# Map Ligand
mouse_gene_pair1["Human Ligand Ensembl ID"] = mouse_gene_pair1["Ligand Ensembl ID"].map(mouse_to_human_map)

# Map Receptor
mouse_gene_pair1["Human Receptor Ensembl ID"] = mouse_gene_pair1["Receptor Ensembl ID"].map(mouse_to_human_map)
def extract_hgnc_id(col):
    """Use regular expression to extract the HGNC ID after 'HGNC:'."""
    match = re.search(r'HGNC:(\d+)', col)
    if match:
        return 'HGNC:' +str(match.group(1))
    return None
    
mouse_gene_pair1['Ligand HGNC ID'] = mouse_gene_pair1['Ligand HGNC ID'].apply(extract_hgnc_id)
mouse_gene_pair1['Receptor HGNC ID'] = mouse_gene_pair1['Receptor HGNC ID'].apply(extract_hgnc_id)
df_merged =pd.read_csv("data/df_merged_with_mouse_ensembl.tsv",sep="\t")

# Step 0: Add 'orig_row' once, keep it clean
mouse_gene_pair1_indexed = mouse_gene_pair1.reset_index(drop=False).rename(columns={"index": "orig_row"})

### === LIGAND MERGE === ###
df_ligand = df_merged.add_prefix("Ligand_")
ligand_merge = mouse_gene_pair1_indexed.merge(
    df_ligand,
    left_on="Human Ligand Ensembl ID",
    right_on="Ligand_human_ensembl_gene_id",
    how="left"
)

# Ensure orig_row is single column (sometimes merge creates duplicates with suffix)
if isinstance(ligand_merge.columns, pd.MultiIndex):
    ligand_merge.columns = ligand_merge.columns.get_level_values(0)

if ligand_merge.columns.duplicated().any():
    ligand_merge = ligand_merge.loc[:, ~ligand_merge.columns.duplicated()]

def resolve_ligand_row(group):
    match = group[group["Ligand_mouse_ensembl_gene_id"] == group["Ligand Ensembl ID"]]
    if len(match) == 1:
        return match
    elif len(match) > 1:
        return match.iloc[[0]]
    else:
        return group.iloc[[0]]

ligand_final = (
    ligand_merge.groupby("orig_row", group_keys=False)
    .apply(resolve_ligand_row)
    .reset_index(drop=True)
)

### === RECEPTOR MERGE === ###
df_receptor = df_merged.add_prefix("Receptor_")
receptor_merge = ligand_final.merge(
    df_receptor,
    left_on="Human Receptor Ensembl ID",
    right_on="Receptor_human_ensembl_gene_id",
    how="left"
)

# Same cleanup for receptor_merge
if isinstance(receptor_merge.columns, pd.MultiIndex):
    receptor_merge.columns = receptor_merge.columns.get_level_values(0)

if receptor_merge.columns.duplicated().any():
    receptor_merge = receptor_merge.loc[:, ~receptor_merge.columns.duplicated()]

def resolve_receptor_row(group):
    match = group[group["Receptor_mouse_ensembl_gene_id"] == group["Receptor Ensembl ID"]]
    if len(match) == 1:
        return match
    elif len(match) > 1:
        return match.iloc[[0]]
    else:
        return group.iloc[[0]]

final_result = (
    receptor_merge.groupby("orig_row", group_keys=False)
    .apply(resolve_receptor_row)
    .reset_index(drop=True)
    .drop(columns=["orig_row"])
)

assert len(final_result) == len(mouse_gene_pair1), f"Row mismatch: {len(final_result)} != {len(mouse_gene_pair1)}"
final_result
