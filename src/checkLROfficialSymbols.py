import requests
import pandas as pd
import sys
import os

sys.path.append(os.path.abspath("src"))  # Add src directory to path
from createDataTable import gene_pair0

species = "Pig"
species_lower = species.lower()
id = "ENSEMBL"
species_name = {
        "Mouse": "mmusculus",
        "Rat": "rnorvegicus",
        "Zebrafish":"drerio" ,
        "Chimp":"ptroglodytes",
        "Chicken":"ggallus",
        "Pig":"sscrofa",
        "Cow":"btaurus",
        "Dog":"clfamiliaris",
        "Horse":"ecaballus",
        "Sheep":"oarambouillet",
        "Marmoset": "cjacchus" ,
        "Macaque": "mmulatta",
        "Frog": "xtropicalis"
    }.get(species, "Unknown species")

ref = pd.read_csv(f"data/hsapiens_ID_biomart_{species_name}_centric.csv")
ref = ref[["ensembl_gene_id", "external_gene_name", "external_synonym"]]

mapping_orth_symbol = dict(zip(ref['ensembl_gene_id'], ref['external_gene_name']))

gene_pair_orth['Ligand Official Symbol'] = 'NA'
gene_pair_orth['Ligand Official Symbol'] =  gene_pair_orth.apply(
    lambda row: mapping_orth_symbol.get(row[f'{id} ligand'], row['Ligand Official Symbol'])
    if pd.notna(row[f'{id} ligand']) else row['Ligand Official Symbol'],
    axis=1
)
gene_pair_orth['Receptor Official Symbol'] = 'NA'
gene_pair_orth['Receptor Official Symbol'] =  gene_pair_orth.apply(
    lambda row: mapping_orth_symbol.get(row[f'{id} receptor'], row['Receptor Official Symbol'])
    if pd.notna(row[f'{id} receptor']) else row['Receptor Official Symbol'],
    axis=1
)
gene_pair_orth["same_as_off_lig"] = (
    gene_pair_orth[f"{species}_ligand"] == gene_pair_orth["Ligand Official Symbol"]
)

gene_pair_orth["same_as_off_rec"] = (
    gene_pair_orth[f"{species}_receptor"] == gene_pair_orth["Receptor Official Symbol"]
)

gene_pair_orth.to_csv(f"data/{species}_{id}_ID_check.csv")