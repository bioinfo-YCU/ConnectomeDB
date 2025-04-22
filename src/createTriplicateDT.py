## Function to prepare functional annotation datatable
import sys, os
from itables import init_notebook_mode
import pandas as pd
from itables import show
from itables import options
from IPython.display import HTML, display
import numpy as np
from bs4 import BeautifulSoup
from createDataTable import gene_pair, gene_pair0
import warnings
import fetchGSheet 

mapping_ID = dict(zip(gene_pair0['Human LR Pair'], gene_pair0['Interaction ID']))
gene_pair_PMID = fetchGSheet.gene_pair.dropna(axis=1, how='all')
gene_pair_PMID= gene_pair_PMID[["LR pair", "PMID", "original source"]]
# Replace values in the column based on the mapping
gene_pair_PMID["Interaction ID"] = gene_pair_PMID['LR pair'].replace(mapping_ID)
df_pub = pd.read_csv("data/pubmed_results.csv", usecols=[0,1,3,4,5])
gene_pair_PMID["PMID"] = gene_pair_PMID["PMID"].astype(str)
df_pub["PMID"] = df_pub["PMID"].astype(str)
gene_pair_trip = pd.merge(gene_pair_PMID, df_pub, how='left', on='PMID')
gene_pair_trip["Year"] = pd.to_numeric(gene_pair_trip["Year"], errors="coerce").astype("Int64")
gene_pair_trip = gene_pair_trip.merge(gene_pair, how='left', left_on='Interaction ID', right_on=gene_pair.columns[0])
gene_pair_trip = gene_pair_trip.drop(columns=["Interaction ID", "LR pair", gene_pair.columns[8], gene_pair.columns[11]])
gene_pair_trip = gene_pair_trip.drop_duplicates()
gene_pair_trip = gene_pair_trip.reset_index(drop=True)  
def generate_perplexity_link_pathway(row): ## EDIT TOMORROW [note on Apr 22, 18:34]
    if pd.isna(row["Related Pathway"]) or row["Related Pathway"] == "unknown":
        query = f"What-biological-pathway-is-the-{row['Human LR Pair']}-associated-with"
    else:
        query = f"What-is-the-role-of-the-ligand-and-receptor-pair-{row['Human LR Pair']}-in-{row['Related Pathway']}"
    
    return (
        f'<a href="https://www.perplexity.ai/search?q={query}" target="_blank">'
        f'<img src="https://img.icons8.com/?size=30&id=0NbBuNOxUwps&format=png&color=000000" alt="Perplexity AI" /></a>'
    )

gene_pair_pathway["Perplexity"] = gene_pair_pathway.apply(generate_perplexity_link_pathway, axis=1)

# Create the links to the HTML cards
gene_pair_pathway["Human LR Pair"] = [
    f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/{lrPairOrig}.html">{lrPair}</a>'
    for lrPairOrig, lrPair in zip(gene_pair_pathway["Human LR Pair"], gene_pair_pathway["Human LR Pair"])
]