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
import string

def make_ids_unique(series):
    return [
        f"{id_val}{letter}"
        #if count > 1 else id_val -- UNCOMMENT SO UNIQUE ONES ARE ALSO WITH APPENDED LETTER ["A"]
        for id_val, count, letter in zip(
            series,
            series.groupby(series).transform('count'),
            series.groupby(series).cumcount().map(lambda i: string.ascii_uppercase[i] if i < 26 else f"_{i}")
        )
    ]


mapping_ID = dict(zip(gene_pair0['Human LR Pair'], gene_pair0['Interaction ID']))
gene_pair_PMID = fetchGSheet.gene_pair.dropna(axis=1, how='all')
gene_pair_PMID= gene_pair_PMID[["LR pair", "PMID", "original source"]]
# Mapping for replacements
mapping = dict(zip(fetchGSheet.src_info['original source'], fetchGSheet.src_info['shortname']))
# Replace values in the column based on the mapping
gene_pair_PMID['Database Source'] = gene_pair_PMID['original source'].replace(mapping)

# Replace values in the column based on the mapping
gene_pair_PMID["Interaction ID"] = gene_pair_PMID['LR pair'].replace(mapping_ID)
df_pub = pd.read_csv("data/pubmed_results.csv", usecols=[0,1,3,4,5])
gene_pair_PMID["PMID"] = gene_pair_PMID["PMID"].astype(str)
df_pub["PMID"] = df_pub["PMID"].astype(str)
gene_pair_trip = pd.merge(gene_pair_PMID, df_pub, how='left', on='PMID')
gene_pair_trip["Year"] = pd.to_numeric(gene_pair_trip["Year"], errors="coerce").astype("Int64")
gene_pair_trip = gene_pair_trip.merge(gene_pair, how='left', left_on='Interaction ID', right_on=gene_pair.columns[0])
gene_pair_trip = gene_pair_trip.drop(columns=["Interaction ID", gene_pair.columns[2], gene_pair.columns[5],gene_pair.columns[6]])
gene_pair_trip = gene_pair_trip.drop_duplicates()
gene_pair_trip = gene_pair_trip.reset_index(drop=True)  
def generate_perplexity_link_pmid(row): 
    query = f"What-is-the-biological-relevance-of-the-ligand-and-receptor-pair-{row['LR pair']}-based-on-Pubmed-ID-{row['PMID']}"
    return (
        f'<a href="https://www.perplexity.ai/search?q={query}" target="_blank">'
        f'<img src="https://img.icons8.com/?size=30&id=0NbBuNOxUwps&format=png&color=000000" alt="Perplexity AI" /></a>'
    )

gene_pair_trip["Perplexity"] = gene_pair_trip.apply(generate_perplexity_link_pmid, axis=1)
# Add
gene_pair_trip = gene_pair_trip.drop(columns=["LR pair", "original source"])
first_columns=[gene_pair_trip.columns[6], gene_pair_trip.columns[7], 'Perplexity', 'Database Source', 'PMID']


# Make ID unique
gene_pair_trip = gene_pair_trip.sort_values(by='Year', ascending=True)
gene_pair_trip[gene_pair_trip.columns[6]] = make_ids_unique(gene_pair_trip[gene_pair_trip.columns[6]])
gene_pair_trip = gene_pair_trip.sort_values(by='Year', ascending=False)
gene_pair_trip = gene_pair_trip[first_columns + [col for col in gene_pair_trip.columns if col not in first_columns]]
gene_pair_trip = gene_pair_trip.reset_index(drop=True)  