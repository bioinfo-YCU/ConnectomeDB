## Function to prepare functional annotation datatable
import sys, os
from itables import init_notebook_mode
import pandas as pd
from itables import show
from itables import options
from IPython.display import HTML, display
import numpy as np
from bs4 import BeautifulSoup
from createDataTable import gene_pair0, gene_pair, top_pathway_df
import warnings

gene_pair_annot = gene_pair0[["Interaction ID", "Human LR Pair", "Cancer-related", "Ligand symbol and aliases",  "Receptor symbol and aliases"]]
df= pd.read_csv("data/disease_annotations_per_pair.csv")
df_cat=pd.read_csv("data/disease_categories.csv")
mapping = dict(zip(df_cat['Disease Name'], df_cat['Category']))
# Replace values in the column based on the mapping
df["Disease Type"] = df['disease'].replace(mapping)
gene_pair_annot = gene_pair_annot.merge(df, how='left', left_on='Human LR Pair', right_on='interaction')
gene_pair_annot = gene_pair_annot.drop(columns=["interaction"])
# PROGENy Pathway retrieved via LIANA+
df= pd.read_csv("data/pathway_annotations_per_pair.csv") 
gene_pair_annot = gene_pair_annot.merge(df, how='left', left_on='Human LR Pair', right_on='interaction')
gene_pair_annot = gene_pair_annot.drop(columns=["interaction", "weight"])

gene_pair_annot = gene_pair_annot.rename(columns={
                                     "disease": "Disease", 
                                     "source": "PROGENy Pathway"}
                            )
# Bring in KEGG Pathways from AU side
gene_pair_annot = gene_pair_annot.merge(top_pathway_df, how='left', left_on='Human LR Pair', right_on='LR Pair')

# reorder
gene_pair_annot = gene_pair_annot[["Interaction ID", "Human LR Pair", "Disease", "Disease Type", "Cancer-related", "KEGG Pathway ID", "KEGG Pathway", "KEGG relationship", "PROGENy Pathway", "Ligand symbol and aliases",  "Receptor symbol and aliases"]]


gene_pair_annot["Disease"] = gene_pair_annot["Disease"].apply(
    lambda x: "unknown" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x)
gene_pair_annot["Disease Type"] = gene_pair_annot["Disease Type"].apply(
    lambda x: "unknown" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x)
gene_pair_annot["PROGENy Pathway"] = gene_pair_annot["PROGENy Pathway"].apply(
    lambda x: "unknown" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x)
gene_pair_annot["KEGG Pathway"] = gene_pair_annot["KEGG Pathway"].apply(
    lambda x: "unknown" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x)

gene_pair_annot["KEGG Pathway ID"] = gene_pair_annot["KEGG Pathway ID"].apply(
    lambda x: "unknown" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x)

gene_pair_annot["KEGG relationship"] = gene_pair_annot["KEGG relationship"].apply(
    lambda x: "unknown" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x)

gene_pair_annot = gene_pair_annot.reset_index(drop=True)

# Separate Disease and Pathway and then rm duplicates
gene_pair_disease = gene_pair_annot[["Interaction ID", "Human LR Pair", "Disease", "Disease Type", "Cancer-related", "Ligand symbol and aliases",  "Receptor symbol and aliases"]]
gene_pair_disease = gene_pair_disease.drop_duplicates()
gene_pair_disease=gene_pair_disease.reset_index(drop=True)  

def generate_perplexity_link(row):
    if pd.isna(row["Disease"]) or row["Disease"] == "unknown":
        query = f"What-disease-is-the-{row['Human LR Pair']}-associated-with"
    else:
        query = f"What-is-the-role-of-the-ligand-and-receptor-pair-{row['Human LR Pair']}-in-{row['Disease']}"
    
    return (
        f'<a href="https://www.perplexity.ai/search?q={query}" target="_blank">'
        f'<img src="https://img.icons8.com/?size=30&id=0NbBuNOxUwps&format=png&color=000000" alt="Perplexity AI" /></a>'
    )

gene_pair_disease["Perplexity"] = gene_pair_disease.apply(generate_perplexity_link, axis=1)

# Create the links to the HTML cards
gene_pair_disease["Human LR Pair"] = [
    f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/{lrPairOrig}.html">{lrPair}</a>'
    for lrPairOrig, lrPair in zip(gene_pair_disease["Human LR Pair"], gene_pair_disease["Human LR Pair"])
]
gene_pair_pathway = gene_pair_annot[["Interaction ID", "Human LR Pair",  "KEGG Pathway ID", "KEGG Pathway", "KEGG relationship", "PROGENy Pathway", "Ligand symbol and aliases",  "Receptor symbol and aliases"]]
gene_pair_pathway = gene_pair_pathway.drop_duplicates()
gene_pair_pathway=gene_pair_pathway.reset_index(drop=True)  
def generate_perplexity_link_pathway(row):
    if pd.isna(row["KEGG Pathway"]) or row["KEGG Pathway"] == "unknown":
        if pd.isna(row["PROGENy Pathway"]) or row["PROGENy Pathway"] == "unknown":
            query = f"What-biological-pathway-is-the-{row['Human LR Pair']}-associated-with"
        else:
            query = f"What-is-the-role-of-the-ligand-and-receptor-pair-{row['Human LR Pair']}-in-{row['PROGENy Pathway']}"
    else:
        query = f"What-is-the-role-of-the-ligand-and-receptor-pair-{row['Human LR Pair']}-in-{row['KEGG Pathway']}"
    
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