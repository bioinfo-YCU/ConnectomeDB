## Function to create horizontal bar plots of each gene in Human Taxon --expression log(x+1) transformed with cell types as y-axis

import requests
import pandas as pd
import numpy as np
import sys
import os

sys.path.append(os.path.abspath("src"))  # Add src directory to path
from createDataTable import gene_pair0

# Get all unique genes
ligand_list = gene_pair0["Ligand"].tolist()
receptor_list = gene_pair0["Receptor"].tolist()
unique_genes = list(set(ligand_list + receptor_list))  # Combine and remove duplicates

connectomeDB = pd.read_table("data/connectome_j.tsv", sep="\t")
connectomeDB = connectomeDB[connectomeDB["Taxon"]== "Human"]
connectomeDB = connectomeDB.drop(columns=["Localization", "Taxon"] + [col for col in connectomeDB.columns if col.startswith("F5_")])
# log(x+1) transform
connectomeDB.iloc[:, 1:] = np.log1p(connectomeDB.iloc[:, 1:])
# Reshape 
connectomeDB_long = connectomeDB.melt(id_vars=["ApprovedSymbol"], 
                                      var_name="cellTypes", value_name="expr_val")
cellCat = pd.read_csv("data/cell_categories.csv")
connectomeDB_long = connectomeDB_long.merge(cellCat, how='left', left_on='cellTypes', right_on='cellType')
connectomeDB_long = connectomeDB_long.drop(columns=["cellType"])

intersection = pd.Series(list(set(connectomeDB_long['cellTypes']).intersection(set(cellCat['cellType']))))
intersection

diff_df = pd.Series(list(set(connectomeDB_long['cellTypes']).difference(set(cellCat['cellType']))))
diff_df

import matplotlib.pyplot as plt
import numpy as np

def plot_gene_expression(df):
    colors = {
        "mesenchymal": "#1f77b4",  # Blue  
        "epithelial": "#ff7f0e",   # Orange  
        "hematopoietic": "#2ca02c",  # Green  
        "endothelial": "#d62728",  # Red  
        "nervous system": "#9467bd",  # Purple  
        "other": "#8c564b",  # Brown  
        "missing": "#7f7f7f",  # Gray  
        np.nan: "#e377c2"  # Pink for NaN values
    }

    # Define sorting order for cell categories
    category_order = {cat: i for i, cat in enumerate(colors.keys())}

    for gene, sub_df in df.groupby("ApprovedSymbol"):
        # Sort by category first, then by expression value (highest first)
        sub_df = sub_df.copy()
        sub_df["category_order"] = sub_df["cellCategory"].map(category_order).fillna(len(category_order))
        sub_df = sub_df.sort_values(["category_order", "expr_val"], ascending=[True, False])

        fig, ax = plt.subplots(figsize=(8, min(10, max(4, len(sub_df) * 0.3))))  # Limit vertical size

        # Assign colors
        bar_colors = sub_df["cellCategory"].map(colors).fillna("gray")

        # Plot
        ax.barh(sub_df["cellTypes"], sub_df["expr_val"], color=bar_colors)

        ax.set_title(gene, fontsize=14)
        ax.set_xlabel("Expression Value (log(x+1) transformed", fontsize=12)
        ax.set_ylabel("")
        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=6)

        # Legend: Only show categories that exist in the plot
        unique_categories = sub_df["cellCategory"].dropna().unique()
        legend_handles = [plt.Rectangle((0, 0), 1, 1, color=colors[cat]) for cat in unique_categories]
        ax.legend(legend_handles, unique_categories, title="Cell Category", fontsize=9, title_fontsize=11, loc="upper right")

        plt.tight_layout()
        plt.savefig(f"data/gene_expr_plots/{gene}.png", dpi=150, bbox_inches="tight")
        plt.close()

plot_gene_expression(connectomeDB_long)
