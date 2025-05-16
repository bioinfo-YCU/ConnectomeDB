## Function to create Expression plots per gene to be placed inside the PDF/HTML card
import requests
import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import plotly.graph_objects as go

sys.path.append(os.path.abspath("src"))  # Add src directory to path
from createDataTable import gene_pair0

# Input file
input_file = "data/ExpressionGenes.txt"
ligand_list = gene_pair0["Ligand"].tolist()
receptor_list = gene_pair0["Receptor"].tolist()
unique_genes = list(set(ligand_list + receptor_list))
output_dir = "data/gene_expr_plots"
os.makedirs(output_dir, exist_ok=True)

connectomeDB = pd.read_table(input_file, sep="\t")

# Drop metadata and unwanted columns
if "Taxon" in connectomeDB.columns:
    connectomeDB = connectomeDB.drop(columns=["Localization", "Taxon"] + 
                                     [col for col in connectomeDB.columns if col.startswith("F5_")])

# Filter to genes of interest
intersection = set(connectomeDB['ApprovedSymbol']).intersection(unique_genes)
connectomeDB = connectomeDB[connectomeDB["ApprovedSymbol"].isin(intersection)]

# Isolate expression data and normalize per gene (row)
# expr_data = connectomeDB.drop(columns=["ApprovedSymbol"])
# expr_data_scaled = expr_data.div(expr_data.max(axis=1).replace(0, np.nan), # axis=0).fillna(0)

# Reattach gene symbols
# connectomeDB_scaled = pd.concat([connectomeDB["ApprovedSymbol"], # expr_data_scaled], axis=1)

# Reshape to long format
connectomeDB_long = connectomeDB.melt(id_vars="ApprovedSymbol", 
                                              var_name="cellTypes", 
                                              value_name="expr_val")

# Merge in categories
cellCat = pd.read_csv("data/cell_categories.csv")
connectomeDB_long = connectomeDB_long.merge(cellCat, how='left', left_on='cellTypes', right_on='cellType')
connectomeDB_long = connectomeDB_long.drop(columns=["cellType"])
connectomeDB_long = connectomeDB_long.drop_duplicates()

def plot_gene_expression(df):
    # Define cell category colors
    colors = {
        "mesenchymal": "#377EB8",
        "epithelial": "#E41A1C",
        "hematopoietic": "#4DAF4A",
        "endothelial": "#984EA3",
        "nervous system": "#FF7F00",
        "other": "#D4A76A",
        "missing": "#B0B0B0"
    }

    # Sorting order by category
    category_order = {cat: i for i, cat in enumerate(colors.keys())}

    for gene, sub_df in df.groupby("ApprovedSymbol"):
        sub_df = sub_df.copy()
        sub_df["category_order"] = sub_df["cellCategory"].map(category_order).fillna(len(category_order))

        # Sort by category and expression value
        sub_df = sub_df.sort_values(["category_order", "expr_val"], ascending=[True, False])
        sub_df = sub_df.reset_index(drop=True)

        # Fixed list of all cell types for consistent row order
        all_cell_types = sub_df["cellTypes"].tolist()
        all_categories = sub_df["cellCategory"].tolist()

        fig = go.Figure()

        # For each category, create a trace with all cell types
        for cat, color in colors.items():
            y_vals = []
            x_vals = []
            hover_vals = []

            for i, row in sub_df.iterrows():
                if row["cellCategory"] == cat:
                    y_vals.append(row["cellTypes"])
                    x_vals.append(row["expr_val"])
                    hover_vals.append(row["cellCategory"])
                else:
                    y_vals.append(row["cellTypes"])
                    x_vals.append(0)  # dummy bar to preserve alignment
                    hover_vals.append(None)  # no hover

            fig.add_trace(go.Bar(
                y=y_vals,
                x=x_vals,
                orientation='h',
                marker=dict(color=color),
                name=cat,
                customdata=np.array(hover_vals).reshape(-1, 1),
                hovertemplate='<b>%{y}</b><br>Expression Value: %{x}<br>Category: %{customdata[0]}<extra></extra>',
                showlegend=True
            ))

        # Layout
        fig.update_layout(
            autosize=True,
            width=450,
            title="",
            xaxis=dict(title="Expr value (TPM)"),
            yaxis=dict(
                tickmode='array',
                tickvals=np.arange(len(all_cell_types)),
                ticktext=all_cell_types,
                tickangle=0,
                tickfont=dict(size=6),
            ),
            showlegend=True,
            legend_title="Cell Category",
            legend=dict(
                orientation="v",
                yanchor="top",
                y=0.5,
                xanchor="left",
                x=1.05,
                font=dict(size=10)
            ),
            margin=dict(t=0, b=50, l=150, r=50),
            height=min(900, max(450, len(all_cell_types) * 30)),
            barmode="stack",
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
        )

        fig.write_html(f"{output_dir}/{gene}.html")


plot_gene_expression(connectomeDB_long)
