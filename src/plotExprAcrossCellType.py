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
expr_data = connectomeDB.drop(columns=["ApprovedSymbol"])
expr_data_scaled = expr_data.div(expr_data.max(axis=1).replace(0, np.nan), axis=0).fillna(0)

# Reattach gene symbols
connectomeDB_scaled = pd.concat([connectomeDB["ApprovedSymbol"], expr_data_scaled], axis=1)

# Reshape to long format
connectomeDB_long = connectomeDB_scaled.melt(id_vars="ApprovedSymbol", 
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

    # Define sorting order
    category_order = {cat: i for i, cat in enumerate(colors.keys())}

    for gene, sub_df in df.groupby("ApprovedSymbol"):
        sub_df = sub_df.copy()
        sub_df["category_order"] = sub_df["cellCategory"].map(category_order).fillna(len(category_order))

        # Sort bars top to bottom
        sub_df = sub_df.sort_values(["category_order", "expr_val"], ascending=[True, False])
        y_order = sub_df["cellTypes"].tolist()

        # Apply y-axis order explicitly
        sub_df = sub_df.set_index("cellTypes").loc[y_order].reset_index()

        # Get categories in top-down order, then reverse it for the legend to match
        ordered_categories = sub_df["cellCategory"].dropna().unique().tolist()

        fig = go.Figure()

        # Plot bars by category (reversed for legend alignment)
        for category in ordered_categories:
            category_data = sub_df[sub_df["cellCategory"] == category]
            color = colors.get(category, "#B0B0B0")

            fig.add_trace(go.Bar(
                y=category_data["cellTypes"],
                x=category_data["expr_val"],
                orientation='h',
                marker=dict(color=color),
                hovertemplate='<b>%{y}</b><br>Expression Value: %{x}',
                name=category,
                showlegend=True,
            ))

        # Layout
        num_bars = len(sub_df)
        fig.update_layout(
            autosize=True,
            width=450,
            title="",
            xaxis=dict(
                title="Expr value (TPM) scaled per gene",
                # type="log",  # Set x-axis to log scale
                # tickmode='array',  # Manually control tick marks
                #tickvals=[1, 10, 100, 1000, 10000, 100000, 1000000],  # Define log scale ticks
                # ticktext=['1', '10', '100', '1000', '10000', '1000000'],  # Corresponding tick labels
            ),
            yaxis=dict(
                tickmode='array',
                tickvals=np.arange(num_bars),
                ticktext=sub_df["cellTypes"],
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
            height=min(900, max(450, num_bars * 30)),
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)',
        )

        fig.write_html(f"{output_dir}/{gene}.html")

plot_gene_expression(connectomeDB_long)
