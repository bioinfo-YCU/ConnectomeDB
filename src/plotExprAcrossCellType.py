## Function to create horizontal bar plots of each gene in Human Taxon --expression log(x+1) transformed with cell types as y-axis

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
input_file="data/ExpressionGenes.txt" #"data/connectome_j.tsv"
# Get all unique genes
ligand_list = gene_pair0["Ligand"].tolist()
receptor_list = gene_pair0["Receptor"].tolist()
unique_genes = list(set(ligand_list + receptor_list))  # Combine and remove duplicates

connectomeDB = pd.read_table(input_file, sep="\t")
# All Taxon for now
#connectomeDB = connectomeDB[connectomeDB["Taxon"]== "Human"]
if "Taxon" in connectomeDB.columns:
    connectomeDB = connectomeDB.drop(columns=["Localization", "Taxon"] + [col for col in connectomeDB.columns if col.startswith("F5_")])

intersection = pd.Series(list(set(connectomeDB['ApprovedSymbol']).intersection(unique_genes)))

connectomeDB = connectomeDB[connectomeDB["ApprovedSymbol"].isin(intersection)]
    
# log(x+1) transform
connectomeDB.iloc[:, 1:] = np.log1p(connectomeDB.iloc[:, 1:])
# Reshape 
connectomeDB_long = connectomeDB.melt(id_vars=["ApprovedSymbol"], 
                                      var_name="cellTypes", value_name="expr_val")
cellCat = pd.read_csv("data/cell_categories.csv")
connectomeDB_long = connectomeDB_long.merge(cellCat, how='left', left_on='cellTypes', right_on='cellType')
connectomeDB_long = connectomeDB_long.drop(columns=["cellType"])
connectomeDB_long= connectomeDB_long.drop_duplicates()

intersection = pd.Series(list(set(connectomeDB_long['cellTypes']).intersection(set(cellCat['cellType']))))


diff_df = pd.Series(list(set(connectomeDB_long['cellTypes']).difference(set(cellCat['cellType']))))


def plot_gene_expression(df):
    # Define the colors for each cell category
    colors = {
        "missing": "#B0B0B0",  # Neutral gray
        "other": "#D4A76A",  # Warm gold
        "mesenchymal": "#377EB8",  # Vibrant blue
        "epithelial": "#E41A1C",  # Bold red
        "hematopoietic": "#4DAF4A",  # Fresh green
        "endothelial": "#984EA3",  # Deep purple
        "nervous system": "#FF7F00",  # Bright orange
    }

    # Define sorting order for cell categories
    category_order = {cat: i for i, cat in enumerate(colors.keys())}

    for gene, sub_df in df.groupby("ApprovedSymbol"):
        # Sort by category first, then by expression value (highest first)
        sub_df = sub_df.copy()
        sub_df["category_order"] = sub_df["cellCategory"].map(category_order).fillna(len(category_order))
        sub_df = sub_df.sort_values(["category_order", "expr_val"], ascending=[True, False])

        num_bars = len(sub_df)

        # Plotly Figure setup
        fig = go.Figure()

        # Loop through each category and create a trace for it
        for category, color in colors.items():
            # Filter data for the current category
            category_data = sub_df[sub_df["cellCategory"] == category]

            # Add the trace for the current category
            fig.add_trace(go.Bar(
                y=category_data["cellTypes"],  # Categories for y-axis
                x=category_data["expr_val"],  # Expression values for x-axis
                orientation='h',  # Horizontal bars
                marker=dict(color=color),
                hovertemplate=
                    '<b>%{y}</b><br>' +  # Cell type (y-axis value)
                    'Expression Value: %{x}',  # Expression value (x-axis value)
                    #'Category: %{text}',  # Custom text (cell category)
                #text=category_data["cellCategory"],  # Pass the cell category as custom text
                name=category,  # Use the category name for the legend
                showlegend=True,  # Ensure the legend is shown for this trace
            ))

        # Update layout settings
        fig.update_layout(
            autosize=True,
            width=450,  # Adjust width to your desired size
            #title_font=dict(size=10), 
            title="",
            xaxis_title="log(x+1) Expr value",
            #yaxis_title="Cell Types",
            yaxis=dict(
                tickmode='array',
                tickvals=np.arange(num_bars),
                ticktext=sub_df["cellTypes"],
                tickangle=0,  # Avoid overlapping labels by setting the angle to 0
                tickfont=dict(size=6),  # Set font size for the labels
            ),
            showlegend=True,
            legend_title="Cell Category",
            legend=dict(
                orientation="v",  # Vertical legend
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.05,  # Position the legend outside of the plot area
                font=dict(size=10)
            ),
            margin=dict(t=0, b=50, l=150, r=50),
            height=min(900, max(450, num_bars * 30)),  # Adjust plot height
            plot_bgcolor='rgba(0,0,0,0)',  # Transparent background
            paper_bgcolor='rgba(0,0,0,0)',  # Transparent paper background
        )

        # Save to HTML file
        fig.write_html(f"data/gene_expr_plots/{gene}.html")


plot_gene_expression(connectomeDB_long)
