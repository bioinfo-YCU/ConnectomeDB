## Function to plot scaled data into a heatmap of cell types x tissues (Tabula sapiens)
import scanpy as sc
import anndata
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import requests
import scanpy as sc
import re
import pickle
sys.path.append(os.path.abspath("src"))  # Add src directory to path
import processTabulaSapiens
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from joblib import Parallel, delayed
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

input_file = "data/grouping_for_heatmap.csv"
scaled_data_file = "data/gene_expr_map_scaled.pkl"
# Enable test mode
test_mode = True
max_genes = 3 if test_mode else None

if os.path.exists(input_file):
    # Do only if the file exists
    cell_cat = pd.read_csv(input_file)
    mapping = dict(zip(cell_cat["cell_type"], cell_cat["category"]))

else:
    print("Cell grouping for heatmap missing")

if os.path.exists(scaled_data_file):
    with open(scaled_data_file, "rb") as f:
    gene_expr_map = pickle.load(f)
else:
    print("need to pull scaled gene expression data from annData")

if test_mode:
    processTabulaSapiens.valid_gene_ids = valid_gene_ids[:max_genes]
colorscale_theme = [
    [0.0, "blue"],
    [0.5, "lightgray"],
    [1.0, "red"]
]

# Map cell_type → cell_class → color
unique_classes = sorted(set(cell_type_ont))
# Define color palette
color_palette = px.colors.qualitative.Set2

# Map: cell_type → cell_class
cell_type_to_class = dict(zip(cell_types, cell_type_ont))

# Sort cell types by cell class and name
sorted_cell_types = sorted(set(cell_types), key=lambda ct: (cell_type_to_class[ct], ct))

# Extract cell class order from sorted y-axis
sorted_classes_from_yaxis = []
seen_classes = set()
for ct in sorted_cell_types:
    cls = cell_type_to_class[ct]
    if cls not in seen_classes:
        sorted_classes_from_yaxis.append(cls)
        seen_classes.add(cls)

# Assign color per class using y-axis order
class_to_color = {
    cls: color_palette[i % len(color_palette)]
    for i, cls in enumerate(sorted_classes_from_yaxis)
}

# Assign color to cell types based on class
cell_type_to_color = {ct: class_to_color[cell_type_to_class[ct]] for ct in set(cell_types)}


n_cell = str(len(np.unique(cell_types))) 
n_cellGroup = str(len(np.unique(cell_type_ont)))
n_tissue = str(len(np.unique(tissues)))

cell_cat = pd.read_csv("data/grouping_for_heatmap.csv")
mapping = dict(zip(cell_cat["cell_type"], cell_cat["category"]))

def generate_heatmap(gene_id, gene_label, expr):
    # --- Prepare data ---
    df = pd.DataFrame({
        "expression": expr,
        "tissue": tissues,
        "cell_type": cell_types,
        "cell_class": cell_type_ont
    })

    df['cell_panel'] = df['cell_class'].replace(mapping)

    # Ordered categorical cell types
    sorted_cell_types = sorted(set(cell_types), key=lambda ct: (cell_type_to_class[ct], ct))
    df["cell_type"] = pd.Categorical(df["cell_type"], categories=sorted_cell_types, ordered=True)

    # --- Global pivot ---
    pivot_df = df.groupby(["tissue", "cell_type"]).agg(
        mean_expression=("expression", "mean"),
        cell_class=("cell_class", "first")
    ).reset_index()

    global_heatmap = pivot_df.pivot(index="cell_type", columns="tissue", values="mean_expression")

    # Global zmin/zmax
    global_zmin = global_heatmap.min().min()
    global_zmax = global_heatmap.max().max()

    # Mapping for customdata
    cell_type_class_map = df.drop_duplicates("cell_type").set_index("cell_type")["cell_class"].to_dict()

    # --- Tissue clustering ---
    clustered_tissues = global_heatmap.T.fillna(0)
    linkage_matrix = linkage(pdist(clustered_tissues, metric='euclidean'), method='average')
    clustered_tissue_order = clustered_tissues.index[leaves_list(linkage_matrix)]
    global_heatmap = global_heatmap[clustered_tissue_order]

    # --- Panels setup ---
    panel_names = sorted(df['cell_panel'].dropna().unique())
    n_panels = len(panel_names)

    # Panel size logic
    row_heights = []
    for panel in panel_names:
        n = len(df[df['cell_panel'] == panel]["cell_type"].unique())
        row_heights.append(n)
    total_rows = sum(row_heights)
    normalized_heights = [h / total_rows for h in row_heights]

    fig = make_subplots(
        rows=n_panels, cols=2,
        shared_xaxes=False,
        shared_yaxes=True,
        column_widths=[0.03, 0.96],
        vertical_spacing=0.02,
        horizontal_spacing=0.005,
        row_heights=normalized_heights,
        specs=[[{"type": "scatter"}, {"type": "heatmap"}]] * n_panels
    )

    for i, panel in enumerate(panel_names, start=1):
        panel_cell_types = df[df['cell_panel'] == panel]["cell_type"].unique()
        panel_cell_types = sorted(panel_cell_types, key=lambda ct: cell_type_class_map.get(ct, ""))
        heatmap_data = global_heatmap.loc[global_heatmap.index.intersection(panel_cell_types)]
        y_labels = heatmap_data.index.astype(str).tolist()

        # For hover
        hover_class_data = np.array([[cell_type_class_map[ct]] * len(heatmap_data.columns) for ct in y_labels])

        # --- Colorbar scatter ---
        fig.add_trace(go.Scatter(
            x=[0] * len(y_labels),
            y=y_labels,
            mode="markers",
            marker=dict(
                color=[cell_type_to_color.get(ct, "gray") for ct in y_labels],
                size=10
            ),
            text=[cell_type_class_map.get(ct, "") for ct in y_labels],
            hovertemplate="Cell class: %{text}<extra></extra>",
            showlegend=False
        ), row=i, col=1)

        # --- Heatmap ---
        # Calculate y position in normalized paper coordinates for the **top** of the current row
        # Default base y position (top of current panel)
        base_y_pos = 1 - sum(normalized_heights[:i-1]) - 0.01
        
        # Apply per-panel tweak
        if i == 3:
            y_pos = base_y_pos - 0.025 # slightly lower
        elif i == 4:
            y_pos = base_y_pos - 0.015
        elif i == 5:
            y_pos = base_y_pos - 0.015
        else:
            y_pos = base_y_pos  # default for other panels



        fig.add_trace(go.Heatmap(
            z=heatmap_data.values,
            x=heatmap_data.columns,
            y=y_labels,
            customdata=hover_class_data,
            colorscale=colorscale_theme,
            zmin=global_zmin,
            zmax=global_zmax,
            showscale=True,
            colorbar=dict(
                title="Mean Expr",
                orientation="h",
                x=-1,
                xanchor="left",
                y=y_pos,
                len=0.8,
                thickness=12
            ),
            hovertemplate="Tissue: %{x}<br>Cell type: %{y}<br>Class: %{customdata}<br>Expr: %{z:.2f}<extra></extra>"
        ), row=i, col=2)
        fig.update_yaxes(
            title_text=panel,
            title_font=dict(size=13, color="black"),
            title_standoff=10,
            row=i,
            col=1
        )

        # Axes formatting
        fig.update_xaxes(showticklabels=False, row=i, col=1)
        fig.update_xaxes(showticklabels=True, tickangle=45, tickfont=dict(size=10), row=i, col=2)
        fig.update_yaxes(autorange="reversed", tickfont=dict(size=9), row=i, col=1)
        fig.update_yaxes(autorange="reversed", tickfont=dict(size=9), row=i, col=2)

    # Final layout with your spacing and background customizations
    fig.update_layout(
        width=600,
        height=total_rows * 20,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=10, r=10, t=50, b=150),
        title=dict(
            text=f"{gene_label} Expression",
            x=0.5,
            xanchor="center",
            yanchor="top",
            font=dict(size=20, color="black")
        ),
        annotations=[
        dict(
            text=f"{gene_label} Expression",
            x=0.5,
            y=-0.03,  # y=0 refers to the bottom of the plotting area
            xref="paper",
            yref="paper",
            xanchor="center",
            yanchor="bottom",
            showarrow=False,
            font=dict(size=20, color="black")
        )
    ]
    )

    # Save
    output_path = f"{output_dir}heatmap/{gene_label}.html"
    fig.write_html(output_path, include_plotlyjs="cdn")
    print(f"[✓] Heatmap written: {output_path}")


def generate_heatmap_wrapper(gene_id):
    try:
        gene_label = gene_label_map[gene_id]
        expr = gene_expr_map[gene_id]
        generate_heatmap(gene_id, gene_label, expr)
    except Exception as e:
        print(f"[!] Error with {gene_id}: {e}")

Parallel(n_jobs=16, backend="threading")(
    delayed(generate_heatmap_wrapper)(processTabulaSapiens.gene_id)
    for gene_id in tqdm(processTabulaSapiens.valid_gene_ids, desc="Generating heatmaps")
)
 
    
        