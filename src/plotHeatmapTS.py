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

    # Only 1 column (heatmap), not 2
    fig = make_subplots(
        rows=n_panels, cols=1,
        shared_xaxes=False,
        shared_yaxes=True,
        vertical_spacing=0.035,
        row_heights=normalized_heights,
        specs=[[{"type": "heatmap"}]] * n_panels
    )

    for i, panel in enumerate(panel_names, start=1):
        panel_cell_types = df[df['cell_panel'] == panel]["cell_type"].unique()
        panel_cell_types = sorted(panel_cell_types, key=lambda ct: cell_type_class_map.get(ct, ""))
        heatmap_data = global_heatmap.loc[global_heatmap.index.intersection(panel_cell_types)]
        y_labels = heatmap_data.index.astype(str).tolist()

        hover_class_data = np.array([[cell_type_class_map[ct]] * len(heatmap_data.columns) for ct in y_labels])

        # --- Heatmap only ---
        base_y_pos = 1 - sum(normalized_heights[:i-1])
        y_pos = base_y_pos
        if i == 2:
            y_pos = base_y_pos - 0.010
        if i == 3:
            y_pos = base_y_pos - 0.035
        elif i == 4:
            y_pos = base_y_pos - 0.023
        elif i == 5:
            y_pos = base_y_pos - 0.026
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
                title="Log(Gene Expression)",
                orientation="h",
                x=-0.7,
                xanchor="left",
                y=y_pos,
                len=0.9,
                thickness=10
            ),
            hovertemplate="Tissue: %{x}<br>Cell type: %{y}<br>Class: %{customdata}<br>Expr: %{z:.2f}<extra></extra>"
        ), row=i, col=1)

        # Axes
        fig.update_xaxes(#title_text=panel,
                         showticklabels=True, 
                         tickangle=270, 
                         tickfont=dict(size=9), row=i, col=1)
        fig.update_yaxes(
            autorange="reversed",
            tickfont=dict(size=9),
            title_text=panel,
            title_font=dict(size=13, color="black"),
            title_standoff=5,
            row=i,
            col=1
        )
        # Compute the center y-position of this panel (in data coordinates)
        mid_index = len(y_labels) // 2
        mid_cell = y_labels[mid_index] if y_labels else ""
        
        # Compute the center y-position of this panel (in data coordinates)
        mid_index = len(y_labels) // 2
        mid_cell = y_labels[mid_index] if y_labels else ""
        
        fig.add_annotation(
            text=f"<b>{panel}</b>",
            xref=f'x{i}',  # panel's x-axis
            yref=f'y{i}',  # panel's y-axis
            x=heatmap_data.columns[len(heatmap_data.columns) // 2],  # center tissue
            y=mid_cell,  # center cell type
            showarrow=False,
            font=dict(size=14, color='black'),
            xanchor="center",
            yanchor="middle"
        )


    # Layout
    fig.update_layout(
        width=585,
        height=total_rows * 20,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=1, r=1, t=80, b=100),
        title=dict(
            text=f"{gene_label} Expression",
            x=0.7,
            y= 0.995,
            xanchor="center",
            yanchor="top",
            font=dict(size=18, color="black")
        ),
        annotations=[
            dict(
                text=f"{gene_label} Expression",
                x=0.5,
                y=-0.03,
                xref="paper",
                yref="paper",
                xanchor="center",
                yanchor="bottom",
                showarrow=False,
                font=dict(size=18, color="black")
            )
        ]
    )
    # Add panel titles globally above each subplot
    for i, panel in enumerate(panel_names, start=1):
        # same as above
        base_y_pos = 1 - sum(normalized_heights[:i-1])
        y_pos = base_y_pos
        if i == 2:
            y_pos = base_y_pos - 0.010
        if i == 3:
            y_pos = base_y_pos - 0.035
        elif i == 4:
            y_pos = base_y_pos - 0.023
        elif i == 5:
            y_pos = base_y_pos - 0.026
    
        fig.add_annotation(
            text=f"<b>{panel}</b>",
            x=0.5,
            y=y_pos,
            xref="paper",
            yref="paper",
            showarrow=False,
            font=dict(size=14),
            xanchor="center",
            yanchor="bottom"
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
 
    
        