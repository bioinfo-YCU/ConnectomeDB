## Function to preprocess (normalize, log) and perform PCA then UMAP to the Tabula Sapiense dataset 
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
from createDataTable import gene_pair0, pop_up_info

output_dir = "data/tabula_sapiens/"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(output_dir +"umap" , exist_ok=True)
os.makedirs(output_dir +"heatmap" , exist_ok=True)

output_file = "data/tissue_dataset.h5ad"

if os.path.exists(output_file):
    # Do only if the file exists
    print("File exists. Proceeding with the task.")

    # Rest of the logic that should only run if file exists
else:
    url = "https://datasets.cellxgene.cziscience.com/9daa676b-07ec-4cea-80aa-daa49200aa64.h5ad"
    #Tabula Sapiens is a benchmark, first-draft human cell atlas of over 1.1M cells from 28 organs of 24 normal human subjects. This work is the product of the Tabula Sapiens Consortium. Taking the organs from the same individual controls for genetic background, age, environment, and epigenetic effects, and allows detailed analysis and comparison of cell types that are shared between tissues.
    # Get file size for progress bar
    response = requests.head(url)
    total_size = int(response.headers.get('Content-Length', 0))
    
    # Stream download with tqdm
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(output_file, 'wb') as f, tqdm(
            total=total_size, unit='B', unit_scale=True, desc=output_file
        ) as pbar:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
                    
# Load with scanpy
adata = sc.read_h5ad(output_file, backed='r')
print(adata)
print(adata.obs.columns) 
print(adata.var_names)    # Gene names

### Looks like the UMAP has been computed ###
# normalize all data and log 
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# sc.pp.pca(adata)
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)

# Build HGNC-to-Ensembl mapping
ensembl_id = dict(zip(pop_up_info['HGNC ID'], pop_up_info['ensembl_gene_id']))

# Define HGNC ID extractor
def extract_hgnc_id(text):
    if pd.isna(text): return None
    match = re.search(r'(HGNC:\d+)', str(text))
    return match.group(1) if match else None

# --- Prepare ligand dataframe ---
ligand_df = gene_pair0[['Ligand', 'Ligand HGNC ID']].copy()
ligand_df.columns = ['gene_symbol', 'hgnc_id']  # Standardize column names

# --- Prepare receptor dataframe ---
receptor_df = gene_pair0[['Receptor', 'Receptor HGNC ID']].copy()
receptor_df.columns = ['gene_symbol', 'hgnc_id']

# Combine both
gene_pair_input = pd.concat([ligand_df, receptor_df], ignore_index=True)

# Extract clean HGNC IDs
gene_pair_input['hgnc_id'] = gene_pair_input['hgnc_id'].apply(extract_hgnc_id)

# Map to Ensembl IDs
gene_pair_input['ensembl_id'] = gene_pair_input['hgnc_id'].map(ensembl_id)

# Drop duplicates and NaNs if needed
gene_pair_input = gene_pair_input.drop_duplicates().dropna(subset=['ensembl_id'])
# Final output
# print(gene_pair_input.head())

# If you want to save the scaled data again
# --- Precompute Scaled Expression ---
def precompute_scaled_expression(adata, gene_id_list):
    if "scale_data" not in adata.layers:
        raise ValueError("Layer 'scale_data' not found in adata.")

    X = adata.layers["scale_data"]  # shape: cells x genes
    var_names_stripped = [strip_version(g) for g in adata.var_names]
    gene_indices = {gene_id: i for i, gene_id in enumerate(var_names_stripped)}

    gene_expr_map = {}
    missing_genes = []

    for gene_id in tqdm(gene_id_list, desc="Precomputing scaled expression"):
        if gene_id not in gene_indices:
            missing_genes.append(gene_id)
            continue
        i = gene_indices[gene_id]
        x = X[:, i]  # ✅ not transposed
        expr = x.toarray().flatten() if hasattr(x, "toarray") else x.flatten()
        gene_expr_map[gene_id] = expr

    print(f"[✓] Precomputed expression for {len(gene_expr_map)} genes.")
    if missing_genes:
        print(f"[!] Skipped {len(missing_genes)} genes not found in adata.var_names.")

    return gene_expr_map, missing_genes


# uncomment for now so it is not accidentally ran
# gene_expr_map, _ = precompute_scaled_expression(adata, gene_id_list_stripped)
# If you want to save the scaled data again
# with open("data/gene_expr_map_scaled.pkl", "wb") as f:
#     pickle.dump(gene_expr_map, f)

X = adata.layers["scale_data"]  # shape: cells x genes
var_names_stripped = [strip_version(g) for g in adata.var_names]
gene_indices = {gene_id: i for i, gene_id in enumerate(var_names_stripped)}

# --- Normalize Ensembl IDs and build label map ---
def strip_version(ensembl_id):
    return ensembl_id.split('.')[0]

# Rebuild gene_label_map using stripped Ensembl IDs
gene_label_map = {
    strip_version(row["ensembl_id"]): row["gene_symbol"]
    for _, row in gene_pair_input.iterrows()
}

# Strip version from gene_id_list
gene_id_list = [strip_version(g) for g in gene_id_list]

# --- Filter to valid genes ---
valid_gene_ids = [g for g in gene_id_list if g in gene_indices]
rows = [gene_indices[g] for g in valid_gene_ids]
missing_genes = [g for g in gene_id_list if g not in gene_indices]
print(f"[!] Skipped {len(missing_genes)} genes not found in adata.var_names.")
#len(valid_gene_ids)

