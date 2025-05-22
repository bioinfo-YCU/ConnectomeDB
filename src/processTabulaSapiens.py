## Function to preprocess (normalize, log) and perform PCA then UMAP to the Tabula Sapiense dataset 
import scanpy as sc
import anndata
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import requests
import scanpy as sc
import re
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

