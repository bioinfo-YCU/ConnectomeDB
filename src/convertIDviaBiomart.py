## Function to convert MGI/RGD IDs/Other Species to gene symbols via biomart
import os
import sys
from biomart import BiomartServer
import pandas as pd

# Add source directory to the path
sys.path.append(os.path.abspath("src"))
from createDataTable import pop_up_info_lim

# Species-specific parameters
species_id_prefix = "RGD"  # "RGD" for rat, "MGI" for mouse
dataset_name = "rnorvegicus_gene_ensembl"  # For rat, "mmusculus_gene_ensembl" for mouse
gene_symbol_field = "rgd_symbol"  # "external_gene_name" for mouse

# Access Biomart server and dataset
biomart_server_url = "http://www.ensembl.org/biomart"
server = BiomartServer(biomart_server_url)
dataset = server.datasets[dataset_name]

# Fetch species IDs from the dataset
species_ids = pop_up_info_lim[species_id_prefix + ' ID'].unique()
species_ids = [id for id in species_ids if id.startswith(species_id_prefix + ":")]

# Clean the IDs if they are from RGD
if species_id_prefix == "RGD":
    species_ids = [value.replace("RGD:", "") for value in species_ids]

# Function to query Biomart in chunks
def query_in_chunks(ids, chunk_size=100):
    gene_mapping = {}
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i + chunk_size]
        try:
            response = dataset.search({
                'filters': {species_id_prefix.lower() + '_id': chunk},
                'attributes': [species_id_prefix.lower() + '_id', gene_symbol_field]
            })
            # Parse the response and add to gene_mapping
            for line in response.iter_lines(decode_unicode=True):
                species_id, gene_name = line.split("\t")
                gene_mapping[species_id] = gene_name
        except Exception as e:
            print(f"Error processing chunk {i // chunk_size + 1}: {e}")
    return gene_mapping

# Query the dataset
gene_mapping = query_in_chunks(species_ids)

# Display results
for species_id, gene_name in gene_mapping.items():
    print(f"{species_id}: {gene_name}")

# Convert the gene_mapping dictionary to a DataFrame and save as CSV
df = pd.DataFrame.from_dict(gene_mapping, orient='index', columns=[species_id_prefix + ' name']) \
    .reset_index() \
    .rename(columns={'index': species_id_prefix + ' ID'})

# Save the DataFrame to a CSV file
df.to_csv(f"data/{species_id_prefix}_ID_biomart.csv", index=False)