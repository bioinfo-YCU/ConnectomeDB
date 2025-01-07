## Function to convert MGI/RGD IDs/Other Species to gene symbols via biomart
import os
import sys
from biomart import BiomartServer
import pandas as pd

# Add source directory to the path
sys.path.append(os.path.abspath("src"))
from createDataTable import pop_up_info_lim

# Species-specific parameters
species_id_prefix = "ZFIN"  # "RGD" for rat, "MGI" for mouse
dataset_name = "drerio_gene_ensembl" # "drerio_gene_ensembl" for zebra fish, "rnorvegicus_gene_ensembl" for rat, "mmusculus_gene_ensembl" for mouse
gene_symbol_field = "zfin_symbol"  # "external_gene_name" for mouse

# Access Biomart server and dataset
biomart_server_url = "http://www.ensembl.org/biomart"
server = BiomartServer(biomart_server_url)
dataset = server.datasets[dataset_name]

# Fetch species IDs from the dataset
species_ids = pop_up_info_lim[species_id_prefix + ' ID'].unique()
species_ids = [id for id in species_ids if id.startswith(species_id_prefix + ":")]

