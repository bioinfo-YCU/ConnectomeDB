## Function to convert Humans (HGNC symbol) to Other species ID via biomart
import os
import sys
from biomart import BiomartServer
import pandas as pd

# Add source directory to the path
sys.path.append(os.path.abspath("src"))
from createDataTable import hgnc_id

# Species-specific parameters
species_id_prefix = "ZFIN"  # "RGD" for rat, "MGI" for mouse
dataset_name = "drerio_gene_ensembl" # "drerio_gene_ensembl" for zebra fish, "rnorvegicus_gene_ensembl" for rat, "mmusculus_gene_ensembl" for mouse
gene_id_field = "zfin_id"  # "external_gene_name" for mouse
gene_symbol_field = "zfin_symbol"  # "external_gene_name" for mouse
output_filename="data/hgnc_to_"+species_id_prefix+"_mapping.csv"

# Define function for conversion
def convert_hgnc_to_zfin(hgnc_ids, output_file=None):
    """
    Converts a list of HGNC IDs (human) to Other Species (e.g. ZFIN IDs (zebrafish)) using Ensembl Biomart.
    
    Args:
        hgnc_ids (list): A list of HGNC IDs to be converted.
        output_file (str): Path to save the conversion results (optional).
    
    Returns:
        pd.DataFrame: DataFrame containing input HGNC IDs and their corresponding ZFIN IDs.
    """
    # Biomart server configuration
    biomart_server_url = "http://www.ensembl.org/biomart"
    server = BiomartServer(biomart_server_url)
    
    # Select zebrafish dataset
    dataset_name = dataset_name
    dataset = server.datasets[dataset_name]

    # Query attributes
    attributes = [
        "ensembl_gene_id",  # Ensembl Gene ID
        "hgnc_id",          # HGNC ID (human)
         gene_id_field,      # Species to convert to
         gene_symbol_field
    ]
    
    # Build query
    response = dataset.search({
        "filters": {
            "hgnc_id": hgnc_ids
        },
        "attributes": attributes,
    })

    # Parse response to DataFrame
    results = pd.read_csv(response, sep="\t", header=None, names=attributes)

    # Save results if output_file is provided
    if output_file:
        results.to_csv(output_file, index=False)

    return results

# Example usage
if __name__ == "__main__":
    # Call function and print results
    result_df = convert_hgnc_to_zfin(hgnc_id, output_filename)
    print(result_df)
