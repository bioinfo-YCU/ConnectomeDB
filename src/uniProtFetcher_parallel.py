import pandas as pd
import requests
from io import StringIO
import os
from multiprocessing import Pool, cpu_count

# Specify a list of species here to process in parallel
species_to_process = ["mouse", "rat", "human", "zebrafish", "chimpanzee", "chicken", "pig", "cow", "dog", "horse", "marmoset",   "macaque"]

                      # === INTERNAL MAPPINGS ===
species_info = {
    "mouse":        {"taxid": "10090", "code": "mmusculus"},
    "rat":          {"taxid": "10116", "code": "rnorvegicus"},
    "zebrafish":    {"taxid": "7955",  "code": "drerio"},
    "chimpanzee":   {"taxid": "9598",  "code": "ptroglodytes"},
    "chicken":      {"taxid": "9031",  "code": "ggallus"},
    "pig":          {"taxid": "9823",  "code": "sscrofa"},
    "cow":          {"taxid": "9913",  "code": "btaurus"},
    "dog":          {"taxid": "9615",  "code": "clfamiliaris"},
    "horse":        {"taxid": "9796",  "code": "ecaballus"},
    "sheep":        {"taxid": "9940",  "code": "oarambouillet"},
    "marmoset":     {"taxid": "9483",  "code": "cjacchus"},
    "macaque":      {"taxid": "9544",  "code": "mmulatta"},
    "human":        {"taxid": "9606",  "code": "hsapiens"},
}

# Function to process a single species
def process_species(species_input):
    """
    Fetches UniProt data for a given species, processes it, and saves to a CSV file.
    """
    if species_input not in species_info:
        print(f"Species '{species_input}' not supported. Skipping.")
        return

    taxid = species_info[species_input]["taxid"]

    # === UNIPROT API SETTINGS ===
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    query = f"reviewed:true AND taxonomy_id:{taxid}"
    fields = "accession,gene_names,protein_name,organism_name"
    format = "tsv"

    params = {
        "query": query,
        "fields": fields,
        "format": format,
        "size": 500  # Max allowed is 500 per request; use paging if needed
    }

    # === API REQUEST ===
    print(f"Querying UniProt for Swiss-Prot entries for {species_input} (taxid={taxid})...")
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise an HTTPError for bad responses (4xx or 5xx)
    except requests.exceptions.RequestException as e:
        print(f"Failed to fetch data for {species_input}: {e}")
        return

    # === PARSE AND SAVE ===
    df = pd.read_csv(StringIO(response.text), sep="\t")

    # === SPLIT MULTIPLE GENE NAMES INTO ROWS ===
    df["Gene Names"] = df["Gene Names"].fillna("")  # Just in case
    df_expanded = df.assign(**{
        "Gene Names": df["Gene Names"].str.split()
    }).explode("Gene Names").reset_index(drop=True)

    # Ensure output directory exists
    os.makedirs("data", exist_ok=True)

    # Save to file
    out_file = f"data/uniprotMapping_{species_input}.csv"
    df_expanded.to_csv(out_file, index=False)

    print(f"Saved {len(df)} entries for {species_input} to {out_file}")

if __name__ == "__main__":
    # === USER INPUT for species to process ===

    # Or to process all species in your mapping:
    # species_to_process = list(species_info.keys())

    # Determine the number of processes to use
    # It's often good practice to use cpu_count() or cpu_count() - 1
    num_processes = cpu_count()
    print(f"Starting parallel processing with {num_processes} processes...")

    # Create a pool of workers
    with Pool(processes=num_processes) as pool:
        # Map the process_species function to the list of species
        pool.map(process_species, species_to_process)

    print("All species processing complete.")