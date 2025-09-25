# concat all species DB in one long table (JSON)
import os
import glob
import pandas as pd

# Define the relative path to your JSON files
s_directory = "CDB2025v1/"
p_directory = "downloads/"
directory = os.path.join(p_directory, s_directory)

# Desired order for concatenation
custom_order = ["all", "human", "mouse", "chimp", "macaque", "marmoset", 
                "rat", "pig", "cow", "dog", "horse", "sheep", "chicken", 
                "frog", "zebrafish"]

# Species-specific ID suffixes, adapted for JSON column names with underscores
id_suffixes = ["HGNC_ID", "MGI_ID", "XX_ID", "RGD_ID", "XEN_ID", "ZFIN_ID"]

# Processing only the "JSON" subfolder
subfolder = "JSON"
sub_dir_path = os.path.join(directory, subfolder)

# Check if the JSON folder exists and create it if necessary
if not os.path.exists(sub_dir_path):
    os.makedirs(sub_dir_path, exist_ok=True)
    print(f"Created the directory: {sub_dir_path}")
    
print(f"Processing {subfolder} folder...")

# Get a list of all JSON files in the directory
all_files = glob.glob(os.path.join(sub_dir_path, "*.json"))

# Exclude the output file from the list of files to be processed
all_files = [f for f in all_files if os.path.basename(f) != "all_species.json"]

print(f"Found {len(all_files)} JSON file(s) to process.")

if not all_files:
    print(f"No JSON files were found in {sub_dir_path}. Please check your file paths and names.")
else:
    dfs = []
    file_dict = {}
    
    # Map species to their corresponding file paths
    for file in all_files:
        filename = os.path.basename(file)
        name_no_ext = os.path.splitext(filename)[0]
        # NEW LOGIC: Extract species based on the "name_gene_pair" pattern
        species_part = name_no_ext.replace('_gene_pair', '')
        # Filter for known species to avoid adding unexpected files
        if species_part in custom_order:
            species = species_part
            file_dict[species] = file
        else:
            print(f"Warning: Filename '{filename}' does not match the expected species format. Skipping.")
            continue
            
    # Process files in the predefined custom order
    for species in custom_order:
        if species in file_dict:
            file = file_dict[species]
            
            try:
                print(f"Processing file: {os.path.basename(file)}...")
                # Read the JSON file into a pandas DataFrame
                df = pd.read_json(file)
                
                # Normalize (flatten) any nested JSON data into a flat table
                df = pd.json_normalize(df.to_dict('records'))
                
                # Prepare a dictionary for renaming columns
                rename_dict = {}
                for suffix in id_suffixes:
                    ligand_col = f"Ligand_{suffix}"
                    if ligand_col in df.columns:
                        rename_dict[ligand_col] = "Ligand Species ID"
                    
                    receptor_col = f"Receptor_{suffix}"
                    if receptor_col in df.columns:
                        rename_dict[receptor_col] = "Receptor Species ID"
                
                # Rename the columns if any matching ones are found
                if rename_dict:
                    df = df.rename(columns=rename_dict)
                
                # Add a new column to identify the species for each row
                df["Species"] = species
                dfs.append(df)
                
            except Exception as e:
                print(f"Error processing file {os.path.basename(file)}: {e}")
                continue
            
    if dfs:
        # Concatenate all DataFrames into a single one
        combined_df = pd.concat(dfs, ignore_index=True)
        output_file = os.path.join(sub_dir_path, "all_species.json")
        
        # Print the absolute path so you know exactly where to find the file
        print(f"Saving combined file to: {os.path.abspath(output_file)}")
        
        # Save the final DataFrame as a new JSON file
        combined_df.to_json(output_file, orient="records", indent=4)
        
        print(f"Created {output_file} with {len(combined_df)} rows")
    else:
        print("No files were successfully processed. The 'dfs' list is empty.")