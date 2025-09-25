# concat all species DB in one long table (csv/excel)

import os
import glob
import pandas as pd

# Parent subdirectory containing CSV/ and Excel/ folders
s_directory = "CDB2025v1/"
p_directory = "downloads/"
directory = os.path.join(p_directory, s_directory)

# Desired order
custom_order = ["all", "human", "mouse", "chimp", "macaque", "marmoset", 
                "rat", "pig", "cow", "dog", "horse", "sheep", "chicken", 
                "frog", "zebrafish"]

# Species-specific ID suffixes
id_suffixes = ["HGNC ID", "MGI ID", "XX ID", "RGD ID", "XEN ID", "ZFIN ID"]

# Detect subfolders "CSV" and "Excel"
for subfolder in ["CSV", "Excel"]:
    sub_dir_path = os.path.join(directory, subfolder)
    if not os.path.exists(sub_dir_path):
        continue  # skip if folder doesn't exist

    print(f"Processing {subfolder} folder...")

    # Detect files
    csv_files = glob.glob(os.path.join(sub_dir_path, "*.csv"))
    excel_files = glob.glob(os.path.join(sub_dir_path, "*.xlsx"))
    
    # Ignore existing "all_species" outputs
    csv_files = [f for f in csv_files if os.path.basename(f) != "all_species.csv"]
    excel_files = [f for f in excel_files if os.path.basename(f) != "all_species.xlsx"]

    all_files = csv_files + excel_files
    if not all_files:
        print(f"No files found in {sub_dir_path}")
        continue

    dfs = []
    file_dict = {}
    
    # Map species -> file
    for file in all_files:
        filename = os.path.basename(file)
        name_no_ext = os.path.splitext(filename)[0]
        species = name_no_ext.split("_")[-1]
        file_dict[species] = file

    # Process files in custom order
    for species in custom_order:
        if species in file_dict:
            file = file_dict[species]
            if file.endswith(".csv"):
                df = pd.read_csv(file)
            else:  # Excel
                df = pd.read_excel(file)
            
            # Rename only species-specific ID columns
            rename_dict = {}
            for col in df.columns:
                for suffix in id_suffixes:
                    if col == f"Ligand {suffix}":
                        rename_dict[col] = "Ligand Species ID"
                    elif col == f"Receptor {suffix}":
                        rename_dict[col] = "Receptor Species ID"
            df = df.rename(columns=rename_dict)

            # Add Species column
            df["Species"] = species
            dfs.append(df)

    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)

        # Save/overwrite
        if subfolder == "CSV":
            output_file = os.path.join(sub_dir_path, "all_species.csv")
            combined_df.to_csv(output_file, index=False)
        else:  # Excel
            output_file = os.path.join(sub_dir_path, "all_species.xlsx")
            combined_df.to_excel(output_file, index=False)

        print(f"Created {output_file} with {len(combined_df)} rows")
