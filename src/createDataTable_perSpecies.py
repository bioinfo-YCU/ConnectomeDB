## Function to prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the database qmds
import sys, os
from itables import init_notebook_mode
import pandas as pd
from itables import show
from itables import options
from IPython.display import HTML, display
import numpy as np
from bs4 import BeautifulSoup
from createDataTable import gene_pair, gene_pair000, human_columns, lrPairsCount
import warnings

# Suppress SettingWithCopyWarning
warnings.simplefilter("ignore", category=UserWarning)

# for mouse-specific, fill out missing info



# Other vertebrates
species_list = [
    "ptroglodytes", "ggallus", "sscrofa", "btaurus", 
    "clfamiliaris", "ecaballus", "oarambouillet", "cjacchus",
    "mmulatta"
]

def process_species(gene_pair_df, gene_pair000_df, species, id_prefix, ligand_index, receptor_index):
    """
    Processes ligand-receptor interactions for a given species.

    Parameters:
        gene_pair_df (pd.DataFrame): The main gene pair DataFrame.
        gene_pair000_df (pd.DataFrame): The filtered gene pair DataFrame.
        species (str): The species name (e.g., "Mouse", "Rat", "Zebrafish").
        id_prefix (str): The identifier prefix for the species (e.g., "MGI" for Mouse, "RGD" for Rat).
        ligand_index (int): Index for selecting the ligand column.
        receptor_index (int): Index for selecting the receptor column.

    Returns:
        pd.DataFrame: Processed DataFrame for the species.
    """
    # Identify relevant columns for the species
    if species in ["Mouse", "Rat", "Zebrafish"]:
        species_columns = [col for col in gene_pair_df.columns if id_prefix in col or species in col]
    else:
        species_columns = [col for col in gene_pair_df.columns if species in col]

    # Filter rows where all species-specific columns are not empty
    species_gene_pair = gene_pair000_df[(gene_pair000_df[species_columns].map(str.strip) != "").all(axis=1)]
    
    # Rename columns to remove species name
    species_gene_pair.columns = [
        col.replace(f"{species} ", "").strip() if "Ligand" in col or "Receptor" in col else col
        for col in species_gene_pair.columns
    ]

    # Identify ligand and receptor columns dynamically
    ligand_col = [col for col in species_gene_pair.columns if "Ligand&nbsp;" in col][ligand_index]
    receptor_col = [col for col in species_gene_pair.columns if "Receptor&nbsp;" in col][receptor_index]
    ligand_Location = [col for col in species_gene_pair.columns if "Ligand Location" in col][0]
    receptor_Location = [col for col in species_gene_pair.columns if "Receptor Location" in col][0]

    if id_prefix == "Ensembl":
        # Extract species names dynamically
        all_species = {"Chimpanzee", "Pig", "Dog", "Cow", 
                       "Chicken", "Horse", "Sheep", "Marmoset", "Macaque"}
        # Function to clean column names by removing HTML tags
        def clean_column_name(col):
            # Remove HTML tags using BeautifulSoup
            cleaned_name = BeautifulSoup(col, "html.parser").get_text()
            return cleaned_name.strip()
        
        # Step 1: Identify columns that contain any species from all_species
        cleaned_columns = [clean_column_name(col) for col in species_gene_pair.columns]
        
        # Identify columns where the species from all_species is in the cleaned name
        columns_to_remove = [
            species_gene_pair.columns[i] for i, cleaned_name in enumerate(cleaned_columns)
            if any(species in cleaned_name for species in all_species)
        ]
        
        # Step 2: Remove those columns from the DataFrame
        species_gene_pair = species_gene_pair.drop(columns=columns_to_remove)
        
        
            # Function to format ligand-receptor pairs
    def format_lr_pair(row):
        if row[ligand_Location] in ['secreted', '']:
            return f"{row[ligand_col]} <span style='font-size: 15px;'>○</span> <span style='font-size: 24px;'>⤚</span> {row[receptor_col]}"
        elif row[receptor_Location] == 'plasma membrane':
            return f"{row[ligand_col]} <span style='font-size: 24px;'>⤙</span> <span style='font-size: 24px;'>⤚</span> {row[receptor_col]}"
        else:
            return f"{row[ligand_col]} \u2192 {row[receptor_col]}"

    # Apply the formatting function
    species_gene_pair1 = species_gene_pair.copy()
    species_gene_pair1.loc[:, f"{species} LR Pair"] = species_gene_pair1.apply(format_lr_pair, axis=1)

    # Identify relevant columns for the species
    species_columns = [col for col in species_gene_pair1.columns if id_prefix in col]
    new_order = [human_columns[0]]+ [f"{species} LR Pair", ligand_col, receptor_col] + species_columns + human_columns[1:]
    species_gene_pair1 = species_gene_pair1[new_order].reset_index(drop=True)
    # Turn Interaction ID into clickable links
    # species_gene_pair1[species_gene_pair1.columns[0]] = species_gene_pair1[species_gene_pair1.columns[0]].apply(
    #     lambda x: f"<a href='https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/database/filter/{x}.html'>{x}</a>"
    # )
    col_to_rename = [col for col in species_gene_pair1.columns if "Human LR Pair" in col][0]
    species_gene_pair1.rename(columns={col_to_rename: "LR Pair Card"}, inplace=True)

    return species_gene_pair1


# Process each species
# Mouse
mouse_gene_pair1 = process_species(gene_pair, gene_pair000, "Mouse", "MGI", 1, 1)
MouselrPairsCount = len(mouse_gene_pair1["Mouse LR Pair"]) # non-unique for now
HumanMouseLRPairsPer = (MouselrPairsCount/lrPairsCount)*100
HumanMouseLRPairsPer = round(HumanMouseLRPairsPer, 2)
# Round up to the nearest 0.5%
HumanMouseLRPairsPer = ((HumanMouseLRPairsPer * 2 + 1) // 1) / 2

mouse_gene_pair1 = mouse_gene_pair1.reset_index(drop=True)  
# Rat
rat_gene_pair1 = process_species(gene_pair, gene_pair000, "Rat", "RGD", 2, 2)
# Zebrafish
zebrafish_gene_pair1 = process_species(gene_pair, gene_pair000, "Zebrafish", "ZFIN", 3, 3)

# Chimpanzee
chimpanzee_gene_pair1 = process_species(gene_pair, gene_pair000, "Chimpanzee", "Ensembl", 4, 4)

# Chicken
chicken_gene_pair1 = process_species(gene_pair, gene_pair000, "Chicken", "Ensembl", 5, 5)

# Pig
pig_gene_pair1 = process_species(gene_pair, gene_pair000, "Pig", "Ensembl", 6, 6)

# Cow
cow_gene_pair1 = process_species(gene_pair, gene_pair000, "Cow", "Ensembl", 7, 7)

# Dog
dog_gene_pair1 = process_species(gene_pair, gene_pair000, "Dog", "Ensembl", 8, 8)

# Horse
horse_gene_pair1 = process_species(gene_pair, gene_pair000, "Horse", "Ensembl", 9, 9)

# Sheep
sheep_gene_pair1 = process_species(gene_pair, gene_pair000, "Sheep", "Ensembl", 10, 10)

# Marmoset
marmoset_gene_pair1 = process_species(gene_pair, gene_pair000, "Marmoset", "Ensembl", 11, 11)

# Rhesus macaque
macaque_gene_pair1 = process_species(gene_pair, gene_pair000, "Macaque", "Ensembl", 12, 12)

