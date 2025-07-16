## Function to prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the database qmds
import sys, os
from itables import init_notebook_mode
import pandas as pd
from itables import show
from itables import options
from IPython.display import HTML, display
import numpy as np
import re
from bs4 import BeautifulSoup
from createDataTable import gene_pair, gene_pair000, human_columns, lrPairsCount
import warnings


# Suppress SettingWithCopyWarning
warnings.simplefilter("ignore", category=UserWarning)

def cleanup_species_info(species_info, species_name):
    """
    Simple cleanup: replace empty associated_gene_name with homolog_ensembl_gene
    """
    # Find the associated gene name column
    associated_gene_col = None
    homolog_ensembl_col = None
    
    for col in species_info.columns:
        if species_name in col and 'homolog' in col and 'associated_gene_name' in col:
            associated_gene_col = col
        elif species_name in col and 'homolog_ensembl_gene' in col:
            homolog_ensembl_col = col
    
    if associated_gene_col and homolog_ensembl_col:
        # Replace empty or NaN values in associated_gene_name with homolog_ensembl_gene
        mask = (species_info[associated_gene_col] == "") | (species_info[associated_gene_col].isna())
        species_info.loc[mask, associated_gene_col] = species_info.loc[mask, homolog_ensembl_col]
    
    return species_info
    
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
        species_gene_pair = species_gene_pair.drop(columns=columns_to_remove)
        
    
    # Apply the formatting function
    species_gene_pair1 = species_gene_pair.copy()

    # Identify ligand and receptor columns dynamically
    ligand_col = [col for col in species_gene_pair1.columns if "Ligand&nbsp;" in col][ligand_index]
    receptor_col = [col for col in species_gene_pair1.columns if "Receptor&nbsp;" in col][receptor_index]
    ligand_Location = [col for col in species_gene_pair1.columns if "Ligand Location" in col][0]
    receptor_Location = [col for col in species_gene_pair1.columns if "Receptor Location" in col][0]
    # Identify relevant columns for the species
    species_columns = [col for col in species_gene_pair1.columns if id_prefix in col]
    new_order = [human_columns[0]]+ [ligand_col, receptor_col] + species_columns + human_columns[1:]
    species_gene_pair1 = species_gene_pair1[new_order].reset_index(drop=True)

    # Apply the formatting function
    #species_gene_pair1 = species_gene_pair.copy()
            
    def extract_hgnc_id(col):
        """Use regular expression to extract the HGNC ID after 'HGNC:'."""
        match = re.search(r'HGNC:(\d+)', col)
        if match:
            return "HGNC:"+match.group(1)
        return None
    #clean HGNC
    

    species_name = {
        "Mouse": "mmusculus",
        "Rat": "rnorvegicus",
        "Zebrafish":"drerio" ,
        "Chimpanzee":"ptroglodytes",
        "Chicken":"ggallus",
        "Pig":"sscrofa",
        "Cow":"btaurus",
        "Dog":"clfamiliaris",
        "Horse":"ecaballus",
        "Sheep":"oarambouillet",
        "Marmoset": "cjacchus" ,
        "Macaque": "mmulatta"   
    }.get(species, "Unknown species")
    
    # Load species-specific data
    species_info = pd.read_csv(f"data/{species_name}_ID_biomart.csv")
    # replace empty ligand/receptor symbols with ens id for now
    species_info = cleanup_species_info(species_info, species_name)
    # Keep relevant columns - use species code, not species_name
    species_mapping = {
        "mmusculus": "mgi_id",
        "rnorvegicus": "rgd_id", 
        "drerio": "zfin_id_id"
    }
 # Fix: Use species (not species_name) for lookup
    species_id = species_mapping.get(species_name, f"{species_name}_homolog_ensembl_gene")
    
    species_info = species_info.dropna(subset=['hgnc_id'])
    species_info = species_info.dropna(subset=[species_id])
    # Merge with ligand data
    ligand_hgnc = [col for col in species_gene_pair1.columns if "Ligand HGNC ID" in col][0]
    species_gene_pair1["Lig HGNC ID"] = species_gene_pair1[ligand_hgnc].apply(extract_hgnc_id)
    species_gene_pair1 = species_gene_pair1.merge(species_info, how='left', 
                               left_on='Lig HGNC ID', right_on='hgnc_id',
                               suffixes=('', '_lig'))

    # Turn Interaction ID into clickable links
    # species_gene_pair1[species_gene_pair1.columns[0]] = species_gene_pair1[species_gene_pair1.columns[0]].apply(
    #     lambda x: f"<a href='https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/database/filter/{x}.html'>{x}</a>"
    # )
        # Rename columns for ligand info
    rename_dict = {}
    if f"{species_name}_homolog_goc_score" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_goc_score"] = f"{species} Ligand GOC score"
    if f"{species_name}_homolog_wga_coverage" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_wga_coverage"] = f"{species} Ligand WGA coverage"   
    if f"{species_name}_homolog_perc_id" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_perc_id"] = f"{species} Ligand % Identity"
    if f"{species_name}_homolog_perc_id_r1" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_perc_id_r1"] = f"{species} Ligand Target % Identity"           
    if f"{species_name}_homolog_orthology_confidence" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_orthology_confidence"] = f"{species} Ligand Orthology Confidence"           
    if f"{species_name}_homolog_associated_gene_name" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_associated_gene_name"] = f"{species} Ligand"
    if "mgi_description" in species_gene_pair1.columns:
        rename_dict["mgi_description"] = f"{species} Ligand Name"
    if "description" in species_gene_pair1.columns:
        rename_dict["description"] = f"{species} Ligand Name"
    if f"{species_name}_homolog_ensembl_gene" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_ensembl_gene"] = f"{species} Ligand Ensembl ID"
    if species_id in species_gene_pair1.columns and species_id != f"{species_name}_homolog_ensembl_gene":
        rename_dict[species_id] = f"{species} Ligand ID"
    species_gene_pair1 = species_gene_pair1.rename(columns=rename_dict)
    # Merge with receptor data
    receptor_hgnc = [col for col in species_gene_pair1.columns if "Receptor HGNC ID" in col][0]
    species_gene_pair1["Rec HGNC ID"] = species_gene_pair1[receptor_hgnc].apply(extract_hgnc_id)
    species_gene_pair1 = species_gene_pair1.merge(species_info, how='left', 
                               left_on='Rec HGNC ID', right_on='hgnc_id',
                               suffixes=('', '_rec'))
    # Rename columns for ligand info
    rename_dict = {}
    if f"{species_name}_homolog_goc_score" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_goc_score"] = f"{species} Receptor GOC score"
    if f"{species_name}_homolog_wga_coverage" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_wga_coverage"] = f"{species} Receptor WGA coverage"   
    if f"{species_name}_homolog_perc_id" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_perc_id"] = f"{species} Receptor % Identity"
    if f"{species_name}_homolog_perc_id_r1" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_perc_id_r1"] = f"{species} Receptor Target % Identity"           
    if f"{species_name}_homolog_orthology_confidence" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_orthology_confidence"] = f"{species} Receptor Orthology Confidence"           
    if f"{species_name}_homolog_associated_gene_name" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_associated_gene_name"] = f"{species} Receptor"
    if "mgi_description" in species_gene_pair1.columns:
        rename_dict["mgi_description"] = f"{species} Receptor Name"
    if "description" in species_gene_pair1.columns:
        rename_dict["description"] = f"{species} Receptor Name"
    if f"{species_name}_homolog_ensembl_gene" in species_gene_pair1.columns:
        rename_dict[f"{species_name}_homolog_ensembl_gene"] = f"{species} Receptor Ensembl ID"
    if species_id in species_gene_pair1.columns and species_id != f"{species_name}_homolog_ensembl_gene":
        rename_dict[species_id] = f"{species} Receptor ID"
        
    species_gene_pair1 = species_gene_pair1.rename(columns=rename_dict)
    # Drop hgnc_id columns
    cols_to_drop = [col for col in species_gene_pair1.columns if col in ['hgnc_id', 'hgnc_id_rec', 'hgnc_symbol_rec',
                                                                         'hgnc_id_lig', 'hgnc_symbol_lig',  'Rec HGNC ID','ensembl_gene_id']]
    if cols_to_drop:
        origDF = species_gene_pair1.drop(columns=cols_to_drop)

    # Drop columns where all values are NaN
    species_gene_pair1 = species_gene_pair1.dropna(axis=1, how='all')
    col_to_rename = [col for col in species_gene_pair1.columns if "Human LR Pair" in col][0]
    species_gene_pair1.rename(columns={col_to_rename: "LR Pair Card"}, inplace=True)
    # Update ligand_col values only if species-specific replacement is available
    # Only override if the cleaned species-specific columns exist now
    new_ligand_col = f"{species} Ligand"
    new_receptor_col = f"{species} Receptor"
    
    if new_ligand_col in species_gene_pair1.columns:
        species_gene_pair1[ligand_col] = [
            new if pd.notna(new) and str(new).strip() else orig
            for new, orig in zip(species_gene_pair1[new_ligand_col], species_gene_pair1[ligand_col])
        ]

    if new_receptor_col in species_gene_pair1.columns:
        species_gene_pair1[receptor_col] = [
            new if pd.notna(new) and str(new).strip() else orig
            for new, orig in zip(species_gene_pair1[new_receptor_col], species_gene_pair1[receptor_col])
        ]

    if species in ["Mouse", "Rat", "Zebrafish"]:
        ligand_mgi_id_col = [col for col in species_gene_pair1.columns if f"Ligand {id_prefix} ID" in col][0]
        receptor_mgi_id_col = [col for col in species_gene_pair1.columns if f"Receptor {id_prefix} ID" in col][0]
    
        ligand_source = species_gene_pair1[f"{species} Ligand ID"]
        receptor_source = species_gene_pair1[f"{species} Receptor ID"]
    
        # URL builders by species
        def build_link(val, species):
            if pd.isna(val) or not str(val).strip() or str(val).strip().lower() in {"none", "nan"}:
                return ""
        
            val = str(val).strip()
            
            try:
                val = str(int(float(val)))  # Remove trailing ".0"
            except ValueError:
                pass
        
            if species == "Mouse":
                return f'<a href="https://www.informatics.jax.org/marker/{val}" target="_blank">{val}</a>'
            elif species == "Rat":
                return f'<a href="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=RGD:{val}" target="_blank">RGD:{val}</a>'
            elif species == "Zebrafish":
                return f'<a href="https://zfin.org/{val}" target="_blank">{val}</a>'
            else:
                return val
    
        ligand_links = [build_link(val, species) for val in ligand_source]
        receptor_links = [build_link(val, species) for val in receptor_source]
    
        # Replace only if new link is non-empty
        species_gene_pair1[ligand_mgi_id_col] = [
            link if link else orig
            for link, orig in zip(ligand_links, species_gene_pair1[ligand_mgi_id_col])
        ]
    
        species_gene_pair1[receptor_mgi_id_col] = [
            link if link else orig
            for link, orig in zip(receptor_links, species_gene_pair1[receptor_mgi_id_col])
        ]
    
        # Drop original columns
        species_gene_pair1 = species_gene_pair1.drop(columns=[f"{species} Ligand ID", f"{species} Receptor ID"])
    else:
        # Update ligand_ens values only if species-specific replacement is available
    
        # Only override if the cleaned species-specific columns exist now
        new_ligand_ens = f"{species} Ligand Ensembl ID"
        new_receptor_ens = f"{species} Receptor Ensembl ID"
        ligand_ens = [col for col in species_gene_pair1.columns if "Ligand Ensembl ID" in col][0]
        receptor_ens= [col for col in species_gene_pair1.columns if "Receptor Ensembl ID" in col][0]
        if new_ligand_ens in species_gene_pair1.columns:
            species_gene_pair1[ligand_ens] = [
                new if pd.notna(new) and str(new).strip() else orig
                for new, orig in zip(species_gene_pair1[new_ligand_ens], species_gene_pair1[ligand_ens])
            ]
    
        if new_receptor_ens in species_gene_pair1.columns:
            species_gene_pair1[receptor_ens] = [
                new if pd.notna(new) and str(new).strip() else orig
                for new, orig in zip(species_gene_pair1[new_receptor_ens], species_gene_pair1[receptor_ens])
            ]
    
    # ligand_col = f"{species} Ligand"   
    # receptor_col = f"{species} Receptor"  
    def format_lr_pair(row):
        if row[ligand_Location] in ['secreted', '']:
            return f"{row[ligand_col]} <span style='font-size: 15px;'>○</span> <span style='font-size: 24px;'>⤚</span> {row[receptor_col]}"
        elif row[receptor_Location] == 'plasma membrane':
            return f"{row[ligand_col]} <span style='font-size: 24px;'>⤙</span> <span style='font-size: 24px;'>⤚</span> {row[receptor_col]}"
        else:
            return f"{row[ligand_col]} {row[receptor_col]}" #\u2192 was the arrow but removed for now
            
    species_gene_pair1.loc[:, f"{species} LR Pair"] = species_gene_pair1.apply(format_lr_pair, axis=1)
       # Identify relevant columns for the species
    species_lr_pair_col = f"{species} LR Pair"
    species_columns = [
        col for col in species_gene_pair1.columns
        if (id_prefix in col or species in col)
        and col not in [ligand_col, receptor_col, species_lr_pair_col]  # <- singular
    ]

    new_order = [human_columns[0]]+ ["LR Pair Card"] + [f"{species} LR Pair", ligand_col, receptor_col] + species_columns +  human_columns[2:]
    species_gene_pair1 = species_gene_pair1[new_order].reset_index(drop=True)
    species_gene_pair1 = species_gene_pair1.loc[:, ~species_gene_pair1.columns.duplicated()]
    species_gene_pair1.columns = [
        col if col == species_lr_pair_col else col.replace(species, "").strip()
        for col in species_gene_pair1.columns
    ]
    if id_prefix == "Ensembl":
        species_gene_pair1 = species_gene_pair1.drop(columns=['Ligand Ensembl ID',
                                                              'Receptor Ensembl ID'])
    species_gene_pair1 = species_gene_pair1.drop(columns=['Ligand',
                                                              'Receptor'])

    return species_gene_pair1


# Process each species
# Mouse
mouse_gene_pair1 = process_species(gene_pair, gene_pair000, "Mouse", "MGI", 1, 1)
MouselrPairsCount = len(pd.unique(mouse_gene_pair1.iloc[:, 0]))# unique
HumanMouseLRPairsPer = (MouselrPairsCount/lrPairsCount)*100
HumanMouseLRPairsPer = round(HumanMouseLRPairsPer, 2)
# Round up to the nearest 0.5%
HumanMouseLRPairsPer = ((HumanMouseLRPairsPer * 2 + 1) // 1) / 2
### adding the mouse-specific annotations ####
mouse_gene_pair1 = mouse_gene_pair1.reset_index(drop=True) 
mouse_rat_info = pd.read_csv("data/mouse_name_mapping.csv")
ligand_mgi_id_col = [col for col in mouse_gene_pair1.columns if f"Ligand MGI ID" in col][0]
receptor_mgi_id_col = [col for col in mouse_gene_pair1.columns if f"Receptor MGI ID" in col][0]   
mapping_mouse_name = dict(zip(mouse_rat_info['MGI ID'], mouse_rat_info['MGI description']))
mapping_mouse_ens = dict(zip(mouse_rat_info['MGI ID'], mouse_rat_info['Gene stable ID']))
# extract mgi
def extract_mgi_id(col):
    """Use regular expression to extract the HGNC ID after 'HGNC:'."""
    match = re.search(r'MGI:(\d+)', col)
    if match:
        return 'MGI:' +str(match.group(1))
    return None
    

mouse_gene_pair1['Ligand MGI ID'] = mouse_gene_pair1[ligand_mgi_id_col].apply(extract_mgi_id)
mouse_gene_pair1['Receptor MGI ID'] = mouse_gene_pair1[receptor_mgi_id_col].apply(extract_mgi_id)
# Apply the mapping to 'Ligand Name'
mouse_gene_pair1['Ligand Name'] = mouse_gene_pair1.apply(
    lambda row: mapping_mouse_name.get(row['Ligand MGI ID'], row['Ligand Name'])
    if pd.isna(row['Ligand Name']) or str(row['Ligand Name']).strip() == '' else row['Ligand Name'],
    axis=1
)

mouse_gene_pair1["Ligand Name"] = [
    f'<span title="{ligand_name}">{ligand_symbol}</span>'
    for ligand_name, ligand_symbol in zip(mouse_gene_pair1["Ligand Name"], 
                                              mouse_gene_pair1["Ligand Name"])
]

 # Apply the mapping to 'Receptor Name'
mouse_gene_pair1['Receptor Name'] = mouse_gene_pair1.apply(
    lambda row: mapping_mouse_name.get(row['Receptor MGI ID'], row['Receptor Name'])
    if pd.isna(row['Receptor Name']) or str(row['Receptor Name']).strip() == '' else row['Receptor Name'],
    axis=1
)

mouse_gene_pair1["Receptor Name"] = [
    f'<span title="{receptor_name}">{receptor_symbol}</span>'
    for receptor_name, receptor_symbol in zip(mouse_gene_pair1["Receptor Name"], 
                                              mouse_gene_pair1["Receptor Name"])
]

# Apply the mapping to 'Ligand Ensembl ID'
mouse_gene_pair1['Ligand Ensembl ID'] = mouse_gene_pair1.apply(
    lambda row: mapping_mouse_ens.get(row['Ligand MGI ID'], row['Ligand Ensembl ID'])
    if pd.isna(row['Ligand Ensembl ID']) or str(row['Ligand Ensembl ID']).strip() == '' else row['Ligand Ensembl ID'],
    axis=1
)

 # Apply the mapping to 'Receptor Ensembl ID'
mouse_gene_pair1['Receptor Ensembl ID'] = mouse_gene_pair1.apply(
    lambda row: mapping_mouse_ens.get(row['Receptor MGI ID'], row['Receptor Ensembl ID'])
    if pd.isna(row['Receptor Ensembl ID']) or str(row['Receptor Ensembl ID']).strip() == '' else row['Receptor Ensembl ID'],
    axis=1
)
mouse_gene_pair1 = mouse_gene_pair1.drop(columns=['Ligand MGI ID',
                                                  'Receptor MGI ID'])
################################################################

# Rat
rat_gene_pair1 = process_species(gene_pair, gene_pair000, "Rat", "RGD", 2, 2)
# Zebrafish
zebrafish_gene_pair1 = process_species(gene_pair, gene_pair000, "Zebrafish", "ZFIN", 3, 3)

# Chimpanzee
chimpanzee_gene_pair1 = process_species(gene_pair, gene_pair000, "Chimpanzee", "Ensembl", 4, 4)

# Chicken
chicken_gene_pair1 = process_species(gene_pair, gene_pair000, "Chicken", "Ensembl", 4, 4)

# Pig
pig_gene_pair1 = process_species(gene_pair, gene_pair000, "Pig", "Ensembl", 4, 4)

# Cow
cow_gene_pair1 = process_species(gene_pair, gene_pair000, "Cow", "Ensembl", 4, 4)

# Dog
dog_gene_pair1 = process_species(gene_pair, gene_pair000, "Dog", "Ensembl", 4, 4)

# Horse
horse_gene_pair1 = process_species(gene_pair, gene_pair000, "Horse", "Ensembl",4, 4)

# Sheep
sheep_gene_pair1 = process_species(gene_pair, gene_pair000, "Sheep", "Ensembl", 4, 4)

# Marmoset
marmoset_gene_pair1 = process_species(gene_pair, gene_pair000, "Marmoset", "Ensembl", 4, 4)

# Rhesus macaque
macaque_gene_pair1 = process_species(gene_pair, gene_pair000, "Macaque", "Ensembl", 4, 4)

