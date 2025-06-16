import sys
import os
import pandas as pd
import yaml
import re

# Assuming createDataTable_perSpecies module is in your path or can be imported.
sys.path.append(os.path.abspath("src"))
import createDataTable_perSpecies
from createDataTable import human_gene_pair

def update_connectomedb_qmd(qmd_file_path: str, lr_pair_data: list, species_name: str, species: str, ortholog):
    """
    Updates the YAML front matter of a Quarto Markdown (.qmd) file
    with dynamic LR pair count and species name.

    Args:
        qmd_file_path (str): The path to the .qmd file to be updated.
        lr_pair_data (list): A list containing the LR pair data (e.g.,
                             createDataTable_perSpecies.mouse_gene_pair1["Mouse LR Pair"]).
                             The length of this list is used for the count.
        species_name (str): The common name of the species (e.g., "Mus musculus").
    """
    print(f"--- Updating {qmd_file_path} for {species_name} ---")

    # --- 1. Compute dynamic values
    lr_pairs_count = len(lr_pair_data)
    formatted_count = f"{lr_pairs_count:,}"

    if ortholog:
        Ortholog="Ortholog"
    else:
        Ortholog=""

    # --- 2. Build YAML contents
    header_script = (
        "<script src='../js/keepDropdownMenuGold.js'></script>"
        if species in ["Human", "Mouse"]
        else "<script src='../../js/keepDropdownMenuGold.js'></script>"
    )
    
    yaml_data = {
        "title": f"{{{{< fa database >}}}} ConnectomeDB2025: {species} – *{species_name}* </span> {Ortholog} LR Pairs",
        "execute": {"echo": False},
        "format": {
            "html": {"table": False}
        },
        "header-includes": header_script
    }
    
    yaml_block = "---\n" + yaml.dump(yaml_data, sort_keys=False) + "---"
    
    # Read the QMD template
    try:
        with open(qmd_file_path, 'r', encoding='utf-8') as f:
            template = f.read()
    except FileNotFoundError:
        print(f"Error: The file '{qmd_file_path}' was not found.")
        return
    except Exception as e:
        print(f"Error reading file '{qmd_file_path}': {e}")
        return
    
    # Use a lambda to safely replace the YAML block without triggering backslash escapes
    new_qmd = re.sub(r"(?s)^---.*?---", lambda m: yaml_block, template, count=1)


    # --- 4. Write out and render
    try:
        with open(qmd_file_path, 'w', encoding='utf-8') as f:
            f.write(new_qmd)
        print(f"Successfully updated '{qmd_file_path}' for {species_name}.")
    except Exception as e:
        print(f"Error writing to file '{qmd_file_path}': {e}")

# --- Example Usage ---

# The order should be (qmd_file_path, lr_pair_data, species_name,species,ortholog= True)
#  'Human LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/human.qmd",
    lr_pair_data=human_gene_pair[human_gene_pair.columns[0]],
    species_name="Homo sapiens",
    species = "Human",
    ortholog = False
)

#  'Mouse LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/mouseOrth.qmd",
    lr_pair_data=createDataTable_perSpecies.mouse_gene_pair1["Mouse LR Pair"],
    species_name="Mus musculus",
    species = "Mouse",
    ortholog = True
)

#  'Rat LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/ratOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.rat_gene_pair1["Rat LR Pair"],
    species_name="Rattus norvegicus",
    species = "Rat",
    ortholog = True
)

#  'Chimp LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/chimpOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.chimpanzee_gene_pair1["Chimpanzee LR Pair"],
    species_name="Pan troglodytes",
    species = "Chimp",
    ortholog = True
)

#  'Macaque LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/macaqueOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.macaque_gene_pair1["Macaque LR Pair"],
    species_name="Macaca mulatta",
    species = "Macaque",
    ortholog = True
)

#  'Marmoset LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/marmosetOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.marmoset_gene_pair1["Marmoset LR Pair"],
    species_name="Callithrix jacchus",
    species = "Marmoset",
    ortholog = True
)

#  'Sheep LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/sheepOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.sheep_gene_pair1["Sheep LR Pair"],
    species_name="Ovis aries",
    species = "Sheep",
    ortholog = True
)

#  'Horse LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/horseOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.horse_gene_pair1["Horse LR Pair"],
    species_name="Equus caballus",
    species = "Horse",
    ortholog = True
)

#  'Dog LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/dogOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.dog_gene_pair1["Dog LR Pair"],
    species_name="Canis lupus familiaris",
    species = "Dog",
    ortholog = True
)

#  'Cow LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/cowOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.cow_gene_pair1["Cow LR Pair"],
    species_name="Bos taurus",
    species = "Cow",
    ortholog = True
)

#  'Pig LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/pigOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.pig_gene_pair1["Pig LR Pair"],
    species_name="Sus scrofa",
    species = "Pig",
    ortholog = True
)

#  'Zebrafish LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/zebrafishOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.zebrafish_gene_pair1["Zebrafish LR Pair"],
    species_name="Danio rerio",
    species = "Zebrafish",
    ortholog = True
)

#  'Chicken LR Pair' data
update_connectomedb_qmd(
    qmd_file_path="database/other/chickenOrth.qmd", 
    lr_pair_data=createDataTable_perSpecies.chicken_gene_pair1["Chicken LR Pair"],
    species_name="Gallus gallus",
    species = "Chicken",
    ortholog = True
)


