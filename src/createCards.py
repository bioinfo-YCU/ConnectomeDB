## Function to create Ligand-Receptor pair cards

import os
import jinja2
import sys
import pandas as pd
import numpy as np

sys.path.append(os.path.abspath("src"))  
import fetchGSheet
from createDataTable import gene_pair0

# Paths
TEMPLATE_PATH = 'HTML/cardTemplate.html'
OUTPUT_DIR = 'data/cards/'

def load_template(template_path):
    """Load Jinja2 template from a file."""
    with open(template_path, 'r') as file:
        return jinja2.Template(file.read())

# Function to convert the links
def convert_hgnc_url(col):
   # Extract the HGNC ID from the original URL
    hgnc_id = col.split("HGNC:")[1].split('"')[0]  # Extract the ID number (e.g., "31702")
    # Extract the visible text inside the <a> tag
    visible_text = col.split(">")[1].split("<")[0]  # Extract "HGNC:31702"
    # Construct the new link, keeping the text intact
    new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}" target="_blank">{visible_text}</a>'
    return new_link


def convert_hgnc_url_disease(col):
   # Extract the HGNC ID from the original URL
    hgnc_id = col.split("HGNC:")[1].split('"')[0]  # Extract the ID number (e.g., "31702")
    # Extract the visible text inside the <a> tag
    visible_text = col.split(">")[1].split("<")[0]  # Extract "HGNC:31702"
    # Construct the new link, keeping the text intact
    new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#diseases" target="_blank">{visible_text}</a>'
    return new_link
    
def convert_hgnc_url_exp(col):
   # Extract the HGNC ID from the original URL
    hgnc_id = col.split("HGNC:")[1].split('"')[0]  # Extract the ID number (e.g., "31702")
    # Extract the visible text inside the <a> tag
    visible_text = col.split(">")[1].split("<")[0]  # Extract "HGNC:31702"
    # Construct the new link, keeping the text intact
    new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#expression" target="_blank">{visible_text}</a>'
    return new_link

def prepare_dataframes(gene_pair0):
    """Prepare interaction, ligand, and receptor dataframes."""
    DBlength = len(gene_pair0)
    gene_pair0["Interaction ID"] = [f"CDB{str(i).zfill(4)}" for i in range(1, DBlength + 1)]
    gene_pair0["Interaction Type"] = [
        f'{ligandLocation} {ligand} ligand binds to {receptorLocation} {receptor} receptor'
        for ligand, ligandLocation, receptor, receptorLocation in zip(
            gene_pair0["Ligand"], gene_pair0["Ligand location"],
            gene_pair0["Receptor"], gene_pair0["Receptor location"]
        )
    ]
    interaction_card = gene_pair0[["Interaction ID", "Human LR Pair", "Interaction Type", "Perplexity", "PMID support"]]
    interaction_card["Perplexity"] = interaction_card["Perplexity"].str.replace('size=30', 'size=80')

    pop_up_info_lim = fetchGSheet.pop_up_info[
        ["Approved symbol", "Alias symbol", "Previous symbol", "Date symbol changed"]
    ].drop_duplicates(subset="Approved symbol", keep="first")
    
    ligand_card = gene_pair0[["Human LR Pair", "Ligand", "Ligand name", "Ligand HGNC ID", "Ligand MGI ID", "Ligand RGD ID", "Ligand location"]].merge(
        pop_up_info_lim, how='left', left_on='Ligand', right_on='Approved symbol'
    ).drop_duplicates(subset='Human LR Pair', keep="first").drop(columns=["Ligand", "Approved symbol"])

    ligand_card_1 = ligand_card[["Human LR Pair", "Alias symbol", "Date symbol changed", "Ligand name"]] 
    ligand_card_2 = ligand_card[["Human LR Pair", "Ligand HGNC ID", "Ligand location"]] 
    # convert links
    ligand_card_2["HGNC gene card"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url)
    ligand_card_2["Disease relevance"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url_disease)
    ligand_card_2["Expression Profile"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url_exp)
    ligand_card_2 = ligand_card_2[["Human LR Pair", "Ligand HGNC ID", "HGNC gene card", "Disease relevance", "Expression Profile", "Ligand location"]]       

          
    receptor_card = gene_pair0[["Human LR Pair", "Receptor", "Receptor name", "Receptor HGNC ID", "Receptor MGI ID", "Receptor RGD ID", "Receptor location"]].merge(
        pop_up_info_lim, how='left', left_on='Receptor', right_on='Approved symbol'
    ).drop_duplicates(subset='Human LR Pair', keep="first").drop(columns=["Receptor", "Approved symbol"])
    
    receptor_card_1 = receptor_card[["Human LR Pair", "Alias symbol", "Date symbol changed", "Receptor name"]] 
    receptor_card_2 = receptor_card[["Human LR Pair", "Receptor HGNC ID", "Receptor location"]] 
    receptor_card_2["HGNC gene card"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url)
    receptor_card_2["Disease relevance"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_disease)
    receptor_card_2["Expression Profile"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_exp)
    receptor_card_2 = receptor_card_2[["Human LR Pair", "Receptor HGNC ID", "HGNC gene card", "Disease relevance", "Expression Profile", "Receptor location"]]       

    return interaction_card, ligand_card_1, ligand_card_2,receptor_card_1, receptor_card_2

def generate_html_files(template, interaction_card, ligand_card_1, receptor_card_1, ligand_card_2, receptor_card_2, output_dir):
    """Generate HTML files for each Human LR Pair."""
    column_values = interaction_card["Human LR Pair"].dropna().unique()
    os.makedirs(output_dir, exist_ok=True)
    
    for value in column_values:
        value1, value2 = value.split()
        row0 = interaction_card[interaction_card['Human LR Pair'] == value]
        row1 = ligand_card_1[ligand_card_1['Human LR Pair'] == value]
        row2 = receptor_card_1[receptor_card_1['Human LR Pair'] == value]
        row3 = ligand_card_2[ligand_card_2['Human LR Pair'] == value]
        row4 = receptor_card_2[receptor_card_2['Human LR Pair'] == value]
        

        table0_data = row0.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row0.empty else {}
        table1_data = row1.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row1.empty else {}
        table2_data = row2.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row2.empty else {}
        table3_data = row3.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row3.empty else {}
        table4_data = row4.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row4.empty else {}

        rendered_content = template.render(
            value1=value1,
            value2=value2,
            table0_data=table0_data,
            table1_data=table1_data,
            table2_data=table2_data,
            table3_data=table3_data,
            table4_data=table4_data
        )
        
        output_file = os.path.join(output_dir, f"{value1} {value2}.html")
        with open(output_file, 'w') as file:
            file.write(rendered_content)

# Main execution
if __name__ == "__main__":
    template = load_template(TEMPLATE_PATH)
    interaction_card, ligand_card_1, receptor_card_1, ligand_card_2, receptor_card_2 = prepare_dataframes(gene_pair0)
    generate_html_files(template, interaction_card, ligand_card_1, receptor_card_1, ligand_card_2, receptor_card_2, OUTPUT_DIR)
