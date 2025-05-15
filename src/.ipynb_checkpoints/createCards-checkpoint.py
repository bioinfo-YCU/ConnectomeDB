## Function to create Ligand-Receptor pair cards

import os
import jinja2
import sys
import pandas as pd
import numpy as np
import time
import base64
import re


sys.path.append(os.path.abspath("src"))  
import fetchGSheet
from createDataTable import pop_up_info, gene_pair0

# Paths
TEMPLATE_PATH = 'HTML/cardTemplate.html'
OUTPUT_DIR = 'data/cards/'
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def load_template(template_path):
    """Load Jinja2 template from a file."""
    with open(template_path, 'r') as file:
        return jinja2.Template(file.read())

def encode_image(image_path):
    """Encode an image to base64."""
    try:
        with open(image_path, "rb") as image_file:
            return base64.b64encode(image_file.read()).decode('utf-8')
    except FileNotFoundError:
        return None

# Function to extract the HGNC ID from the anchor tag URL
def extract_hgnc_id(col):
    # Use regular expression to extract the HGNC ID after "HGNC:"
    match = re.search(r'HGNC:(\d+)', col)
    if match:
        return match.group(1)  # Return the HGNC ID (number part)
    return None  # Return None if the format doesn't match or it's not a string

# Updated functions to convert the HGNC link (using extract_hgnc_id)
def convert_hgnc_url(col):
    hgnc_id = extract_hgnc_id(col)  # Extract the HGNC ID
    if hgnc_id:
        visible_text = "genecard.org"
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}" target="_blank">{visible_text}</a>'
        return new_link
    return None

def convert_hgnc_url_disease(col):
    hgnc_id = extract_hgnc_id(col)  # Extract the HGNC ID
    if hgnc_id:
        visible_text = "see here"
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#diseases" target="_blank">{visible_text}</a>'
        return new_link
    return None

def convert_hgnc_url_exp(col):
    hgnc_id = extract_hgnc_id(col)  # Extract the HGNC ID
    if hgnc_id:
        visible_text = "see here"
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#expression" target="_blank">{visible_text}</a>'
        return new_link
    return None

def prepare_dataframes(gene_pair0):
    """Prepare interaction, ligand, and receptor dataframes."""
    # DBlength = len(gene_pair0)
    # gene_pair0["Interaction ID"] = [f"CDB{str(i).zfill(4)}" for i in range(1, DBlength + 1)]
    gene_pair0["Interaction Type"] = [
        f'{ligand} {ligandLocation} ligand binds to {receptor} {receptorLocation} receptor'
        for ligand, ligandLocation, receptor, receptorLocation in zip(
            gene_pair0["Ligand"], gene_pair0["Ligand location"],
            gene_pair0["Receptor"], gene_pair0["Receptor location"]
        )
    ]
    interaction_card = gene_pair0[["Interaction ID", "Human LR Pair", "Interaction Type", "Perplexity", "PMID", "KEGG Pathway", "Cancer-related", "Disease Type"]]
    interaction_card["Perplexity"] = interaction_card["Perplexity"].str.replace('size=30', 'size=80')

    pop_up_info_lim = pop_up_info[
        ["Approved symbol", "Alias symbol", "Previous symbol", "Date symbol changed"]
    ].drop_duplicates(subset="Approved symbol", keep="first")
    
    ligand_card = gene_pair0[["Human LR Pair", "Ligand", "Ligand name", "Ligand HGNC ID", "Ligand MGI ID", "Ligand RGD ID", "Ligand location"]].merge(
        pop_up_info_lim, how='left', left_on='Ligand', right_on='Approved symbol'
    ).drop_duplicates(subset='Human LR Pair', keep="first").drop(columns=["Ligand", "Approved symbol"])

    ligand_card_1 = ligand_card[["Human LR Pair", "Alias symbol", "Date symbol changed", "Ligand name"]] 
    ligand_card_2 = ligand_card[["Human LR Pair", "Ligand HGNC ID", "Ligand location"]] 
    # Convert links
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

    return interaction_card, ligand_card_1, ligand_card_2, receptor_card_1, receptor_card_2

def generate_html_files(template, interaction_card, ligand_card_1, receptor_card_1, ligand_card_2, receptor_card_2, output_dir):
    """Generate HTML files for each Human LR Pair."""
    column_values = interaction_card["Human LR Pair"].dropna().unique()
    os.makedirs(output_dir, exist_ok=True)

    # Encode the plotlegend image to base64
    plotlegend_image_path = "data/image/plotlegend.webp"
    plotlegend_base64 = encode_image(plotlegend_image_path)  # Convert WebP to base64

    for value in column_values:
        value1, value2 = value.split()
        row0 = interaction_card[interaction_card['Human LR Pair'] == value]
        row1 = ligand_card_1[ligand_card_1['Human LR Pair'] == value]
        row2 = receptor_card_1[receptor_card_1['Human LR Pair'] == value]
        row3 = ligand_card_2[ligand_card_2['Human LR Pair'] == value]
        row4 = receptor_card_2[receptor_card_2['Human LR Pair'] == value]

        # Check if the HTML files exist
        ligand_image_path = f'data/gene_expr_plots/{value1}.html'
        receptor_image_path = f'data/gene_expr_plots/{value2}.html'
        
        if os.path.exists(ligand_image_path):
            with open(ligand_image_path, "r") as html_file:
                ligand_image = html_file.read()  # Read the HTML content
        else:
            ligand_image = "Plot does not exist"
        
        if os.path.exists(receptor_image_path):
            with open(receptor_image_path, "r") as html_file:
                receptor_image = html_file.read()  # Read the HTML content
        else:
            receptor_image = "Plot does not exist"


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
            table4_data=table4_data,
            ligand_image=ligand_image,
            receptor_image=receptor_image,
            plotlegend_base64=plotlegend_base64 
        )
        
        output_file = os.path.join(output_dir, f"{value1} {value2}.html")
        with open(output_file, 'w') as file:
            #time.sleep(0.5)
            file.write(rendered_content)

# Main execution
if __name__ == "__main__":
    template = load_template(TEMPLATE_PATH)
    # if only one replace gene_pair0 to e.g. gene_pair0[gene_pair0["Human LR Pair"] == "PLG IGF2R"]
    interaction_card, ligand_card_1, receptor_card_1, ligand_card_2, receptor_card_2 = prepare_dataframes(gene_pair0)
    generate_html_files(template, interaction_card, ligand_card_1, receptor_card_1, ligand_card_2, receptor_card_2, OUTPUT_DIR)
