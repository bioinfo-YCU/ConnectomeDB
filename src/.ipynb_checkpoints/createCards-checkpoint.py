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
from createDataTable import pop_up_info, gene_pair0, generate_perplexity_links
from createFunctionalAnnotTable import gene_pair_annot_ligand, gene_pair_annot_receptor


# Add Disease (specific) to cards
df= pd.read_csv("data/disease_annotations_per_pair.csv")
df = df.groupby('interaction')['disease'].apply(', '.join).reset_index()
# Map and set to unknown for now
mapping = dict(zip(df['interaction'],df['disease']))
gene_pair0["Disease"] = gene_pair0['Human LR Pair'].map(mapping).fillna("unknown")

gene_pair0 = generate_perplexity_links(
    gene_pair0,
    pathway_col="Disease",
    default_query_template="What-diseases-is-the-ligand-receptor-pair-{pair}-associated-with"
)

# if only one replace gene_pair0 to e.g. 
gene_pair_input = gene_pair0[gene_pair0["Human LR Pair"].isin(["VEGFA KDR", "ADAM17 IL6R"])]
#gene_pair_input = gene_pair0 


# add Ligand/Receptor group info
agg_func = lambda x: ', '.join(sorted(set(map(str, x))))
gene_pair_annot_ligand = gene_pair_annot_ligand.groupby('Ligand HGNC ID').agg(agg_func).reset_index()
ligand_mapping = dict(zip(gene_pair_annot_ligand['Ligand HGNC ID'],gene_pair_annot_ligand['Ligand group']))
gene_pair_annot_receptor = gene_pair_annot_receptor.groupby('Receptor HGNC ID').agg(agg_func).reset_index()
receptor_mapping = dict(zip(gene_pair_annot_receptor['Receptor HGNC ID'],gene_pair_annot_receptor['Receptor group']))

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
        visible_text = "GeneCards"
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}" target="_blank">{visible_text}</a>'
        return new_link
    return None

def convert_hgnc_url_disease(col):
    hgnc_id = extract_hgnc_id(col)  # Extract the HGNC ID
    if hgnc_id:
        visible_text = "MalaCards"
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#diseases" target="_blank">{visible_text}</a>'
        return new_link
    return None

def convert_hgnc_url_exp(col):
    hgnc_id = extract_hgnc_id(col)  # Extract the HGNC ID
    if hgnc_id:
        visible_text = "mRNA expression in normal human tissues"
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#expression" target="_blank">{visible_text}</a>'
        return new_link
    return None

def prepare_dataframes(gene_pair_input):
    """Prepare interaction, ligand, and receptor dataframes."""
    # DBlength = len(gene_pair_input)
    # gene_pair_input["Interaction ID"] = [f"CDB{str(i).zfill(4)}" for i in range(1, DBlength + 1)]
    gene_pair_input["Interaction Type"] = [
        f'{ligand} {ligandLocation} ligand binds to {receptor} {receptorLocation} receptor'
        for ligand, ligandLocation, receptor, receptorLocation in zip(
            gene_pair_input["Ligand"], gene_pair_input["Ligand Location"],
            gene_pair_input["Receptor"], gene_pair_input["Receptor Location"]
        )
    ]
    interaction_card = gene_pair_input[["Interaction ID", "Human LR Pair", "Interaction Type", "Perplexity", "PMID", "KEGG Pathway",  "PROGENy Pathway", "Cancer-related", "Disease Type", "Disease"]]
    interaction_card["Perplexity"] = interaction_card["Perplexity"].str.replace('size=30', 'size=80')

    pop_up_info_lim = pop_up_info[
        ["Approved symbol", "Alias symbol", "Previous symbol", "Date symbol changed"]
    ].drop_duplicates(subset="Approved symbol", keep="first")

    def format_symbol_aliases(old_symbol, aliases):
        # Filter out "N/A" values
        parts = [p for p in (old_symbol, aliases) if p != "N/A"]
        # Return just the symbol if no valid aliases or old symbols
        return f"{', '.join(parts)}" if parts else aliases

    pop_up_info_lim['Other Symbols'] =pop_up_info_lim.apply(
        lambda row: format_symbol_aliases(row["Previous symbol"], row["Alias symbol"]),
        axis=1
    )
    
    ligand_card = gene_pair_input[["Human LR Pair", "Ligand", "Ligand name", "Ligand HGNC ID", "Ligand MGI ID", "Ligand RGD ID", "Ligand Location"]].merge(
        pop_up_info_lim, how='left', left_on='Ligand', right_on='Approved symbol'
    ).drop_duplicates(subset='Human LR Pair', keep="first").drop(columns=["Ligand", "Approved symbol"])

    ligand_card_1 = ligand_card[["Human LR Pair", "Other Symbols", "Ligand name"]] 
    ligand_card_2 = ligand_card[["Human LR Pair", "Ligand HGNC ID", "Ligand Location"]] 
    # Convert links
    ligand_card_2["HGNC gene card"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url)
    ligand_card_2["Disease relevance"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url_disease)
    ligand_card_2["Expression Profile"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url_exp)
    # Add ligand group
    ligand_card_2["Lineage group"] = ligand_card_2['Ligand HGNC ID'].map(ligand_mapping).fillna("none")
    # Add extra line to balance out cards if necessary
    def add_spacing(text):
        if len(text) <= 60:
            return f"{text}<br><br>"  # extra line for short text
        else:
            return f"{text}"      # no extra line for long text


    #ligand_card_2["Lineage group"] = ligand_card_2["Lineage group"].apply(add_spacing)
    ligand_card_2 = ligand_card_2[["Human LR Pair", "Ligand HGNC ID", "HGNC gene card", "Lineage group",  "Ligand Location", "Disease relevance", "Expression Profile"]]       

    
    receptor_card = gene_pair_input[["Human LR Pair", "Receptor", "Receptor name", "Receptor HGNC ID", "Receptor MGI ID", "Receptor RGD ID", "Receptor Location"]].merge(
        pop_up_info_lim, how='left', left_on='Receptor', right_on='Approved symbol'
    ).drop_duplicates(subset='Human LR Pair', keep="first").drop(columns=["Receptor", "Approved symbol"])
    
    receptor_card_1 = receptor_card[["Human LR Pair", "Other Symbols", "Receptor name"]] 
    receptor_card_2 = receptor_card[["Human LR Pair", "Receptor HGNC ID", "Receptor Location"]] 
    receptor_card_2["HGNC gene card"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url)
    receptor_card_2["Disease relevance"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_disease)
    receptor_card_2["Expression Profile"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_exp)
    # Add Receptor group
    receptor_card_2["Lineage group"] = receptor_card_2['Receptor HGNC ID'].map(receptor_mapping).fillna("none")
    # Add extra line to balance out cards if necessary
    #receptor_card_2["Lineage group"] = receptor_card_2["Lineage group"].apply(add_spacing)
    receptor_card_2 = receptor_card_2[["Human LR Pair", "Receptor HGNC ID",  "HGNC gene card", "Lineage group", "Receptor Location","Disease relevance", "Expression Profile" ]]       

    return interaction_card, ligand_card_1, ligand_card_2, receptor_card_1, receptor_card_2

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

        # make the pair cards here
        def convert_pair_url(df_pairs):
            df_pairs["Human LR Pair"] = [
                f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/{lrpair}.html" target="_blank" '
                f'role="button" title="Open {lrpair} card" class="btn btn-outline-primary" '
                f'style="background-color: #3498db; color: white; border-color: #2980b9; font-size: 14px; '
                f'padding: 2px 3px; margin: 2px; text-decoration: none; border-radius: 2px;">'
                f'{lrpair}</a>'
                if pd.notna(lrpair) and lrpair.strip() else ""
                for lrpair in df_pairs["Human LR Pair"]
            ]
            return df_pairs

        # Add other ligand receptor pairs (ligand card)
        ligand_pairs = gene_pair0[gene_pair0['Ligand'] == value1]
        # print out all Human LR Pair values except for value
        ligand_pairs = ligand_pairs[ligand_pairs["Human LR Pair"] != value]
        ligand_pairs = ligand_pairs[["Human LR Pair"]]
        # Apply the transformation directly (no axis=1!)
        ligand_pairs = convert_pair_url(ligand_pairs)
        # Aggregate into one value separated by space
        ligand_pairs = ' '.join([btn for btn in ligand_pairs["Human LR Pair"] if btn])

        # Add other ligand receptor pairs (receptor card)
        receptor_pairs = gene_pair0[gene_pair0['Receptor'] == value2]
        # print out all Human LR Pair values except for value
        receptor_pairs = receptor_pairs[receptor_pairs["Human LR Pair"] != value]
        receptor_pairs = receptor_pairs[["Human LR Pair"]]
        # Apply the transformation directly (no axis=1!)
        receptor_pairs = convert_pair_url(receptor_pairs)
        # Aggregate into one value separated by space
        receptor_pairs = ' '.join([btn for btn in receptor_pairs["Human LR Pair"] if btn])

        # Check if the HTML files exist
        ligand_image_path = f'data/tabula_sapiens/heatmap/{value1}.html'
        receptor_image_path = f'data/tabula_sapiens/heatmap/{value2}.html'
        
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
            ligand_pairs = ligand_pairs,
            receptor_pairs = receptor_pairs
            #plotlegend_base64=plotlegend_base64 
        )
        
        output_file = os.path.join(output_dir, f"{value1} {value2}.html")
        with open(output_file, 'w') as file:
            #time.sleep(0.5)
            file.write(rendered_content)

# Main execution
if __name__ == "__main__":
    template = load_template(TEMPLATE_PATH)
    interaction_card, ligand_card_1, receptor_card_1, ligand_card_2, receptor_card_2 = prepare_dataframes(gene_pair_input)
    generate_html_files(template, interaction_card, ligand_card_1, receptor_card_1, ligand_card_2, receptor_card_2, OUTPUT_DIR)
