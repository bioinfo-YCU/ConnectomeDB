## Function to create Ligand-Receptor pair cards with PMID details on top

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
from createDataTable import pop_up_info, gene_pair0, generate_perplexity_links, gene_pair00
from createFunctionalAnnotTable import gene_pair_annot_ligand, gene_pair_annot_receptor

# Paths
TEMPLATE_PATH = 'HTML/cardwithPMIDTemplate.html'
OUTPUT_DIR = 'data/cards/'
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


#### Bring in PMID details ####

# Load PubMed data
pubmed_data = pd.read_csv("data/pubmed_results.csv")
pubmed_data["Year"] = pubmed_data["Year"].astype(str).str.replace(".0", 
                                                                  "", 
                                                                  regex=False).astype(int)

pubmed_data["PMID"] = pubmed_data["PMID"].astype(str)

# add llm results
bio_keywords = pd.read_csv("data/llm_results.csv")


# Replace spaces in "Human LR Pair" with a placeholder
gene_pair00["Human LR Pair"] = gene_pair00["Human LR Pair"].str.replace(" ", "——")

gene_pair000 = gene_pair00.merge(bio_keywords, how='left', left_on="Human LR Pair", right_on='Human LR Pair')
gene_pair000["Relevance Keywords"] = gene_pair000["Relevance Keywords"].astype(str)
gene_pair000["Human LR Pair"]  = gene_pair000["Human LR Pair"].astype(str)
pubmed_data = pubmed_data.reset_index(drop=True)  # Remove the index


#################################################################################################
#### Bring in original card info ####
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

# add external link icon

icon_html = '<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left:4px;"></i></a>'
columns_to_update = [
    "KEGG Pathway", "PROGENy Pathway", "Cancer-related",
    "Disease Type", "Disease"
]

for col in columns_to_update:
    gene_pair_input[col] = gene_pair_input[col].str.replace(
        "</a>", icon_html, regex=False
    )

# add Ligand/Receptor group info
agg_func = lambda x: ', '.join(sorted(set(map(str, x))))
gene_pair_annot_ligand = gene_pair_annot_ligand.groupby('Ligand HGNC ID').agg(agg_func).reset_index()
ligand_mapping = dict(zip(gene_pair_annot_ligand['Ligand HGNC ID'],gene_pair_annot_ligand['Ligand group']))
gene_pair_annot_receptor = gene_pair_annot_receptor.groupby('Receptor HGNC ID').agg(agg_func).reset_index()
receptor_mapping = dict(zip(gene_pair_annot_receptor['Receptor HGNC ID'],gene_pair_annot_receptor['Receptor group']))


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
        visible_text = 'GeneCards <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>'
        new_link = (
            f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}" '
            f'target="_blank" style="color: #0000EE; text-decoration: underline;">{visible_text}</a>'
        )        
        return new_link
    return None

def convert_hgnc_url_disease(col):
    hgnc_id = extract_hgnc_id(col)  # Extract the HGNC ID
    if hgnc_id:
        visible_text = 'MalaCards <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>'
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#diseases" target="_blank">{visible_text}</a>'
        return new_link
    return None

def convert_hgnc_url_exp(col):
    hgnc_id = extract_hgnc_id(col)  # Extract the HGNC ID
    if hgnc_id:
        visible_text = 'mRNA expression in normal human tissues <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>'
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

    ligand_card_1 = ligand_card[["Human LR Pair", "Ligand name", "Other Symbols" ]] 
    ligand_card_2 = ligand_card[["Human LR Pair", "Ligand HGNC ID", "Ligand Location"]] 
    # Convert links
    ligand_card_2["HGNC gene card"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url)
    ligand_card_2["Disease relevance"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url_disease)
    ligand_card_2["Expression Profile"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url_exp)
    # Add ligand group
    ligand_card_2["Lineage group"] = ligand_card_2['Ligand HGNC ID'].map(ligand_mapping).fillna("none")
    icon_html = '<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left:4px;"></i></a>'
    columns_to_update = ["Ligand HGNC ID"]
    
    for col in columns_to_update:
        ligand_card_2[col] = ligand_card_2[col].str.replace(
            "</a>", icon_html, regex=False
        )
    # Add extra line to balance out cards if necessary
    def add_spacing(text):
        if len(text) <= 60:
            return f"{text}<br><br>"  # extra line for short text
        else:
            return f"{text}"      # no extra line for long text


    #ligand_card_2["Lineage group"] = ligand_card_2["Lineage group"].apply(add_spacing)
    ligand_card_2 = ligand_card_2[["Human LR Pair", "Ligand HGNC ID", "HGNC gene card", "Ligand Location", "Lineage group", "Disease relevance", "Expression Profile"]]       

    
    receptor_card = gene_pair_input[["Human LR Pair", "Receptor", "Receptor name", "Receptor HGNC ID", "Receptor MGI ID", "Receptor RGD ID", "Receptor Location"]].merge(
        pop_up_info_lim, how='left', left_on='Receptor', right_on='Approved symbol'
    ).drop_duplicates(subset='Human LR Pair', keep="first").drop(columns=["Receptor", "Approved symbol"])
    
    receptor_card_1 = receptor_card[["Human LR Pair", "Receptor name", "Other Symbols"]] 
    receptor_card_2 = receptor_card[["Human LR Pair", "Receptor HGNC ID", "Receptor Location"]] 
    receptor_card_2["HGNC gene card"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url)
    receptor_card_2["Disease relevance"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_disease)
    receptor_card_2["Expression Profile"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_exp)
    # Add Receptor group
    receptor_card_2["Lineage group"] = receptor_card_2['Receptor HGNC ID'].map(receptor_mapping).fillna("none")
    icon_html = '<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left:4px;"></i></a>'
    columns_to_update = ["Receptor HGNC ID"]
    
    for col in columns_to_update:
        receptor_card_2[col] = receptor_card_2[col].str.replace(
            "</a>", icon_html, regex=False
        )
    # Add extra line to balance out cards if necessary
    #receptor_card_2["Lineage group"] = receptor_card_2["Lineage group"].apply(add_spacing)
    receptor_card_2 = receptor_card_2[["Human LR Pair", "Receptor HGNC ID",  "HGNC gene card", "Receptor Location", "Lineage group", "Disease relevance", "Expression Profile" ]]       

    return interaction_card, ligand_card_1, ligand_card_2, receptor_card_1, receptor_card_2

import os

def generate_combined_html_files(
    template_env,
    interaction_card,
    ligand_card_1,
    receptor_card_1,
    ligand_card_2,
    receptor_card_2,
    pubmed_data,
    gene_pair_df,
    output_dir
):
    os.makedirs(output_dir, exist_ok=True)
    column_values = interaction_card["Human LR Pair"].dropna().unique()

    def convert_pair_url(df_pairs):
        df_pairs["Human LR Pair"] = [
            f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/{lrpair}.html" target="_blank" '
            f'title="Open {lrpair} card" style="color: #0000EE; text-decoration: underline;">'
            f'{lrpair}</a>'
            if pd.notna(lrpair) and lrpair.strip() else ""
            for lrpair in df_pairs["Human LR Pair"]
        ]
        return df_pairs

    for value in column_values:
        value1, value2 = value.split()

        # Data rows for each section
        row0 = interaction_card[interaction_card['Human LR Pair'] == value]
        row1 = ligand_card_1[ligand_card_1['Human LR Pair'] == value]
        row2 = receptor_card_1[receptor_card_1['Human LR Pair'] == value]
        row3 = ligand_card_2[ligand_card_2['Human LR Pair'] == value]
        row4 = receptor_card_2[receptor_card_2['Human LR Pair'] == value]

        # Related ligand pairs
        ligand_pairs = gene_pair_df[gene_pair_df['Ligand'] == value1]
        ligand_pairs = ligand_pairs[ligand_pairs["Human LR Pair"] != value]
        ligand_pairs = convert_pair_url(ligand_pairs[["Human LR Pair"]])
        ligand_pairs = ' ・ '.join([btn for btn in ligand_pairs["Human LR Pair"] if btn])

        # Related receptor pairs
        receptor_pairs = gene_pair_df[gene_pair_df['Receptor'] == value2]
        receptor_pairs = receptor_pairs[receptor_pairs["Human LR Pair"] != value]
        receptor_pairs = convert_pair_url(receptor_pairs[["Human LR Pair"]])
        receptor_pairs = ' ・ '.join([btn for btn in receptor_pairs["Human LR Pair"] if btn])

        # Heatmaps
        ligand_image_path = f'data/tabula_sapiens/heatmap/{value1}.html'
        receptor_image_path = f'data/tabula_sapiens/heatmap/{value2}.html'
        ligand_image = open(ligand_image_path).read() if os.path.exists(ligand_image_path) else "Plot does not exist"
        receptor_image = open(receptor_image_path).read() if os.path.exists(receptor_image_path) else "Plot does not exist"

        # Tables
        table0_data = row0.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row0.empty else {}
        table1_data = row1.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row1.empty else {}
        table2_data = row2.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row2.empty else {}
        table3_data = row3.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row3.empty else {}
        table4_data = row4.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row4.empty else {}

        # PubMed Tabs
        row_main = gene_pair_df[gene_pair_df["Human LR Pair"] == value]
        if not row_main.empty:
            row_main = row_main.iloc[0]
            pmids = row_main.get("PMID", "")
            keywords = row_main.get("Relevance Keywords", "")
            cards = row_main.get("Cards", "")

            sources = [pmid.strip() for pmid in pmids.split(',') if pmid.strip()]
            tab_headers = []
            tab_contents = []

            for i, pmid in enumerate(sources):
                pubmed_row = pubmed_data[pubmed_data["PMID"] == pmid]
                if not pubmed_row.empty:
                    title = pubmed_row["Title"].values[0]
                    abstract = pubmed_row["Abstract"].values[0]
                    journal = pubmed_row["Journal"].values[0]
                    year = pubmed_row["Year"].values[0]
                else:
                    title = "No Title Found"
                    abstract = "No Abstract Found"
                    journal = "Journal Unknown"
                    year = "Year Unknown"

                active_class = "active" if i == 0 else ""
                tab_headers.append(f'<button class="tablinks {active_class}" onclick="openTab(event, \'tab{pmid}\')">{pmid}</button>')
                tab_contents.append(f"""
                    <div id="tab{pmid}" class="tabcontent {active_class}">
                        <h2>{title}</h2>
                        <p><strong>{journal}, {year}; <a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">For more details, see PubMed</a></strong></p>
                        <p>{abstract}</p>
                    </div>
                """)
        else:
            tab_headers, tab_contents, keywords, cards = [], [], "", ""

        rendered_content = template_env.render(
            value1=value1,
            value2=value2,
            table0_data=table0_data,
            table1_data=table1_data,
            table2_data=table2_data,
            table3_data=table3_data,
            table4_data=table4_data,
            ligand_image=ligand_image,
            receptor_image=receptor_image,
            ligand_pairs=ligand_pairs,
            receptor_pairs=receptor_pairs,
            tab_headers=''.join(tab_headers),
            tab_contents=''.join(tab_contents),
            cards=cards,
            keywords=keywords
        )

        output_file = os.path.join(output_dir, f"{value1} {value2}.html")
        with open(output_file, 'w') as file:
            file.write(rendered_content)
