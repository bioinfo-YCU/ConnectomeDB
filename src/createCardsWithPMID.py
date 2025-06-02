import os
import jinja2
import sys
import pandas as pd
import numpy as np
import time
import base64
import re

# Add the src directory to the path for importing modules
sys.path.append(os.path.abspath("src"))

# Import necessary modules from your existing src files
# Ensure createDataTable and createFunctionalAnnotTable are in your 'src' directory
from createDataTable import pop_up_info, gene_pair0, generate_perplexity_links, gene_pair00
from createFunctionalAnnotTable import gene_pair_annot_ligand, gene_pair_annot_receptor

# --- Paths ---
MERGED_TEMPLATE_PATH = 'HTML/mergedCardWithPMIDTemplate.html'
OUTPUT_DIR = 'data/cards/' # New output directory for combined files

# --- Load and Preprocess Data (Combined from both scripts) ---

# Load PubMed data (from createPMIDpages.py)
pubmed_data = pd.read_csv("data/pubmed_results.csv")
pubmed_data["Year"] = pubmed_data["Year"].astype(str).str.replace(".0", "", regex=False).astype(int)
pubmed_data["PMID"] = pubmed_data["PMID"].astype(str)
pubmed_data = pubmed_data.reset_index(drop=True)

# Load LLM results (from createPMIDpages.py)
bio_keywords = pd.read_csv("data/llm_results.csv")

# --- Prepare gene_pair00 for PMID section (from createPMIDpages.py) ---
# gene_pair00 is used for PMID and Keywords, so it needs the '——' placeholder
# Ensure gene_pair00 is a copy to avoid SettingWithCopyWarning later
gene_pair00_copy = gene_pair00.copy()
gene_pair00_copy["Human LR Pair"] = gene_pair00_copy["Human LR Pair"].str.replace(" ", "——")

# Merge with LLM results
gene_pair000 = gene_pair00_copy.merge(bio_keywords, how='left', left_on="Human LR Pair", right_on='Human LR Pair')
gene_pair000["Relevance Keywords"] = gene_pair000["Relevance Keywords"].astype(str)
gene_pair000["Human LR Pair"] = gene_pair000["Human LR Pair"].astype(str) # Ensure string type

# --- Prepare gene_pair0 for Card section (from createCards.py) ---
# gene_pair0 is used for card details, it should retain spaces for splitting gene names
# Ensure gene_pair0 is a copy to avoid SettingWithCopyWarning later
gene_pair0_copy = gene_pair0.copy()

# Add Disease (specific) to cards
df_disease = pd.read_csv("data/disease_annotations_per_pair.csv")
df_disease = df_disease.groupby('interaction')['disease'].apply(', '.join).reset_index()
mapping_disease = dict(zip(df_disease['interaction'], df_disease['disease']))
gene_pair0_copy["Disease"] = gene_pair0_copy['Human LR Pair'].map(mapping_disease).fillna("unknown")

gene_pair0_copy = generate_perplexity_links(
    gene_pair0_copy,
    pathway_col="Disease",
    default_query_template="What-diseases-is-the-ligand-receptor-pair-{pair}-associated-with"
)

gene_pair0_copy["Interaction ID"] = gene_pair0_copy["Interaction ID"].apply(
    lambda x: f"<a href='https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/database/filter/{x}.html'>{x}</a>"
)
# Add external link icon
icon_html = '<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left:4px;"></i></a>'
columns_to_update = [
    "KEGG Pathway", "PROGENy Pathway", "Cancer-related",
    "Disease Type", "Disease"
]
for col in columns_to_update:
    gene_pair0_copy[col] = gene_pair0_copy[col].str.replace(
        "</a>", icon_html, regex=False
    )

# Add Ligand/Receptor group info
agg_func = lambda x: ', '.join(sorted(set(map(str, x))))
gene_pair_annot_ligand = gene_pair_annot_ligand.groupby('Ligand HGNC ID').agg(agg_func).reset_index()
ligand_mapping = dict(zip(gene_pair_annot_ligand['Ligand HGNC ID'], gene_pair_annot_ligand['Ligand group']))

gene_pair_annot_receptor = gene_pair_annot_receptor.groupby('Receptor HGNC ID').agg(agg_func).reset_index()
receptor_mapping = dict(zip(gene_pair_annot_receptor['Receptor HGNC ID'], gene_pair_annot_receptor['Receptor group']))


# --- Helper Functions (Combined and adjusted) ---

def load_template(template_path):
    """Load Jinja2 template from a file."""
    with open(template_path, 'r', encoding='utf-8') as file:
        return jinja2.Template(file.read())

def encode_image(image_path):
    """Encode an image to base64. (Not used in this version, but kept for reference)"""
    try:
        with open(image_path, "rb") as image_file:
            return base64.b64encode(image_file.read()).decode('utf-8')
    except FileNotFoundError:
        return None

def extract_hgnc_id(col):
    """Use regular expression to extract the HGNC ID after 'HGNC:'."""
    match = re.search(r'HGNC:(\d+)', col)
    if match:
        return match.group(1)
    return None

def convert_hgnc_url(col):
    hgnc_id = extract_hgnc_id(col)
    if hgnc_id:
        visible_text = 'GeneCards <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>'
        new_link = (
            f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}" '
            f'target="_blank" style="color: #0000EE; text-decoration: underline;">{visible_text}</a>'
        )
        return new_link
    return None

def convert_hgnc_url_disease(col):
    hgnc_id = extract_hgnc_id(col)
    if hgnc_id:
        visible_text = 'MalaCards <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>'
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#diseases" target="_blank">{visible_text}</a>'
        return new_link
    return None

def convert_hgnc_url_exp(col):
    hgnc_id = extract_hgnc_id(col)
    if hgnc_id:
        visible_text = 'mRNA expression in normal human tissues <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>'
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#expression" target="_blank">{visible_text}</a>'
        return new_link
    return None

def prepare_card_dataframes(gene_pair_input_df):
    """Prepare interaction, ligand, and receptor dataframes for the card section."""
    # Ensure gene_pair_input_df is a copy to avoid SettingWithCopyWarning
    gene_pair_input_df = gene_pair_input_df.copy()

    gene_pair_input_df["Interaction Type"] = [
        f'{ligand} {ligandLocation} ligand binds to {receptor} {receptorLocation} receptor'
        for ligand, ligandLocation, receptor, receptorLocation in zip(
            gene_pair_input_df["Ligand"], gene_pair_input_df["Ligand Location"],
            gene_pair_input_df["Receptor"], gene_pair_input_df["Receptor Location"]
        )
    ]
    interaction_card = gene_pair_input_df[["Interaction ID", "Human LR Pair", "Interaction Type", "Perplexity", "PMID", "KEGG Pathway",  "PROGENy Pathway", "Cancer-related", "Disease Type", "Disease"]]
    interaction_card["Perplexity"] = interaction_card["Perplexity"].str.replace('size=30', 'size=80')

    pop_up_info_lim = pop_up_info[
        ["Approved symbol", "Alias symbol", "Previous symbol", "Date symbol changed"]
    ].drop_duplicates(subset="Approved symbol", keep="first")

    def format_symbol_aliases(old_symbol, aliases):
        parts = [p for p in (old_symbol, aliases) if p != "N/A"]
        return f"{', '.join(parts)}" if parts else aliases

    pop_up_info_lim['Other Symbols'] = pop_up_info_lim.apply(
        lambda row: format_symbol_aliases(row["Previous symbol"], row["Alias symbol"]),
        axis=1
    )

    ligand_card = gene_pair_input_df[["Human LR Pair", "Ligand", "Ligand name", "Ligand HGNC ID", "Ligand MGI ID", "Ligand RGD ID", "Ligand Location"]].merge(
        pop_up_info_lim, how='left', left_on='Ligand', right_on='Approved symbol'
    ).drop_duplicates(subset='Human LR Pair', keep="first").drop(columns=["Ligand", "Approved symbol"])

    ligand_card_1 = ligand_card[["Human LR Pair", "Ligand name", "Other Symbols" ]]
    ligand_card_2 = ligand_card[["Human LR Pair", "Ligand HGNC ID", "Ligand Location"]]
    ligand_card_2["HGNC gene card"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url)
    ligand_card_2["Disease relevance"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url_disease)
    ligand_card_2["Expression Profile"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url_exp)
    ligand_card_2["Lineage group"] = ligand_card_2['Ligand HGNC ID'].map(ligand_mapping).fillna("none")
    icon_html_card = '<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left:4px;"></i></a>' # Use a different name to avoid conflict
    for col in ["Ligand HGNC ID"]:
        ligand_card_2[col] = ligand_card_2[col].str.replace(
            "</a>", icon_html_card, regex=False
        )
    ligand_card_2 = ligand_card_2[["Human LR Pair", "Ligand HGNC ID", "HGNC gene card", "Ligand Location", "Lineage group", "Disease relevance", "Expression Profile"]]


    receptor_card = gene_pair_input_df[["Human LR Pair", "Receptor", "Receptor name", "Receptor HGNC ID", "Receptor MGI ID", "Receptor RGD ID", "Receptor Location"]].merge(
        pop_up_info_lim, how='left', left_on='Receptor', right_on='Approved symbol'
    ).drop_duplicates(subset='Human LR Pair', keep="first").drop(columns=["Receptor", "Approved symbol"])

    receptor_card_1 = receptor_card[["Human LR Pair", "Receptor name", "Other Symbols"]]
    receptor_card_2 = receptor_card[["Human LR Pair", "Receptor HGNC ID", "Receptor Location"]]
    receptor_card_2["HGNC gene card"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url)
    receptor_card_2["Disease relevance"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_disease)
    receptor_card_2["Expression Profile"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_exp)
    receptor_card_2["Lineage group"] = receptor_card_2['Receptor HGNC ID'].map(receptor_mapping).fillna("none")
    for col in ["Receptor HGNC ID"]:
        receptor_card_2[col] = receptor_card_2[col].str.replace(
            "</a>", icon_html_card, regex=False
        )
    receptor_card_2 = receptor_card_2[["Human LR Pair", "Receptor HGNC ID",  "HGNC gene card", "Receptor Location", "Lineage group", "Disease relevance", "Expression Profile" ]]

    return interaction_card, ligand_card_1, ligand_card_2, receptor_card_1, receptor_card_2


def convert_pair_url(df_pairs):
    """Converts 'Human LR Pair' column to HTML links."""
    # Ensure df_pairs is a copy to avoid SettingWithCopyWarning
    df_pairs = df_pairs.copy()
    df_pairs["Human LR Pair"] = [
        f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/{lrpair.replace(" ", "%20")}.html" target="_blank" '
        f'title="Open {lrpair} card" style="color: #0000EE; text-decoration: underline;">'
        f'{lrpair}</a>'
        if pd.notna(lrpair) and lrpair.strip() else ""
        for lrpair in df_pairs["Human LR Pair"]
    ]
    return df_pairs


def generate_combined_html_files(
    gene_pair_keywords_df, # Corresponds to gene_pair000 (with '——' in LR Pair)
    template,
    interaction_card_df,
    ligand_card_1_df,
    receptor_card_1_df,
    ligand_card_2_df,
    receptor_card_2_df,
    pubmed_data_df,
    gene_pair_main_df, # Corresponds to gene_pair0 (with spaces in LR Pair)
    output_dir
):
    """
    Generate combined HTML pages with PMID details on top and card details at the bottom
    for each Human LR Pair.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over the DataFrame that contains PMID, Keywords, and the LR Pair (with '——')
    for idx, row in gene_pair_keywords_df.iterrows():
        lr_pair_name_hyphen = row["Human LR Pair"] # e.g., "VEGFA——KDR"
        keywords = row["Relevance Keywords"]
        pmids_str = row["PMID"]

        # Convert LR pair name to space-separated for display and file naming
        lr_pair_name_space = lr_pair_name_hyphen.replace("——", " ") # e.g., "VEGFA KDR"
        value1, value2 = lr_pair_name_space.split() # Ligand and Receptor names

        # --- Prepare PMID Section Data ---
        tab_headers = []
        tab_contents = []
        sources = [pmid.strip() for pmid in str(pmids_str).split(',') if pmid.strip()]

        if sources:
            for i, pmid in enumerate(sources):
                pubmed_row = pubmed_data_df[pubmed_data_df["PMID"] == pmid]
                
                if not pubmed_row.empty:
                    title = pubmed_row["Title"].values[0]
                    abstract = pubmed_row["Abstract"].values[0]
                    journal = pubmed_row["Journal"].values[0]
                    year = pubmed_row["Year"].values[0]
                    # Add this line to check abstract length
                    print(f"PMID: {pmid}, Abstract Length (characters): {len(abstract)}, Abstract Length (words): {len(abstract.split())}")
                else:
                    title = "No Title Found"
                    abstract = "No Abstract Found"
                    journal = "Journal Unknown"
                    year = "Year Unknown"
                    print(f"PMID: {pmid}, No abstract found.") # Also print for missing abstracts
        # .                

                active_class = "active" if i == 0 else ""
                tab_headers.append(f'<button class="tablinks {active_class}" onclick="openTab(event, \'tab{pmid}\')">{pmid}</button>')
                
                # --- MODIFIED: Abstract content with wrapper and button ---
                tab_contents.append(f"""
    <div id="tab{pmid}" class="tabcontent {active_class}">
        <h2>{title}</h2>
        <p><strong>{journal}, {year}; <a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">For more details, see PubMed</a></strong></p>
        <div class="abstract-wrapper">
            <p class="abstract-content" id="abstract-content-{pmid}">{abstract}</p>
            <button class="read-more-btn" id="read-more-btn-{pmid}" onclick="toggleAbstract('abstract-content-{pmid}', 'read-more-btn-{pmid}')">Read more</button>
        </div>
    </div>
""")
                # --- END MODIFIED ---

        # --- Prepare Card Section Data ---
        # Filter card dataframes using the space-separated LR pair name
        row0 = interaction_card_df[interaction_card_df['Human LR Pair'] == lr_pair_name_space]
        row1 = ligand_card_1_df[ligand_card_1_df['Human LR Pair'] == lr_pair_name_space]
        row2 = receptor_card_1_df[receptor_card_1_df['Human LR Pair'] == lr_pair_name_space]
        row3 = ligand_card_2_df[ligand_card_2_df['Human LR Pair'] == lr_pair_name_space]
        row4 = receptor_card_2_df[receptor_card_2_df['Human LR Pair'] == lr_pair_name_space]

        table0_data = row0.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row0.empty else {}
        table1_data = row1.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row1.empty else {}
        table2_data = row2.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row2.empty else {}
        table3_data = row3.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row3.empty else {}
        table4_data = row4.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not row4.empty else {}

        # Related ligand pairs (using gene_pair_main_df which has space-separated LR pairs)
        ligand_pairs_df = gene_pair_main_df[gene_pair_main_df['Ligand'] == value1]
        # Filter out the current LR pair (which is space-separated in gene_pair_main_df)
        ligand_pairs_df = ligand_pairs_df[ligand_pairs_df["Human LR Pair"] != lr_pair_name_space]
        ligand_pairs_df = convert_pair_url(ligand_pairs_df[["Human LR Pair"]])
        ligand_pairs_str = ' ・ '.join([btn for btn in ligand_pairs_df["Human LR Pair"] if btn])

        # Related receptor pairs (using gene_pair_main_df which has space-separated LR pairs)
        receptor_pairs_df = gene_pair_main_df[gene_pair_main_df['Receptor'] == value2]
        # Filter out the current LR pair (which is space-separated in gene_pair_main_df)
        receptor_pairs_df = receptor_pairs_df[receptor_pairs_df["Human LR Pair"] != lr_pair_name_space]
        receptor_pairs_df = convert_pair_url(receptor_pairs_df[["Human LR Pair"]])
        receptor_pairs_str = ' ・ '.join([btn for btn in receptor_pairs_df["Human LR Pair"] if btn])

        # Expression Plots (assuming HTML content is read from files)
        # Uncomment the following block if you have these plot files and want to include them
        ligand_image_path = f'data/tabula_sapiens/heatmap/{value1}.html'
        receptor_image_path = f'data/tabula_sapiens/heatmap/{value2}.html'

        ligand_image = "Plot does not exist"
        if os.path.exists(ligand_image_path):
            try:
                with open(ligand_image_path, "r", encoding='utf-8') as html_file:
                    ligand_image = html_file.read()
            except Exception as e:
                print(f"Error reading ligand image HTML for {value1}: {e}")
                ligand_image = "Plot does not exist (Error reading file)"


        receptor_image = "Plot does not exist"
        if os.path.exists(receptor_image_path):
            try:
                with open(receptor_image_path, "r", encoding='utf-8') as html_file:
                    receptor_image = html_file.read()
            except Exception as e:
                print(f"Error reading receptor image HTML for {value2}: {e}")
                receptor_image = "Plot does not exist (Error reading file)"


        # Render the template with all collected data
        rendered_content = template.render(
            GENE_NAME=lr_pair_name_space, # For PMID section title
            KEYWORDS=keywords, # For PMID section keywords
            TAB_HEADERS="".join(tab_headers),
            TAB_CONTENTS="".join(tab_contents),
            value1=value1, # For card section ligand name
            value2=value2, # For card section receptor name
            table0_data=table0_data,
            table1_data=table1_data,
            table2_data=table2_data,
            table3_data=table3_data,
            table4_data=table4_data,
            ligand_image=ligand_image,
            receptor_image=receptor_image,
            ligand_pairs=ligand_pairs_str,
            receptor_pairs=receptor_pairs_str
        )

        # Save the HTML file
        output_file = os.path.join(output_dir, f"{lr_pair_name_space}.html")
        with open(output_file, 'w', encoding='utf-8') as file:
            file.write(rendered_content)

# --- Main Execution Block ---
if __name__ == "__main__":
    # Define test genes - these should be in the 'space' format for gene_pair0
    # and will be converted to '——' for gene_pair000 internally.
    test_genes = ["VEGFA KDR", "ADAM17 IL6R"] # Example genes

    # Filter gene_pair0 for the test genes to be used in prepare_card_dataframes
    # This gene_pair_input should have space-separated LR pairs
    gene_pair_input = gene_pair0_copy[gene_pair0_copy["Human LR Pair"].isin(test_genes)]

    # Prepare card-specific dataframes
    interaction_card, ligand_card_1, ligand_card_2, receptor_card_1, receptor_card_2 = \
        prepare_card_dataframes(gene_pair_input)

    # Load the merged HTML template
    merged_template = load_template(MERGED_TEMPLATE_PATH)

    # Filter gene_pair000 for the test genes (which are in '——' format)
    # Convert test_genes to '——' format for filtering gene_pair000
    test_genes_hyphen = [gene.replace(" ", "——") for gene in test_genes]
    gene_pair_keywords_filtered = gene_pair000[gene_pair000["Human LR Pair"].isin(test_genes_hyphen)]

    # Generate combined HTML files
    generate_combined_html_files(
        gene_pair_keywords_df=gene_pair_keywords_filtered,
        template=merged_template,
        interaction_card_df=interaction_card,
        ligand_card_1_df=ligand_card_1,
        receptor_card_1_df=receptor_card_1,
        ligand_card_2_df=ligand_card_2,
        receptor_card_2_df=receptor_card_2,
        pubmed_data_df=pubmed_data,
        gene_pair_main_df=gene_pair0_copy, # Use the original gene_pair0_copy for related pairs
        output_dir=OUTPUT_DIR
    )
    print(f"Generated combined HTML files in: {OUTPUT_DIR}")
