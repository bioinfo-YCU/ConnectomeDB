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

# Test or all
test = True
test_genes = ["VEGFA ITGB1", "VEGFA KDR", "VEGFA NRP1", "THPO MPL", "FGF1 FGFR3"] # Example genes
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
    ligand_card_2["Gene Group (HGNC)"] = ligand_card_2['Ligand HGNC ID'].map(ligand_mapping).fillna("none")
    icon_html_card = '<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left:4px;"></i></a>' # Use a different name to avoid conflict
    for col in ["Ligand HGNC ID"]:
        ligand_card_2[col] = ligand_card_2[col].str.replace(
            "</a>", icon_html_card, regex=False
        )
    ligand_card_2 = ligand_card_2[["Human LR Pair", "Ligand HGNC ID", "HGNC gene card", "Ligand Location", "Gene Group (HGNC)", "Disease relevance", "Expression Profile"]]
    receptor_card = gene_pair_input_df[["Human LR Pair", "Receptor", "Receptor name", "Receptor HGNC ID", "Receptor MGI ID", "Receptor RGD ID", "Receptor Location"]].merge(
        pop_up_info_lim, how='left', left_on='Receptor', right_on='Approved symbol'
    ).drop_duplicates(subset='Human LR Pair', keep="first").drop(columns=["Receptor", "Approved symbol"])

    receptor_card_1 = receptor_card[["Human LR Pair", "Receptor name", "Other Symbols"]]
    receptor_card_2 = receptor_card[["Human LR Pair", "Receptor HGNC ID", "Receptor Location"]]
    receptor_card_2["HGNC gene card"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url)
    receptor_card_2["Disease relevance"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_disease)
    receptor_card_2["Expression Profile"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_exp)
    receptor_card_2["Gene Group (HGNC)"] = receptor_card_2['Receptor HGNC ID'].map(receptor_mapping).fillna("none")
    for col in ["Receptor HGNC ID"]:
        receptor_card_2[col] = receptor_card_2[col].str.replace(
            "</a>", icon_html_card, regex=False
        )
    receptor_card_2 = receptor_card_2[["Human LR Pair", "Receptor HGNC ID",  "HGNC gene card", "Receptor Location", "Gene Group (HGNC)", "Disease relevance", "Expression Profile" ]]

    return interaction_card, ligand_card_1, ligand_card_2, receptor_card_1, receptor_card_2


def convert_pair_url(df_pairs):
    """Converts 'Human LR Pair' column to HTML links with interaction IDs."""
    df_pairs = df_pairs.copy()

    # Extract clean Interaction ID from HTML
    df_pairs["Clean Interaction ID"] = df_pairs["Interaction ID"].apply(
        lambda x: re.search(r'(CDB\d+)</a>', x).group(1)
        if isinstance(x, str) and re.search(r'(CDB\d+)</a>', x)
        else None
    )

    # Generate new HTML anchor tag
    def make_link(row):
        lrpair = row["Human LR Pair"]
        interaction_id = row["Clean Interaction ID"]
        if pd.notna(lrpair) and lrpair.strip() and pd.notna(interaction_id):
            encoded_lrpair = lrpair.replace(" ", "%20")
            return (
                f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/'
                f'{encoded_lrpair}_{interaction_id}.html" target="_blank" '
                f'title="Open {lrpair} card" style="color: #0000EE; text-decoration: underline;">'
                f'{lrpair}</a>'
            )
        return ""

    df_pairs["Human LR Pair"] = df_pairs.apply(make_link, axis=1)

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
    # --- Create a map of all pages for navigation ---
   # 1. Clean the 'Interaction ID' column in gene_pair_main_df and add it as a new column
    #    This ensures we have a clean ID for sorting and filename generation.
    cleaned_main_df = gene_pair_main_df.copy()
    cleaned_main_df["Clean Interaction ID"] = cleaned_main_df["Interaction ID"].apply(
        lambda x: re.search(r'(CDB\d+)</a>', x).group(1)
        if isinstance(x, str) and re.search(r'(CDB\d+)</a>', x)
        else None
    )
    # Filter out rows where the clean ID could not be extracted (malformed entries)
    cleaned_main_df = cleaned_main_df.dropna(subset=["Clean Interaction ID"])

    # 2. Sort the DataFrame by the clean Interaction ID to ensure correct navigation order.
    #    This sorted DataFrame will be the single source of truth for iteration.
    df_for_rendering_and_navigation = cleaned_main_df.sort_values(
        by="Clean Interaction ID",
        key=lambda series: series.str[3:].astype(int) # Now it's just "CDBXXXXX" so simple slicing works
    ).reset_index(drop=True)


    rendered_pages = []

    # --- First pass: Collect all content and metadata ---
    for idx, row in gene_pair_keywords_df.iterrows():
        lr_pair_name_hyphen = row["Human LR Pair"]  # e.g., "VEGFA——KDR"
        keywords = row["Relevance Keywords"]
        pmids_str = row["PMID"]
    
        # Convert to space-separated format
        lr_pair_name_space = lr_pair_name_hyphen.replace("——", " ")
        value1, value2 = lr_pair_name_space.split()
    
        # --- Filter card dataframes ---
        row0 = interaction_card_df[interaction_card_df['Human LR Pair'] == lr_pair_name_space]
        if row0.empty:
            print(f"[SKIP] No card data found for: {lr_pair_name_space}")
            continue
    
        # Extract and clean interaction ID from the actual table data
        table0_data = row0.drop('Human LR Pair', axis=1).to_dict(orient='records')[0]
        raw_interaction_id_html = table0_data["Interaction ID"]
        match = re.search(r'(CDB\d+)</a>', raw_interaction_id_html)
        if not match:
            print(f"[SKIP] Could not extract clean interaction ID for: {lr_pair_name_space}")
            continue
        clean_interaction_id = match.group(1)
    
        # Construct filename and output path
        filename = f"{value1.strip()} {value2.strip()}_{clean_interaction_id}.html"
        output_file = os.path.join(output_dir, filename)
    
        # Save for navigation
        rendered_pages.append({
            "interaction_id": clean_interaction_id,
            "lr_pair_name_space": lr_pair_name_space,
            "value1": value1,
            "value2": value2,
            "keywords": keywords,
            "pmids_str": pmids_str,
            "table0_data": table0_data,
            "row1": ligand_card_1_df[ligand_card_1_df['Human LR Pair'] == lr_pair_name_space],
            "row2": receptor_card_1_df[receptor_card_1_df['Human LR Pair'] == lr_pair_name_space],
            "row3": ligand_card_2_df[ligand_card_2_df['Human LR Pair'] == lr_pair_name_space],
            "row4": receptor_card_2_df[receptor_card_2_df['Human LR Pair'] == lr_pair_name_space],
            "output_file": output_file,
            "filename": filename,
            "value1": value1,
            "value2": value2,
        })

    # --- Second pass: Render each page with correct navigation ---
    for i, page in enumerate(rendered_pages):
        prev_page_info = rendered_pages[i - 1] if i > 0 else None
        next_page_info = rendered_pages[i + 1] if i < len(rendered_pages) - 1 else None
    
        # --- Prepare PMID section ---
        tab_headers = []
        tab_contents = []
        pmid_list = [pmid.strip() for pmid in str(page["pmids_str"]).split(',') if pmid.strip()]
        for j, pmid in enumerate(pmid_list):
            pubmed_row = pubmed_data_df[pubmed_data_df["PMID"] == pmid]
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
    
            active_class = "active" if j == 0 else ""
            tab_headers.append(f'<button class="tablinks {active_class}" onclick="openTab(event, \'tab{pmid}\')">{pmid}</button>')
            tab_contents.append(f"""
            <div id="tab{pmid}" class="tabcontent {active_class}">
                <h2>{title}</h2>
                <p><strong>{journal}, {year}; <a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">For more details, see PubMed</a></strong></p>
                <div class="abstract-wrapper">
                    <p class="abstract-content" id="abstract-content-{pmid}">{abstract}</p>
                </div>
            </div>
            """)
    
        # --- Prepare other tables ---
        def get_table_data(df):
            return df.drop('Human LR Pair', axis=1).to_dict(orient='records')[0] if not df.empty else {}
    
        table1_data = get_table_data(page["row1"])
        table2_data = get_table_data(page["row2"])
        table3_data = get_table_data(page["row3"])
        table4_data = get_table_data(page["row4"])
    
        # --- Related pairs ---
        ligand_pairs_df = gene_pair_main_df[(gene_pair_main_df['Ligand'] == page["value1"]) &
                                            (gene_pair_main_df["Human LR Pair"] != page["lr_pair_name_space"])]
        receptor_pairs_df = gene_pair_main_df[(gene_pair_main_df['Receptor'] == page["value2"]) &
                                              (gene_pair_main_df["Human LR Pair"] != page["lr_pair_name_space"])]
    
        ligand_pairs_df = convert_pair_url(ligand_pairs_df)
        receptor_pairs_df = convert_pair_url(receptor_pairs_df)
    
        ligand_pairs_str = ' ・ '.join([btn for btn in ligand_pairs_df["Human LR Pair"] if btn])
        receptor_pairs_str = ' ・ '.join([btn for btn in receptor_pairs_df["Human LR Pair"] if btn])
    
        # --- Expression plots ---
        def get_plot_content(path):
            if os.path.exists(path):
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        return f.read()
                except Exception as e:
                    return f"Plot could not be loaded: {e}"
            return "Plot does not exist"
    
        ligand_image = get_plot_content(f'data/tabula_sapiens/heatmap/{page["value1"]}.html')
        receptor_image = get_plot_content(f'data/tabula_sapiens/heatmap/{page["value2"]}.html')
    
        # --- Final render ---
        rendered_html = template.render(
            GENE_NAME=page["lr_pair_name_space"],
            KEYWORDS=page["keywords"],
            TAB_HEADERS="".join(tab_headers),
            TAB_CONTENTS="".join(tab_contents),
            value1=page["value1"],
            value2=page["value2"],
            table0_data=page["table0_data"],
            table1_data=table1_data,
            table2_data=table2_data,
            table3_data=table3_data,
            table4_data=table4_data,
            ligand_image=ligand_image,
            receptor_image=receptor_image,
            ligand_pairs=ligand_pairs_str,
            receptor_pairs=receptor_pairs_str,
            prev_page_info={
                "interaction_id": prev_page_info["interaction_id"],
                "url": f"https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/{prev_page_info['filename']}",
                "lr_pair_name_space": prev_page_info["lr_pair_name_space"]
            } if prev_page_info else None,
            next_page_info={
                "interaction_id": next_page_info["interaction_id"],
                "url": f"https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/{next_page_info['filename']}",
                "lr_pair_name_space": next_page_info["lr_pair_name_space"]
            } if next_page_info else None
        )
    
        # --- Save HTML file ---
        with open(page["output_file"], "w", encoding="utf-8") as f:
            f.write(rendered_html)


# --- Main Execution Block ---
if __name__ == "__main__":
    if test:
        # Define test genes - these should be in the 'space' format for gene_pair0
        # and will be converted to '——' for gene_pair000 internally.
        # Convert test_genes to '——' format for filtering gene_pair000
        test_genes_hyphen = [gene.replace(" ", "——") for gene in test_genes]
        # Filter gene_pair0 for the test genes to be used in prepare_card_dataframes
        # This gene_pair_input should have space-separated LR pairs
        gene_pair_input = gene_pair0_copy[gene_pair0_copy["Human LR Pair"].isin(test_genes)]
        # Filter gene_pair000 for the test genes (which are in '——' format)
        gene_pair_keywords_filtered = gene_pair000[gene_pair000["Human LR Pair"].isin(test_genes_hyphen)]
    else:
        gene_pair_input = gene_pair0_copy
        # Filter gene_pair000 for the test genes (which are in '——' format)
        gene_pair_keywords_filtered = gene_pair000
        print("Making all cards")

    # Prepare card-specific dataframes
    interaction_card, ligand_card_1, ligand_card_2, receptor_card_1, receptor_card_2 = \
        prepare_card_dataframes(gene_pair_input)

    # Load the merged HTML template
    merged_template = load_template(MERGED_TEMPLATE_PATH)


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
