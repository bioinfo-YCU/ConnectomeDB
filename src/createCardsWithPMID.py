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
from fetchGSheet import gene_pair_mouse, conservation, ligand_loc, receptor_loc
from createDataTable import pop_up_info, gene_pair0, site_url, generate_perplexity_links, gene_pair00, is_mouse_specific, grab_mouse_info, dedup_locations
from createFunctionalAnnotTable import gene_pair_annot_ligand, gene_pair_annot_receptor

# Test or all
test = False
test_genes = ["H60a Klrk1", "H60b Klrk1", "VEGFA NRP1", "THPO MPL", "FGF1 FGFR3", "Lair1 Lilrb4a"] # Example genes
# --- Paths ---
MERGED_TEMPLATE_PATH = 'HTML/mergedCard_tabs.html'
OUTPUT_DIR = 'data/cards/' # New output directory for combined files
OUTPUT_DIR_h = 'data/cards/human' 
OUTPUT_DIR_m = 'data/cards/mouse'

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR_h, exist_ok=True)
os.makedirs(OUTPUT_DIR_m, exist_ok=True)

gene_pair_per_row = gene_pair_mouse[["LR_pair_orig", "PMID", "lig_species", "rec_species", "ligand_orig", "receptor_orig", "LR Pair Card"]]

# function for replacing visible text:
def update_link_text_with_symbol(html_str, new_symbol):
    """
    Replace the visible text in an anchor tag with the provided symbol.
    The HGNC ID is extracted from the href and left unchanged.
    """
    # Only proceed if input is valid
    if not isinstance(html_str, str) or not html_str.strip():
        return html_str
    
    return re.sub(r">(.*?)</a>", f">{new_symbol}</a>", html_str) #<i class='fa-solid fa-arrow-up-right-from-square' style='margin-left:4px;'></i></a>

# Create Perplexity link

# Recreate Perplexity link
def create_url_basic(symbol):
          label = "Perplexity (LLM)"
          query = f"What diseases is {symbol} implicated in?"
          encoded_query = query.replace(" ", "%20")
          output = f'<a href="https://www.perplexity.ai/search?q={encoded_query}" target="_blank">{label}</a>' #<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left:4px;"></i>
          return output

# --- Load and Preprocess Data (Combined from both scripts) ---

# Load PubMed data (from createPMIDpages.py)
pubmed_data = pd.read_csv("data/pubmed_results.csv")
pubmed_data["Year"] = pubmed_data["Year"].astype(str).str.replace(".0", "", regex=False).astype(int)
pubmed_data["PMID"] = pubmed_data["PMID"].astype(str)
pubmed_data = pubmed_data.reset_index(drop=True)

### For latest DB, skip
# # Load LLM results (from createPMIDpages.py)
# bio_keywords = pd.read_csv("data/llm_results.csv")

# --- Prepare gene_pair00 for PMID section (from createPMIDpages.py) ---
# gene_pair00 is used for PMID and Keywords, so it needs the '—' placeholder
# Ensure gene_pair00 is a copy to avoid SettingWithCopyWarning later
gene_pair00_copy = gene_pair00.copy()
#gene_pair00_copy["Human LR Pair"] = gene_pair00_copy["Human LR Pair"].str.replace(" ", "—")

# Merge with LLM results
# gene_pair000 = gene_pair00_copy.merge(bio_keywords, how='left', left_on="Human LR Pair", right_on='Human LR Pair')
# gene_pair000["Relevance Keywords"] = gene_pair000["Relevance Keywords"].astype(str)
 # Ensure string type
# gene_pair000["Human LR Pair"] = gene_pair000["Human LR Pair"].astype(str)

### For latest DB,
gene_pair000 = gene_pair00.copy()
gene_pair000["LR Pair Card"] = gene_pair000["LR Pair Card"].astype(str)
# --- Prepare gene_pair0 for Card section (from createCards.py) ---
# gene_pair0 is used for card details, it should retain spaces for splitting gene names
# Ensure gene_pair0 is a copy to avoid SettingWithCopyWarning later
gene_pair0_copy = gene_pair0.copy()

# grab the pairs by interaction id that has human evidence that is absent in human
# mouse_interaction_ids = gene_pair0_copy[gene_pair0_copy['Ligand'].apply(is_mouse_specific)]['Interaction ID'].tolist()
mouse_interaction_ids = gene_pair0_copy["Interaction ID"][gene_pair0_copy["Human evidence"].isin(["absent in human", "not conserved"])]

### For for latest DB, remove for now
# # Add Disease (specific) to cards
# df_disease = pd.read_csv("data/disease_annotations_per_pair.csv")
# df_disease = df_disease.groupby('interaction')['disease'].apply(', '.join).reset_index()
# mapping_disease = dict(zip(df_disease['interaction'], df_disease['disease']))
# gene_pair0_copy["Disease"] = gene_pair0_copy['Human LR Pair'].map(mapping_disease).fillna("unknown")

# gene_pair0_copy = generate_perplexity_links(
#     gene_pair0_copy,
#     pathway_col="Disease",
#     default_query_template="What-diseases-is-the-ligand-receptor-pair-{pair}-associated-with"
# )

### SHOULD BE ACTIVATED ONCE WE DECIDE TO OPEN DB 
# # Hide for now (linking to actual PMID database
# gene_pair0_copy["Interaction ID"] = gene_pair0_copy["Interaction ID"].apply(
#     lambda x: f"<a href='{site_url}database/filter/{x}.html'>{x}</a>"
# )

# Add external link icon
# icon_html = '<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left:4px;"></i></a>'
# columns_to_update = [
#     "KEGG Pathway", "PROGENy Pathway", "Cancer-related",
#     "Disease Type", "Disease"
# ]
# for col in columns_to_update:
#     gene_pair0_copy[col] = gene_pair0_copy[col].str.replace(
#         "</a>", icon_html, regex=False
#     )

# Add missing ligand name for mouse-specific
mouse_info = pd.read_csv("data/MRK_Merged_MGI_DB.tsv", sep="\t", dtype=str)
mouse_info["Aliases"] = mouse_info["Aliases"].str.replace("|", ", ", regex=False)
mapping_mouse_name = dict(zip(mouse_info['MGI Marker Accession ID'], mouse_info['Marker Name']))
mapping_mouse_aliases = dict(zip(mouse_info['MGI Marker Accession ID'], mouse_info['Aliases']))
mapping_mouse_uniprot = dict(zip(mouse_info['MGI Marker Accession ID'], mouse_info['UniProt IDs']))
mapping_mouse_ncbi = dict(zip(mouse_info['MGI Marker Accession ID'], mouse_info['EntrezGene ID']))

# Combine ligand and receptor columns from mouse DB input
gene_pair_mouse_card = gene_pair_mouse[gene_pair_mouse["LR Pair Card"].isin(grab_mouse_info)]
mouse_id_map = pd.concat([
    gene_pair_mouse_card[['MGI ligand', 'ENSEMBL ligand']].rename(columns={
        'MGI ligand': 'MGI ID', 'ENSEMBL ligand': 'ENSEMBL ID'
    }),
    gene_pair_mouse_card[['MGI receptor', 'ENSEMBL receptor']].rename(columns={
        'MGI receptor': 'MGI ID', 'ENSEMBL receptor': 'ENSEMBL ID'
    })
])

# Drop duplicate rows
mapping_mouse_ens = mouse_id_map.drop_duplicates().reset_index(drop=True)

### For the latest DB, remove for now
# mouse_rat_info = pd.read_csv("data/mouse_to_rat_mapping.csv")
# mapping_mr2 = dict(zip(mouse_rat_info['RGD ID'], mouse_rat_info['RGD name']))


# extract mgi
def extract_mgi_id(col):
    """Use regular expression to extract the HGNC ID after 'HGNC:'."""
    match = re.search(r'MGI:(\d+)', col)
    if match:
        return 'MGI:' +str(match.group(1))
    return None

### For latest DB, no need
# gene_pair0_copy['Ligand MGI ID'] = gene_pair0_copy['Ligand MGI ID'].apply(extract_mgi_id)
# gene_pair0_copy['Receptor MGI ID'] = gene_pair0_copy['Receptor MGI ID'].apply(extract_mgi_id)

gene_pair_mouse_card_lim = gene_pair_mouse_card[["LR Pair Card", "Mouse_ligand", "Mouse_receptor", "MGI ligand", "MGI receptor"]]
gene_pair0_copy =gene_pair0_copy.merge(gene_pair_mouse_card_lim, how="left", on="LR Pair Card")
gene_pair0_copy = gene_pair0_copy.rename(columns={
        'MGI ligand': 'Ligand MGI ID',
        'MGI receptor': 'Receptor MGI ID'
    })
### for the latest DB, remove this and use new one instead
# Apply the mapping to 'Ligand Name'
# gene_pair0_copy['Ligand Name'] = gene_pair0_copy.apply(
#     lambda row: mapping_mouse_name.get(row['Ligand MGI ID'], row['Ligand Name'])
#     if pd.isna(row['Ligand Name']) or str(row['Ligand Name']).strip() == '' else row['Ligand Name'],
#     axis=1
# )

#  # Apply the mapping to 'Receptor Name'
# gene_pair0_copy['Receptor Name'] = gene_pair0_copy.apply(
#     lambda row: mapping_mouse_name.get(row['Receptor MGI ID'], row['Receptor Name'])
#     if pd.isna(row['Receptor Name']) or str(row['Receptor Name']).strip() == '' else row['Receptor Name'],
#     axis=1
# )

gene_pair0_copy['Ligand Name'] = gene_pair0_copy.apply(
    lambda row: mapping_mouse_name.get(row['Ligand MGI ID'], row['Ligand Name'])
    if pd.notna(row['Ligand MGI ID']) else row['Ligand Name'],
    axis=1
)
gene_pair0_copy['Receptor Name'] = gene_pair0_copy.apply(
    lambda row: mapping_mouse_name.get(row['Receptor MGI ID'], row['Receptor Name'])
    if pd.notna(row['Receptor MGI ID']) else row['Receptor Name'],
    axis=1
)
# for grouping (flattening file)
agg_func = lambda x: ', '.join(sorted(set(map(str, x))))

## add the conservation info
conservation = conservation[["LR Pair Card", "Direct", "Conserved"]]
conservation = conservation.groupby('LR Pair Card').agg(agg_func).reset_index()
# Sample data (assuming your DataFrame is named df)
def clean_species_column(series):
    def clean_cell(cell):
        if pd.isna(cell) or cell.strip() == "":
            return ""
        # Split by comma, strip whitespace, and remove empty entries
        items = [x.strip() for x in cell.split(",") if x.strip()]
        # Remove duplicates while preserving order
        seen = set()
        deduped = []
        for item in items:
            if item not in seen:
                seen.add(item)
                deduped.append(item)
        return ", ".join(deduped)
    
    return series.apply(clean_cell)

# Apply to Direct and Conserved columns
conservation["Direct"] = clean_species_column(conservation["Direct"])
conservation["Conserved"] = clean_species_column(conservation["Conserved"])
gene_pair0_copy = gene_pair0_copy.merge(conservation, how= 'left', on = "LR Pair Card")


# Replace known null-like strings with NaN
col = gene_pair0_copy["Conserved"].astype(str).str.strip().str.lower()
# Replace known null-like strings with NaN
col = col.replace(["", "na", "nan", "null", "none"], np.nan)
# Fill any remaining NaNs with "none"
gene_pair0_copy["Conserved"] = col.fillna("none")

# Make sure the Ligand and Receptor are the LR Pair
# Split the "LR Pair Card" into two new columns: "Ligand" and "Receptor"
gene_pair0_copy[["Ligand", "Receptor"]] = gene_pair0_copy["LR Pair Card"].str.split(" ", n=1, expand=True)

### Adding ligand loc with source
def append_source_to_location(row, geneloc, loc_col):
    ligand = row[loc_col]
    locations = [loc.strip() for loc in row["location"].split(',')]
    sources = [src.strip() for src in row["source"].split(',')]

    # Deduplicate locations while preserving order
    unique_locations = list(dict.fromkeys(locations))

    # Create a list to hold location and source pairs
    location_with_sources = []
    for loc in unique_locations:
        matching_rows = geneloc[
            (geneloc[loc_col] == ligand) &
            (geneloc["location"].str.contains(loc))
        ]
        all_sources = matching_rows["source"].unique()
        sources_str = ", ".join(sorted(set(all_sources)))
        location_with_sources.append(f"{loc} based on {sources_str}")
    
    return ", ".join(location_with_sources)

#### Ligand
ligand_loc = ligand_loc.dropna(axis=1, how='all')

grouped = ligand_loc.groupby("Ligand").agg({
    "location": lambda x: dedup_locations(', '.join(x)),
    "source": lambda x: ', '.join(x)
}).reset_index()

# Append source information to location
grouped["Ligand Location"] = grouped.apply(
    lambda row: append_source_to_location(row, ligand_loc, loc_col="Ligand"),
    axis=1
)

mapping_loc = dict(zip(grouped['Ligand'], grouped['Ligand Location']))
gene_pair0_copy['Ligand Location'] = gene_pair0_copy['Ligand'].replace(mapping_loc)

#### Receptor
ligand_loc = receptor_loc.dropna(axis=1, how='all')

grouped = receptor_loc.groupby("Receptor").agg({
    "location": lambda x: dedup_locations(', '.join(x)),
    "source": lambda x: ', '.join(x)
}).reset_index()

# Append source information to location
grouped["Receptor Location"] = grouped.apply(
    lambda row: append_source_to_location(row, ligand_loc, loc_col="Receptor"),
    axis=1
)

mapping_loc = dict(zip(grouped['Receptor'], grouped['Receptor Location']))
gene_pair0_copy['Receptor Location'] = gene_pair0_copy['Receptor'].replace(mapping_loc)

#########################################

# Add Ligand/Receptor group info
gene_pair_annot_ligand = gene_pair_annot_ligand.groupby('Ligand HGNC ID').agg(agg_func).reset_index()
ligand_mapping = dict(zip(gene_pair_annot_ligand['Ligand HGNC ID'], gene_pair_annot_ligand['Ligand group']))

gene_pair_annot_receptor = gene_pair_annot_receptor.groupby('Receptor HGNC ID').agg(agg_func).reset_index()
receptor_mapping = dict(zip(gene_pair_annot_receptor['Receptor HGNC ID'], gene_pair_annot_receptor['Receptor group']))

# --- Add ENSP info for Jensen lab DISEASES to pop_info_lim ---
ensp_df= pd.read_csv("data/hgnc_ensp_biomart.csv")
# add linking
ensp_df["ensembl_peptide_id"] = ensp_df["ensembl_peptide_id"].apply(
    lambda x: f"<a href='https://diseases.jensenlab.org/Entity?order=textmining,knowledge,experiments&textmining=10&knowledge=10&experiments=10&type1=9606&type2=-26&id1={x}' target='_blank'>DISEASES (text mining)</a>" #<i class='fa-solid fa-arrow-up-right-from-square' style='margin-left:4px;'></i>
)

ensp_df = ensp_df.groupby('hgnc_id').agg(agg_func).reset_index()
ensp_df= ensp_df[["hgnc_id", "ensembl_peptide_id"]]
pop_up_info= pop_up_info.merge(
        ensp_df, how='left', left_on='HGNC ID', right_on='hgnc_id'
    ).drop_duplicates(subset='hgnc_id', keep="first")
pop_up_info = pop_up_info.rename(columns={"ensembl_peptide_id": "JensenLab DISEASES"})

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

def convert_hgnc_url(row, symbol_col, hgnc_col):
    symbol = row[symbol_col]
    hgnc_id = extract_hgnc_id(row[hgnc_col])
    
    if hgnc_id:
        visible_text = f'{symbol}' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = (
            f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}" '
            f'target="_blank">{visible_text}</a>'
        )
        return new_link
    return None


def convert_mgi_url(row, symbol_col, mgi_col):
    symbol = row[symbol_col]
    mgi_id = row[mgi_col]
    
    if mgi_id:
        visible_text = f'{symbol}' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = (
            f'<a href="https://www.informatics.jax.org/marker/{mgi_id}" '
            f'target="_blank">{visible_text}</a>'
        )
        return new_link
    return None

def convert_mgi_GO_url(row, mgi_col):
    mgi_id = row[mgi_col]
    if mgi_id:
        visible_text = 'Mouse Genome Informatics' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = (
            f'<a href="https://www.informatics.jax.org/go/marker/{mgi_id}" '
            f'target="_blank">{visible_text}</a>'
        )
        return new_link
    return None

def convert_GEO_url(row, ncbi_col):
    ncbi_id = row[ncbi_col]
    if ncbi_id:
        visible_text = 'GEO (RNA and protein)' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = (
            f'<a href="https://www.ncbi.nlm.nih.gov/geoprofiles?LinkName=gene_geoprofiles&from_uid={ncbi_id}" '
            f'target="_blank">{visible_text}</a>'
        )
        return new_link
    return None

def convert_quickGO_url(row, uniprot_col):
    uniprot_id = row[uniprot_col]
    if uniprot_id:
        visible_text = 'quickGO' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = (
            f'<a href="https://www.ebi.ac.uk/QuickGO/annotations?geneProductId={uniprot_id}" '
            f'target="_blank">{visible_text}</a>'
        )
        return new_link
    return None

def convert_ncbi_url(row, symbol_col, ncbi_col):
    symbol = row[symbol_col]
    ncbi_id = row[ncbi_col]
    
    if ncbi_id:
        visible_text = f'{symbol}' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = (
            f'<a href="https://www.ncbi.nlm.nih.gov/gene/{ncbi_id}" '
            f'target="_blank">{visible_text}</a>'
        )
        return new_link
    return None

def convert_ensm_url_exp(row, ensm_col):
    ensm_id = row[ensm_col]
    # Check if ensm_id is NaN or empty
    if pd.isna(ensm_id) or not ensm_id:
        visible_text = f'Expression Atlas (RNA) <i>(no ensemble gene id)</i>'
        return visible_text
    visible_text = 'Expression Atlas (RNA)' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
    new_link = f'<a href="https://www.ebi.ac.uk/gxa/genes/{ensm_id}" target="_blank">{visible_text}</a>'
    return new_link
    return None

def convert_mgi_ensembl_url(row, symbol_col, ensm_col):
    symbol = row[symbol_col]
    ensm_id = row[ensm_col]
    
    # Check if ensm_id is NaN or empty
    if pd.isna(ensm_id) or not ensm_id:
        visible_text = f'{symbol} <i>(not in the official report)</i>'
        return visible_text
    
    # If ensm_id exists, create the link
    visible_text = f'{symbol}' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
    new_link = (
        f'<a href="https://www.ensembl.org/id/{ensm_id}" '
        f'target="_blank">{visible_text}</a>'
    )
    return new_link

def convert_symbol_url_allenBrain(row, symbol_col):
    symbol = row[symbol_col]
    if symbol:
        visible_text = 'Allen Brain Atlas (RNA)'
        new_link = (
            f'<a href="https://celltypes.brain-map.org/rnaseq/searches?{{'
            f'%22exact_match%22:true,%22search_term%22:%22{symbol}%22,'
            f'%22search_type%22:%22gene%22,%22features%22:[],%22tumors%22:[],%22page_num%22:0'
            f'}}" target="_blank">{visible_text}</a>'
        )
        return new_link
    return None

def convert_symbol_url_MCA(row, symbol_col):
    symbol = row[symbol_col]
    if symbol:
        visible_text = 'Mouse Cell Atlas 1.1 (RNA)'
        new_link = (
            f'<a href="https://bis.zju.edu.cn/MCA/assets/img/mca1.1/tsne/gene-pic/{symbol}.jpg" target="_blank">{visible_text}</a>'
        )
        return new_link
    return None

    
# for Malacards
def convert_symbol_url_disease(symbol):
    if symbol:
        visible_text = 'MalaCards' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = f'<a href="https://www.malacards.org/search/results?q={symbol}" target="_blank">{visible_text}</a>'
        return new_link
    return None

# for Human cell atlas
def convert_symbol_url_exp(symbol):
    if symbol:
        visible_text = 'Human Cell Atlas (RNA)' #<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;">
        new_link = f'<a href="https://cellxgene.cziscience.com/gene-expression?ver=2&genes={symbol}" target="_blank">{visible_text}</a>'
        return new_link
    return None

# for GEPIA, cancer exp
def convert_symbol_url_exp_GEPIA(symbol):
    if symbol:
        visible_text = 'GEPIA cancer vs normal (RNA)' #<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"> 
        new_link = f'<a href="http://gepia.cancer-pku.cn/detail.php?gene={symbol}" target="_blank">{visible_text}</a>'
        return new_link
    return None


# for OMIM diseases
def convert_omim_url_disease(omim_id):
    if omim_id:
        visible_text = 'OMIM' #<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = f'<a href="https://omim.org/entry/{omim_id}" target="_blank">{visible_text}</a>'
        return new_link
    return None

pop_up_info['OMIM'] = pop_up_info["omim_id"].apply(convert_omim_url_disease)

# For open targets platform diseases
def convert_ensg_url_disease(ensg_id):
    if ensg_id:
        visible_text = 'Open Targets Platform' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = f'<a href="https://platform.opentargets.org/target/{ensg_id}/associations" target="_blank">{visible_text}</a>'
        return new_link
    return None

pop_up_info['Open Targets Platform'] = pop_up_info["ensembl_gene_id"].apply(convert_ensg_url_disease)

# For the Human protein atlas
def convert_ensg_url_exp(ensg_id):
    if ensg_id:
        visible_text = 'Human Protein Atlas (RNA and protein)' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = f'<a href="https://www.proteinatlas.org/{ensg_id}" target="_blank">{visible_text}</a>'
        return new_link
    return None

pop_up_info['Human Protein Atlas'] = pop_up_info["ensembl_gene_id"].apply(convert_ensg_url_exp)

## Gene (Protein) Ontology
def convert_uniprot_url_GO(uniprot_id):
    if uniprot_id:
        visible_text = 'AmiGO 2' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = f'<a href="https://amigo.geneontology.org/amigo/gene_product/UniProtKB:{uniprot_id}" target="_blank">{visible_text}</a>'
        return new_link
    return None

pop_up_info['AmiGO'] = pop_up_info["uniprot_ids"].apply(convert_uniprot_url_GO)

def convert_uniprot_url_panGO(uniprot_id):
    if uniprot_id:
        visible_text = 'PAN-GO' #  <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = f'<a href="https://functionome.geneontology.org/gene/UniProtKB:{uniprot_id}" target="_blank">{visible_text}</a>'
        return new_link
    return None

pop_up_info['PAN-GO'] = pop_up_info["uniprot_ids"].apply(convert_uniprot_url_panGO)


def convert_hgnc_url_disease(col):
    hgnc_id = extract_hgnc_id(col)
    if hgnc_id:
        visible_text = 'GeneCards Diseases' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#diseases" target="_blank">{visible_text}</a>'
        return new_link
    return None

def convert_hgnc_url_exp(col):
    hgnc_id = extract_hgnc_id(col)
    if hgnc_id:
        visible_text = 'GeneCards (RNA and protein)' # <i class="fa-solid fa-arrow-up-right-from-square" style="margin-left: 4px;"></i>
        new_link = f'<a href="https://www.genecards.org/cgi-bin/carddisp.pl?id_type=hgnc&id={hgnc_id}#expression" target="_blank">{visible_text}</a>'
        return new_link
    return None


def prepare_card_dataframes(gene_pair_input_df, mouse_interaction_ids=None):
    """
    Prepare interaction, ligand, and receptor dataframes for the card section.
    
    Args:
        gene_pair_input_df: Main dataframe with gene pair information
        mouse_interaction_ids: List of interaction IDs that are mouse-specific (optional)
    """
    # Ensure gene_pair_input_df is a copy to avoid SettingWithCopyWarning
    gene_pair_input_df = gene_pair_input_df.copy()
    
    # Determine if we're dealing with mouse-specific interactions
    if mouse_interaction_ids is None:
        mouse_interaction_ids = []
    
    # Check if current rows contain mouse-specific interactions
    gene_pair_input_df['is_mouse_specific'] = gene_pair_input_df['Interaction ID'].isin(mouse_interaction_ids)

    gene_pair_input_df["Interaction Type"] = [
        f'{ligand} {ligandLocation} ligand binds to {receptor} {receptorLocation} receptor; Direct Pair in {direct} and Conserved in {conserved}'
        for ligand, ligandLocation, receptor, receptorLocation, direct, conserved in zip(
            gene_pair_input_df["Ligand"], gene_pair_input_df["Ligand Location"],
            gene_pair_input_df["Receptor"], gene_pair_input_df["Receptor Location"],  gene_pair_input_df["Direct"],  gene_pair_input_df["Conserved"]
        )
    ]
    
    # Create base interaction card
    interaction_card = gene_pair_input_df[["Interaction ID", "LR Pair Card", "Interaction Type", "PMID", "is_mouse_specific", "Direct", "Conserved"]] # "KEGG Pathway", "Perplexity" "PROGENy Pathway", "Cancer-related", "Disease Type", "Disease",
### For the latest DB, skip for now
    # interaction_card["Perplexity"] = interaction_card["Perplexity"].str.replace('size=30', 'size=80')

    pop_up_info_lim = pop_up_info[
        ["Approved symbol", "Alias symbol", "Previous symbol", "Date symbol changed",  'ensembl_gene_id', 'OMIM', 'Open Targets Platform', "JensenLab DISEASES", 'Human Protein Atlas', "AmiGO", "PAN-GO"]
    ].drop_duplicates(subset="Approved symbol", keep="first")
    
    def format_symbol_aliases(old_symbol, aliases):
        parts = [p for p in (old_symbol, aliases) if p != "N/A"]
        return f"{', '.join(parts)}" if parts else aliases

    pop_up_info_lim['Other Symbols'] = pop_up_info_lim.apply(
        lambda row: format_symbol_aliases(row["Previous symbol"], row["Alias symbol"]),
        axis=1
    )
    
    ### Ligand cards
    #"Ligand RGD ID", 
    ligand_card = gene_pair_input_df[["LR Pair Card", "Ligand", "Ligand Name", "Ligand HGNC ID", "Ligand MGI ID", "Ligand Location", "is_mouse_specific"]].merge(
        pop_up_info_lim, how='left', left_on='Ligand', right_on='Approved symbol'
    ).drop_duplicates(subset='LR Pair Card', keep="first").drop(columns=["Approved symbol"])
    # create expression atlas link
    def create_exp_atlas_link(row):
        if row['is_mouse_specific']:
            base_link = convert_ensm_url_exp(row, "ensembl_gene_id")
        else:
            base_link = row["Human Protein Atlas"]
        return base_link
    ### For latest DB, use the mouse mapping for both LR if absent in Human
    # Apply the mapping to 'ensembl_gene_id'
    # ligand_card['ensembl_gene_id'] = ligand_card.apply(
    #     lambda row: mapping_mouse_ens.get(row['Ligand MGI ID'], row['ensembl_gene_id'])
    #     if pd.isna(row['ensembl_gene_id']) or str(row['ensembl_gene_id']).strip() == '' else row['ensembl_gene_id'],
    #     axis=1
    # )
    ligand_card['ensembl_gene_id'] = ligand_card.apply(
    lambda row: mapping_mouse_ens.get(row['Ligand MGI ID'], row['ensembl_gene_id'])
    if pd.notna(row['Ligand MGI ID']) else row['ensembl_gene_id'],
    axis=1
)
    # actual expression atlas creation
    ligand_card['Human Protein Atlas'] = ligand_card.apply(create_exp_atlas_link, axis=1)
    
    # Apply the mapping to 'Other Symbols'
    # ligand_card['Other Symbols'] = ligand_card.apply(
    #     lambda row: mapping_mouse_aliases.get(row['Ligand MGI ID'], row['Other Symbols'])
    #     if pd.isna(row['Other Symbols']) or str(row['Other Symbols']).strip() == '' else row['Other Symbols'],
    #     axis=1
    # )
    ligand_card['Other Symbols'] = ligand_card.apply(
    lambda row: mapping_mouse_aliases.get(row['Ligand MGI ID'], row['ensembl_gene_id'])
    if pd.notna(row['Ligand MGI ID']) else row['Other Symbols'],
    axis=1
)
    
    ligand_card['Other Symbols'] = ligand_card['Other Symbols'].apply(
    lambda x: "N/A" if pd.isna(x) or (isinstance(x, str) and x.strip().lower() == 'nan') else x
)

    ligand_card_1 = ligand_card[["LR Pair Card", "Ligand Name", "Other Symbols", "Ligand Location", "is_mouse_specific"]]
    
    ligand_card_2 = ligand_card[["LR Pair Card", "Ligand HGNC ID", "JensenLab DISEASES", "OMIM", 'Open Targets Platform', 'Human Protein Atlas', "AmiGO", "PAN-GO", "Ligand", "is_mouse_specific", "Ligand MGI ID", "ensembl_gene_id"]]
    
    # Create GeneCards - modify based on mouse-specific status
    def create_ligand_genecards_link(row):
        if row['is_mouse_specific']:
            base_link = convert_mgi_url(row, "Ligand", "Ligand MGI ID")
            if base_link:
                # Add mouse-specific indicator (example modification)
                return base_link.replace('target="_blank">', 'target="_blank" class="mouse-specific">')
            return base_link
        else:
            return convert_hgnc_url(row, "Ligand", "Ligand HGNC ID")
    
    ligand_card_2["GeneCards"] = ligand_card_2.apply(create_ligand_genecards_link, axis=1)

    # # Create a gene ontology for mouse
    def create_ligand_GO_link(row):
        if row['is_mouse_specific']:
            base_link = convert_mgi_GO_url(row, "Ligand MGI ID")
        else:
            base_link = row["PAN-GO"]
        return base_link
    
    ligand_card_2["PAN-GO"] = ligand_card_2.apply(create_ligand_GO_link, axis=1)

    
    # Apply conditional formatting for other ligand card elements
    def create_disease_relevance_link(row):
        if row['is_mouse_specific']:
            # For mouse-specific, you might want to modify the query or add a note
            base_link = convert_symbol_url_disease(row["Ligand"])
            # Example: modify the link text or add additional information
            if base_link:
                return base_link.replace('MalaCards', 'MalaCards (Mouse-specific)')
            return base_link
        else:
            return convert_symbol_url_disease(row["Ligand"])
    
    ligand_card_2["Disease relevance"] = ligand_card_2.apply(create_disease_relevance_link, axis=1)
    ligand_card_2["Human Cell Atlas"] = ligand_card["Ligand"].apply(convert_symbol_url_exp)
    ligand_card_2["GEPIA"] = ligand_card["Ligand"].apply(convert_symbol_url_exp_GEPIA)
    # add allen brain to gepia slot
    def create_brain_atlas_link(row):
        if row['is_mouse_specific']:
            base_link = convert_symbol_url_allenBrain(row, "Ligand")
        else:
            base_link = row["GEPIA"]
        return base_link
    ligand_card_2['GEPIA'] = ligand_card_2.apply(create_brain_atlas_link, axis=1)
    
    ligand_card_2["Expression Profile"] = ligand_card_2["Ligand HGNC ID"].apply(convert_hgnc_url_exp)
    def create_MCA_link(row):
        if row['is_mouse_specific']:
            base_link = convert_symbol_url_MCA(row, "Ligand")
        else:
            base_link = row["Expression Profile"]
        return base_link
    ligand_card_2['Expression Profile'] = ligand_card_2.apply(create_MCA_link, axis=1)
    ligand_card_2["HGNC Gene Group"] = ligand_card_2['Ligand HGNC ID'].map(ligand_mapping).fillna("none")
    
    icon_html_card = '<i class="fa-solid fa-arrow-up-right-from-square" style="margin-left:4px;"></i></a>'
    for col in ["Ligand HGNC ID"]:
        ligand_card_2[col] = ligand_card_2[col].str.replace(
            "</a>", icon_html_card, regex=False
        )
     # Apply the NCBI mapping to 'HGNC Gene Group'
    # ligand_card_2['HGNC Gene Group'] = ligand_card_2.apply(
    #     lambda row: mapping_mouse_ncbi.get(row['Ligand MGI ID'], row['HGNC Gene Group'])
    #     if pd.isna(row['HGNC Gene Group']) or str(row['HGNC Gene Group']).strip() == 'unknown' else row['HGNC Gene Group'],
    #     axis=1
    # )
    ligand_card_2['HGNC Gene Group'] = ligand_card_2.apply(
    lambda row: (
        mapping_mouse_ncbi.get(row['Ligand MGI ID'], row['HGNC Gene Group'])
        if pd.notna(row['Ligand MGI ID']) else row['HGNC Gene Group']
    ),
    axis=1
)

    # Set to "unknown" if still missing or blank after mapping
    ligand_card_2['HGNC Gene Group'] = ligand_card_2['HGNC Gene Group'].apply(
        lambda x: 'unknown' if pd.isna(x) or str(x).strip() == '' else x
    )
    def create_GEO_link(row):
        if row['is_mouse_specific']:
            base_link = convert_GEO_url(row, "HGNC Gene Group")
        else:
            base_link = row["Human Cell Atlas"]
        return base_link
    
    ligand_card_2["Human Cell Atlas"] = ligand_card_2.apply(create_GEO_link, axis=1)

    def create_ligand_ncbi_link(row):
        if row['is_mouse_specific']:
            base_link = convert_ncbi_url(row, "Ligand", "HGNC Gene Group")
            if base_link:
                # Add mouse-specific indicator (example modification)
                return base_link.replace('target="_blank">', 'target="_blank" class="mouse-specific">')
            return base_link
        else:
            return row["HGNC Gene Group"]

    ligand_card_2["HGNC Gene Group"] = ligand_card_2.apply(create_ligand_ncbi_link, axis=1)
    
    # Apply the uniprot mapping to 'AmiGO'
    # ligand_card_2['AmiGO'] = ligand_card_2.apply(
    #     lambda row: mapping_mouse_uniprot.get(row['Ligand MGI ID'], row['AmiGO'])
    #     if pd.isna(row['AmiGO']) or str(row['AmiGO']).strip() == 'unknown' else row['AmiGO'],
    #     axis=1
    # )
    ligand_card_2['AmiGO'] = ligand_card_2.apply(
    lambda row: mapping_mouse_uniprot.get(row['Ligand MGI ID'], row['AmiGO'])
    if pd.notna(row['Ligand MGI ID']) else row['AmiGO'],
    axis=1
)

    def create_quickGO_link(row):
        if row['is_mouse_specific']:
            base_link = convert_quickGO_url(row, "AmiGO")
        else:
            base_link = row["AmiGO"]
        return base_link
    
    ligand_card_2["AmiGO"] = ligand_card_2.apply(create_quickGO_link, axis=1)
    
    ligand_card_2 = ligand_card_2[["LR Pair Card", "Ligand HGNC ID", "GeneCards", "HGNC Gene Group", "Disease relevance", "Human Cell Atlas","GEPIA", "OMIM", "JensenLab DISEASES", "Expression Profile", 'Open Targets Platform', 'Human Protein Atlas', "AmiGO", "PAN-GO", "Ligand", "is_mouse_specific", "ensembl_gene_id"]]

    ###HERE
    
    #ligand_card_2 = create_conditional_dataframes(ligand_card_2)
    
    ### Receptor cards (similar modifications)
    # "Receptor RGD ID"
    receptor_card = gene_pair_input_df[["LR Pair Card", "Receptor", "Receptor Name", "Receptor HGNC ID", "Receptor MGI ID", "Receptor Location", "is_mouse_specific"]].merge(
        pop_up_info_lim, how='left', left_on='Receptor', right_on='Approved symbol'
    ).drop_duplicates(subset='LR Pair Card', keep="first").drop(columns=["Approved symbol"])

    # Apply the mapping to 'ensembl_gene_id'
    # receptor_card['ensembl_gene_id'] = receptor_card.apply(
    #     lambda row: mapping_mouse_ens.get(row['Receptor MGI ID'], row['ensembl_gene_id'])
    #     if pd.isna(row['ensembl_gene_id']) or str(row['ensembl_gene_id']).strip() == '' else row['ensembl_gene_id'],
    #     axis=1
    # )
    receptor_card['ensembl_gene_id'] = receptor_card.apply(
    lambda row: mapping_mouse_ens.get(row['Receptor MGI ID'], row['ensembl_gene_id'])
    if pd.notna(row['Receptor MGI ID']) else row['ensembl_gene_id'],
    axis=1
    )
    # actual expression atlas creation
    receptor_card['Human Protein Atlas'] = receptor_card.apply(create_exp_atlas_link, axis=1)
    # Apply the mapping to 'Other Symbols'
    # receptor_card['Other Symbols'] = receptor_card.apply(
    #     lambda row: mapping_mouse_aliases.get(row['Receptor MGI ID'], row['Other Symbols'])
    #     if pd.isna(row['Other Symbols']) or str(row['Other Symbols']).strip() == '' else row['Other Symbols'],
    #     axis=1
    # )
    receptor_card['Other Symbols'] = receptor_card.apply(
    lambda row: mapping_mouse_aliases.get(row['Receptor MGI ID'], row['ensembl_gene_id'])
    if pd.notna(row['Receptor MGI ID']) else row['Other Symbols'],
    axis=1
)
    receptor_card['Other Symbols'] = receptor_card['Other Symbols'].apply(
    lambda x: "N/A" if pd.isna(x) or (isinstance(x, str) and x.strip().lower() == 'nan') else x
)

    receptor_card_1 = receptor_card[["LR Pair Card", "Receptor Name", "Other Symbols", "Receptor Location", "is_mouse_specific"]]
    receptor_card_2 = receptor_card[["LR Pair Card", "Receptor HGNC ID", "JensenLab DISEASES", "OMIM", 'Open Targets Platform', 'Human Protein Atlas', "AmiGO", "PAN-GO", "Receptor", "is_mouse_specific", "Receptor MGI ID", "ensembl_gene_id"]]
    
    # Create GeneCards for receptors
    def create_receptor_genecards_link(row):
        if row['is_mouse_specific']:
            base_link = convert_mgi_url(row, "Receptor", "Receptor MGI ID")
            if base_link:
                return base_link.replace('target="_blank">', 'target="_blank" class="mouse-specific">')
            return base_link
        else:
            return convert_hgnc_url(row, "Receptor", "Receptor HGNC ID")
    
    receptor_card_2["GeneCards"] = receptor_card_2.apply(create_receptor_genecards_link, axis=1)
    
    def create_receptor_GO_link(row):
        if row['is_mouse_specific']:
            base_link = convert_mgi_GO_url(row, "Receptor MGI ID")
        else:
            base_link = row["PAN-GO"]
        return base_link
    
    receptor_card_2["PAN-GO"] = receptor_card_2.apply(create_receptor_GO_link, axis=1)
    def create_receptor_disease_relevance_link(row):
        if row['is_mouse_specific']:
            base_link = convert_symbol_url_disease(row["Receptor"])
            if base_link:
                return base_link.replace('MalaCards', 'MalaCards (Mouse-specific)')
            return base_link
        else:
            return convert_symbol_url_disease(row["Receptor"])
    
    receptor_card_2["Disease relevance"] = receptor_card_2.apply(create_receptor_disease_relevance_link, axis=1)
    receptor_card_2["Human Cell Atlas"] = receptor_card["Receptor"].apply(convert_symbol_url_exp)
    receptor_card_2["GEPIA"] = receptor_card["Receptor"].apply(convert_symbol_url_exp_GEPIA)
        # add allen brain to gepia slot
    def create_brain_atlas_link_receptor(row):
        if row['is_mouse_specific']:
            base_link = convert_symbol_url_allenBrain(row, "Receptor")
        else:
            base_link = row["GEPIA"]
        return base_link
    receptor_card_2['GEPIA'] = receptor_card_2.apply(create_brain_atlas_link_receptor, axis=1)
    receptor_card_2["Expression Profile"] = receptor_card_2["Receptor HGNC ID"].apply(convert_hgnc_url_exp)
    def create_MCA_link_receptor(row):
        if row['is_mouse_specific']:
            base_link = convert_symbol_url_MCA(row, "Receptor")
        else:
            base_link = row["Expression Profile"]
        return base_link
    receptor_card_2['Expression Profile'] = receptor_card_2.apply(create_MCA_link_receptor, axis=1)
    receptor_card_2["HGNC Gene Group"] = receptor_card_2['Receptor HGNC ID'].map(receptor_mapping).fillna("none")
    
    for col in ["Receptor HGNC ID"]:
        receptor_card_2[col] = receptor_card_2[col].str.replace(
            "</a>", icon_html_card, regex=False
        )
    # Apply the NCBI mapping to 'HGNC Gene Group'
    # receptor_card_2['HGNC Gene Group'] = receptor_card_2.apply(
    #     lambda row: mapping_mouse_ncbi.get(row['Receptor MGI ID'], row['HGNC Gene Group'])
    #     if pd.isna(row['HGNC Gene Group']) or str(row['HGNC Gene Group']).strip() == 'unknown' else row['HGNC Gene Group'],
    #     axis=1
    # )
    receptor_card_2['HGNC Gene Group'] = receptor_card_2.apply(
    lambda row: (
        mapping_mouse_ncbi.get(row['Receptor MGI ID'], row['HGNC Gene Group'])
        if pd.notna(row['Receptor MGI ID']) else row['HGNC Gene Group']
        ),
        axis=1
    )

    # Set to "unknown" if still missing or blank after mapping
    receptor_card_2['HGNC Gene Group'] = receptor_card_2['HGNC Gene Group'].apply(
        lambda x: 'unknown' if pd.isna(x) or str(x).strip() == '' else x
    )
    receptor_card_2["Human Cell Atlas"] = receptor_card_2.apply(create_GEO_link, axis=1)

    def create_ligand_ncbi_link_receptor(row):
        if row['is_mouse_specific']:
            base_link = convert_ncbi_url(row, "Receptor", "HGNC Gene Group")
            if base_link:
                # Add mouse-specific indicator (example modification)
                return base_link.replace('target="_blank">', 'target="_blank" class="mouse-specific">')
            return base_link
        else:
            return row["HGNC Gene Group"]

    receptor_card_2["HGNC Gene Group"] = receptor_card_2.apply(create_ligand_ncbi_link_receptor, axis=1)
 # Apply the uniprot mapping to 'AmiGO'
    # receptor_card_2['AmiGO'] = receptor_card_2.apply(
    #     lambda row: mapping_mouse_uniprot.get(row['Receptor MGI ID'], row['AmiGO'])
    #     if pd.isna(row['AmiGO']) or str(row['AmiGO']).strip() == 'unknown' else row['AmiGO'],
    #     axis=1
    # )
    receptor_card_2['AmiGO'] = receptor_card_2.apply(
    lambda row: mapping_mouse_uniprot.get(row['Receptor MGI ID'], row['AmiGO'])
    if pd.notna(row['Receptor MGI ID']) else row['AmiGO'],
    axis=1
)
    receptor_card_2["AmiGO"] = receptor_card_2.apply(create_quickGO_link, axis=1)
    

    
    receptor_card_2 = receptor_card_2[["LR Pair Card", "Receptor HGNC ID",  "GeneCards", "HGNC Gene Group", "Disease relevance", "Human Cell Atlas","GEPIA", "OMIM", 'Open Targets Platform', "JensenLab DISEASES", "Expression Profile",'Human Protein Atlas', "AmiGO", "PAN-GO", "Receptor", "is_mouse_specific", "ensembl_gene_id"]]
    
    # Rename and update HGNC columns
    ligand_card_2 = ligand_card_2.rename(columns={"Ligand HGNC ID": "HGNC Gene Symbol Report"})
    # ligand_card_2["HGNC Gene Symbol Report"] = ligand_card_2.apply(
    #     lambda row: update_link_text_with_symbol(row["HGNC Gene Symbol Report"], row["Ligand"]),
    #     axis=1
    # )
    def create_ensemblcards_link(row):
        if row['is_mouse_specific']:
            base_link = convert_mgi_ensembl_url(row, "Ligand", "ensembl_gene_id")
            if base_link:
                return base_link.replace('target="_blank">', 'target="_blank" class="mouse-specific">')
            return base_link
        else:
            return update_link_text_with_symbol(row["HGNC Gene Symbol Report"], row["Ligand"])
    
    ligand_card_2["HGNC Gene Symbol Report"] = ligand_card_2.apply(create_ensemblcards_link, axis=1)
    
    # Add perplexity with conditional formatting
    def create_ligand_perplexity_link(row):
        if row['is_mouse_specific']:
            # Modify query for mouse-specific ligands
            symbol = row["Ligand"]
            label = "Perplexity (LLM - Mouse-specific)"
            query = f"What diseases is {symbol} implicated in as shown in mouse models?"
            encoded_query = query.replace(" ", "%20")
            return f'<a href="https://www.perplexity.ai/search?q={encoded_query}" target="_blank">{label}</a>'
        else:
            return create_url_basic(row["Ligand"])
    
    ligand_card_2['Perplexity'] = ligand_card_2.apply(create_ligand_perplexity_link, axis=1)
    ligand_card_2 = ligand_card_2.drop(columns=["Ligand", "is_mouse_specific"])
    
    receptor_card_2 = receptor_card_2.rename(columns={"Receptor HGNC ID": "HGNC Gene Symbol Report"})
    def create_ensemblcards_link_receptor(row):
        if row['is_mouse_specific']:
            base_link = convert_mgi_ensembl_url(row, "Receptor", "ensembl_gene_id")
            if base_link:
                return base_link.replace('target="_blank">', 'target="_blank" class="mouse-specific">')
            return base_link
        else:
            return update_link_text_with_symbol(row["HGNC Gene Symbol Report"], row["Receptor"])
    
    receptor_card_2["HGNC Gene Symbol Report"] = receptor_card_2.apply(create_ensemblcards_link_receptor, axis=1)
    # Add perplexity for receptors
    def create_receptor_perplexity_link(row):
        if row['is_mouse_specific']:
            symbol = row["Receptor"]
            label = "Perplexity (LLM - Mouse-specific)"
            query = f"What diseases is {symbol} implicated in as shown in mouse models?"
            encoded_query = query.replace(" ", "%20")
            return f'<a href="https://www.perplexity.ai/search?q={encoded_query}" target="_blank">{label}</a>'
        else:
            return create_url_basic(row["Receptor"])
    
    receptor_card_2['Perplexity'] = receptor_card_2.apply(create_receptor_perplexity_link, axis=1)
    #receptor_card_2 = create_conditional_dataframes(receptor_card_2)
    

    receptor_card_2 = receptor_card_2.drop(columns=["Receptor", "is_mouse_specific"])
    # Remove is_mouse_specific from final output dataframes (except interaction_card if you want to use it in templates)
    ligand_card_1 = ligand_card_1.drop(columns=["is_mouse_specific"])
    receptor_card_1 = receptor_card_1.drop(columns=["is_mouse_specific"])
    
    return interaction_card, ligand_card_1, ligand_card_2, receptor_card_1, receptor_card_2



def convert_pair_url(df_pairs):
        """Converts 'LR Pair Card' column to HTML links with interaction IDs."""
        df_pairs = df_pairs.copy()
        
            # Extract clean Interaction ID from HTML
        df_pairs["Clean Interaction ID"] = df_pairs["Interaction ID"].apply(
                lambda x: re.search(r'(CDB\d+)</a>', x).group(1)
                if isinstance(x, str) and re.search(r'(CDB\d+)</a>', x)
                else None
            )
        # Generate new HTML anchor tag (with LR_Pair only
        def make_link(row):
            lrpair = row["LR Pair Card"]
            if pd.notna(lrpair) and lrpair.strip():
                encoded_lrpair = lrpair.replace(" ", "-")
                lrpair_dash = lrpair.replace(" ", " — ")
                
                # Determine subdirectory: 'mouse/' if any lowercase, else 'human/'
                subdir = "mouse/" if any(c.islower() for c in lrpair) else "human/"
                
                return (
                    f'<a href="{site_url}cards/'
                    f'{subdir}{encoded_lrpair}.html" target="_blank" '
                    f'title="Open {lrpair} card">'
                    f'{lrpair_dash}</a>'
                )
            return ""
        
        df_pairs["LR Pair Card"] = df_pairs.apply(make_link, axis=1)
        
        return df_pairs

# Find the generate_combined_html_files function and update these specific parts:

def generate_combined_html_files(
    gene_pair_keywords_df, # Corresponds to gene_pair000 (with '—' in LR Pair)
    template,
    interaction_card_df,
    ligand_card_1_df,
    receptor_card_1_df,
    ligand_card_2_df,
    receptor_card_2_df,
    pubmed_data_df,
    gene_pair_main_df, # Corresponds to gene_pair0 (with spaces in LR Pair)
    output_dir, 
    site
):
    """
    Generate combined HTML pages with PMID details on top and card details at the bottom
    for each LR Pair Card.
    """
    
    # --- Create a map of all pages for navigation ---
    cleaned_main_df = gene_pair_main_df.copy()
    
    def extract_clean_id(x):
        if isinstance(x, str):
            # Try HTML format first
            match = re.search(r'(CDB\d+)</a>', x)
            if match:
                return match.group(1)
            # Try plain text format
            match = re.search(r'(CDB\d+)', x)
            if match:
                return match.group(1)
        return None
    
    cleaned_main_df["Clean Interaction ID"] = cleaned_main_df["Interaction ID"].apply(extract_clean_id)
    cleaned_main_df = cleaned_main_df.dropna(subset=["Clean Interaction ID"])

    df_for_rendering_and_navigation = cleaned_main_df.sort_values(
        by="Clean Interaction ID",
        key=lambda series: series.str[3:].astype(int)
    ).reset_index(drop=True)

    rendered_pages = []

    def create_conditional_dataframes(dataframe, is_mouse_column='is_mouse_specific'):
        """
        Split dataframe into mouse-specific and regular dataframes with different column names
        """
        mouse_specific_df = dataframe[dataframe[is_mouse_column] == True].copy()
        regular_df = dataframe[dataframe[is_mouse_column] == False].copy()
        
        if not mouse_specific_df.empty and not regular_df.empty:
            combined_df = pd.concat([regular_df, mouse_specific_df], ignore_index=True)
        elif not mouse_specific_df.empty:
            combined_df = mouse_specific_df
        else:
            combined_df = regular_df
            
        return combined_df

    # --- Helper function to clean 'nan' values in dictionaries ---
    def clean_nan_values_in_dict(data_dict):
        cleaned_dict = {}
        for key, value in data_dict.items():
            if isinstance(value, str):
                # Aggressively clean and check
                stripped_lower_value = value.strip().lower()
                if stripped_lower_value == "nan" or stripped_lower_value == "null" or stripped_lower_value == "":
                    cleaned_dict[key] = "" # Replace with empty string
                else:
                    cleaned_dict[key] = value.strip() # Keep the original value, but trim it
            elif pd.isna(value): # Catches numpy/pandas NaNs (float type)
                cleaned_dict[key] = "" # Replace NaN float with empty string
            else:
                cleaned_dict[key] = value # Keep non-string, non-NaN values as is
        return cleaned_dict

    # --- First pass: Collect all content and metadata ---
    for idx, row in gene_pair_keywords_df.iterrows():
        lr_pair_name_hyphen = row["LR Pair Card"]  # e.g., "VEGFA—KDR"
        pmids_str = row["PMID"]
        
        lr_pair_name_space = lr_pair_name_hyphen.replace("—", " ")
        value1, value2 = lr_pair_name_space.split()
        #keywords = row["Relevance Keywords"]
        # --- Filter card dataframes for the current pair ---
        current_interaction_card = interaction_card_df[interaction_card_df['LR Pair Card'] == lr_pair_name_space].copy()
        current_ligand_card_1 = ligand_card_1_df[ligand_card_1_df['LR Pair Card'] == lr_pair_name_space].copy()
        current_receptor_card_1 = receptor_card_1_df[receptor_card_1_df['LR Pair Card'] == lr_pair_name_space].copy()
        current_ligand_card_2 = ligand_card_2_df[ligand_card_2_df['LR Pair Card'] == lr_pair_name_space].copy()
        current_receptor_card_2 = receptor_card_2_df[receptor_card_2_df['LR Pair Card'] == lr_pair_name_space].copy()

        if current_interaction_card.empty:
            print(f"[SKIP] No card data found for: {lr_pair_name_space}")
            continue
        
        # Determine if the current pair is mouse-specific
        is_current_pair_mouse_specific = current_interaction_card['is_mouse_specific'].iloc[0]
        
        # Apply conditional dataframe transformations based on is_mouse_specific
        if is_current_pair_mouse_specific:
            current_ligand_card_2 = current_ligand_card_2.rename(columns={
                "GeneCards": "MGI",
                "Ligand Name": "Mouse Ligand Name",
                "HGNC Gene Symbol Report": "Ensembl Gene Symbol Report",
                "HGNC Gene Group" : "NCBI"
            })
            current_receptor_card_2 = current_receptor_card_2.rename(columns={
                "GeneCards": "MGI",
                "Receptor Name": "Receptor Ligand Name",
                "HGNC Gene Symbol Report": "Ensembl Gene Symbol Report",
                "HGNC Gene Group" : "NCBI"
            })

        # Extract and clean interaction ID from the actual table data
        table0_data = current_interaction_card.drop('LR Pair Card', axis=1).to_dict(orient='records')[0]
        raw_interaction_id_html = table0_data["Interaction ID"]
        
        match = re.search(r'(CDB\d+)</a>', raw_interaction_id_html)
        if not match:
            match = re.search(r'(CDB\d+)', raw_interaction_id_html)
        
        if not match:
            print(f"[SKIP] Could not extract clean interaction ID for: {lr_pair_name_space}")
            print(f"[DEBUG] Interaction ID content: {raw_interaction_id_html}")
            continue
        clean_interaction_id = match.group(1)
        
        subdir = "mouse" if any(c.islower() for c in lr_pair_name_space) else "human"

        # Make sure the directory exists before writing
        target_dir = os.path.join(output_dir, subdir)
        os.makedirs(target_dir, exist_ok=True)
        
        filename = f"{value1.strip()}-{value2.strip()}.html"
        output_file = os.path.join(target_dir, filename)
        
        # Now drop the 'is_mouse_specific' and 'Ligand'/'Receptor' columns from the individual pair DataFrames
        # since their purpose for conditional renaming has been served.
        # Ensure they exist before dropping, using errors='ignore' for safety.
        current_interaction_card = current_interaction_card.drop(columns=['is_mouse_specific'], errors='ignore')
        current_ligand_card_1 = current_ligand_card_1.drop(columns=['is_mouse_specific', 'LR Pair Card'], errors='ignore')
        current_receptor_card_1 = current_receptor_card_1.drop(columns=['is_mouse_specific', 'LR Pair Card'], errors='ignore')
        current_ligand_card_2 = current_ligand_card_2.drop(columns=['Ligand', 'is_mouse_specific', 'LR Pair Card'], errors='ignore')
        current_receptor_card_2 = current_receptor_card_2.drop(columns=['Receptor', 'is_mouse_specific', 'LR Pair Card'], errors='ignore')

        # Save for navigation
        rendered_pages.append({
            "interaction_id": clean_interaction_id,
            "lr_pair_name_space": lr_pair_name_space,
            "value1": value1,
            "value2": value2,
            #"keywords": keywords,
            "pmids_str": pmids_str,
            "table0_data": table0_data, # This is already the modified interaction card row
            "row1": current_ligand_card_1,
            "row2": current_receptor_card_1,
            "row3": current_ligand_card_2,
            "row4": current_receptor_card_2,
            "row5":current_interaction_card,
            "output_file": output_file,
            "filename": filename,
        })

    # --- Second pass: Render each page with correct navigation ---
    for i, page in enumerate(rendered_pages):
        prev_page_info = rendered_pages[i - 1] if i > 0 else None
        next_page_info = rendered_pages[i + 1] if i < len(rendered_pages) - 1 else None
        
        # --- Prepare PMID section ---
# Find this section in your code (around line 800-850 in the generate_combined_html_files function)
# and replace the PMID section preparation with this:

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
                
                # NEW: Extract species and original gene info from gene_pair_per_row
                species_row = gene_pair_per_row[
                    (gene_pair_per_row["LR Pair Card"] == page["lr_pair_name_space"]) & 
                    (gene_pair_per_row["PMID"] == pmid)
                ]
                
                if not species_row.empty:
                    ligand_species = species_row["lig_species"].values[0]
                    ligand_orig = species_row["ligand_orig"].values[0]
                    receptor_species = species_row["rec_species"].values[0]
                    receptor_orig = species_row["receptor_orig"].values[0]
                else:
                    ligand_species = "Unknown"
                    ligand_orig = "Unknown"
                    receptor_species = "Unknown"
                    receptor_orig = "Unknown"
            else:
                title = "No Title Found"
                abstract = "No Abstract Found"
                journal = "Journal Unknown"
                year = "Year Unknown"
                ligand_species = "Unknown"
                ligand_orig = "Unknown"
                receptor_species = "Unknown"
                receptor_orig = "Unknown"
            
            active_class = "active" if j == 0 else ""
            tab_headers.append(f'<button class="tablinks {active_class}" onclick="openTab(event, \'tab{pmid}\')">{pmid}</button>')
            tab_contents.append(f"""
            <div id="tab{pmid}" class="tabcontent {active_class}">
                <h2 strong>{title}</h2>
                <div style="margin-left: 10px;">{journal}, {year}; <a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">PubMed</a>, <strong>{ligand_species} {ligand_orig} — {receptor_species} {receptor_orig}</strong></div> 
                <div class="abstract-wrapper">
                    <div class="abstract-content" id="abstract-content-{pmid}"><strong>ABSTRACT:</strong> {abstract}</div>
                </div>
            </div>
            """)
        
        # --- Prepare other tables ---
        def get_table_data(df):
            return df.to_dict(orient='records')[0] if not df.empty else {}
        
        table1_data = get_table_data(page["row1"])
        table2_data = get_table_data(page["row2"])
        
        # Apply the cleaning function to table3_data and table4_data here
        table3_data_raw = get_table_data(page["row3"])
        table3_data = clean_nan_values_in_dict(table3_data_raw) 
        
        table4_data_raw = get_table_data(page["row4"])
        table4_data = clean_nan_values_in_dict(table4_data_raw) 

        table5_data_raw = get_table_data(page["row5"])
        table5_data = clean_nan_values_in_dict(table5_data_raw) 

        # --- Related pairs ---
        ligand_pairs_df = gene_pair_main_df[(gene_pair_main_df['Ligand'] == page["value1"]) &
                                             (gene_pair_main_df["LR Pair Card"] != page["lr_pair_name_space"])].copy()
        receptor_pairs_df = gene_pair_main_df[(gene_pair_main_df['Receptor'] == page["value2"]) &
                                               (gene_pair_main_df["LR Pair Card"] != page["lr_pair_name_space"])].copy()
        
        # Apply URL conversion early
        ligand_pairs_df = convert_pair_url(ligand_pairs_df)
        receptor_pairs_df = convert_pair_url(receptor_pairs_df)
        
        # Now generate the joined string
        ligand_pairs_str = ' ・ '.join([btn for btn in ligand_pairs_df["LR Pair Card"] if btn])
        receptor_pairs_str = ' ・ '.join([btn for btn in receptor_pairs_df["LR Pair Card"] if btn])

        
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
        
        # --- UPDATED: Final render with updated navigation URLs ---
        rendered_html = template.render(
            GENE_NAME=page["lr_pair_name_space"],
            #KEYWORDS=page["keywords"],
            TAB_HEADERS="".join(tab_headers),
            TAB_CONTENTS="".join(tab_contents),
            value1=page["value1"],
            value2=page["value2"],
            table0_data=page["table0_data"],
            table1_data=table1_data,
            table2_data=table2_data,
            table3_data=table3_data,  # Use the cleaned table3_data
            table4_data=table4_data,  # Use the cleaned table4_data
            table5_data=table5_data,  # Use the cleaned table5_data
            ligand_image=ligand_image,
            receptor_image=receptor_image,
            ligand_pairs=ligand_pairs_str,
            receptor_pairs=receptor_pairs_str,
            site = site,
            prev_page_info=(
                {
                    "interaction_id": prev_page_info["interaction_id"],
                    "url": f"{site}cards/"
                           f"{'mouse/' if any(c.islower() for c in prev_page_info['lr_pair_name_space']) else 'human/'}"
                           f"{prev_page_info['filename']}",
                    "lr_pair_name_space": prev_page_info["lr_pair_name_space"]
                }
                if prev_page_info else None
            ),
            next_page_info=(
                {
                    "interaction_id": next_page_info["interaction_id"],
                    "url": f"{site}cards/"
                           f"{'mouse/' if any(c.islower() for c in next_page_info['lr_pair_name_space']) else 'human/'}"
                           f"{next_page_info['filename']}",
                    "lr_pair_name_space": next_page_info["lr_pair_name_space"]
                }
                if next_page_info else None
            )
        )

        
        # --- Save HTML file ---
        with open(page["output_file"], "w", encoding="utf-8") as f:
            f.write(rendered_html)


# --- Main Execution Block ---
if __name__ == "__main__":
    if test:
        # Define test genes - these should be in the 'space' format for gene_pair0
        # and will be converted to '—' for gene_pair000 internally.
        # Convert test_genes to '—' format for filtering gene_pair000
        #test_genes_hyphen = [gene.replace(" ", "—") for gene in test_genes]
        # Filter gene_pair0 for the test genes to be used in prepare_card_dataframes
        # This gene_pair_input should have space-separated LR pairs
        gene_pair_input = gene_pair0_copy[gene_pair0_copy["LR Pair Card"].isin(test_genes)]
        # Filter gene_pair000 for the test genes (which are in '—' format)
        # gene_pair_keywords_filtered = gene_pair000[gene_pair000["LR Pair Card"].isin(test_genes_hyphen)]
        gene_pair_keywords_filtered = gene_pair000[gene_pair000["LR Pair Card"].isin(test_genes)]
    
    else:
        gene_pair_input = gene_pair0_copy
        # Filter gene_pair000 for the test genes (which are in '—' format)
        gene_pair_keywords_filtered = gene_pair000
        print("Making all cards")

    # Prepare card-specific dataframes - NOW PASSING mouse_interaction_ids
    interaction_card, ligand_card_1, ligand_card_2, receptor_card_1, receptor_card_2 = prepare_card_dataframes(
        gene_pair_input, 
        mouse_interaction_ids=mouse_interaction_ids
    )

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
        output_dir=OUTPUT_DIR,
        site = site_url
    )
    print(f"Generated combined HTML files in: {OUTPUT_DIR}")
