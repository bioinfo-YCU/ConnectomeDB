## Function to prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the database qmds
import sys, os
from itables import init_notebook_mode
import pandas as pd
from itables import show
from itables import options
from IPython.display import HTML, display
import numpy as np
import fetchGSheet 
import warnings
import urllib.parse
import re

# Suppress SettingWithCopyWarning
warnings.simplefilter("ignore", category=UserWarning)

site_url = "https://connectomedb.org/"
# Other vertebrates
species_list = [
    "mmusculus", "rnorvegicus", "drerio", "ptroglodytes", "ggallus", "sscrofa", "btaurus", 
    "clfamiliaris", "ecaballus", "oarambouillet",
    "cjacchus", "mmulatta", "xtropicalis"
]

# Select only the relevant columns from pop_up_info
cols_to_keep = cols_to_keep = list(range(0, 30)) 
# Step 3: Load file using only the desired columns
df = pd.read_table("data/HGNC_gene_info_full.tsv", usecols=cols_to_keep)
pop_up_info = pd.read_table("data/HGNC_gene_info_full.tsv")
pop_up_info = pop_up_info.rename(columns={"hgnc_id": "HGNC ID", 
                                          "name": "Approved name",
                                          "symbol": "Approved symbol",
                                          #"rgd_id": "RGD ID",
                                          #"mgd_id": "MGI ID", 
                                          "alias_symbol": "Alias symbol", # add to table
                                          "prev_symbol": "Previous symbol", # add to table
                                          "date_symbol_changed": "Date symbol changed"
                                         })

# Keep only first MGI/RGD ID
#pop_up_info["MGI ID"] = pop_up_info["MGI ID"].str.split("|").str[0]
#pop_up_info["RGD ID"] = pop_up_info["RGD ID"].str.split("|").str[0]

pop_up_info["Alias symbol"] = pop_up_info["Alias symbol"].apply(
    lambda x: "N/A" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x
)

pop_up_info["Previous symbol"] = pop_up_info["Previous symbol"].apply(
    lambda x: "N/A" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x
)

# Replace "|" with ", "
pop_up_info["Alias symbol"] = [value.replace("|", ", ") for value in pop_up_info["Alias symbol"]]
pop_up_info["Previous symbol"] = [value.replace("|", ", ") for value in pop_up_info["Previous symbol"]]

pop_up_info["Date symbol changed"] = pop_up_info["Date symbol changed"].apply(
    lambda x: "N/A" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x
)


pop_up_info_lim = pop_up_info[["HGNC ID", "Approved name", "Alias symbol", # "MGI ID", "RGD ID"
                               "Approved symbol", "Previous symbol"]] # rm "Approved symbol" for now
pop_up_info_lim = pop_up_info_lim.drop_duplicates(subset="HGNC ID", keep="first")

# Drop columns where all values are NA in gene_pair
gene_pair = fetchGSheet.gene_pair_human.dropna(axis=1, how='all')
gene_pair = gene_pair[gene_pair['LR_pair_orig'] != '']
# for now set source count as triplicates
sourceCount = len(gene_pair[['LR_pair_orig']])

### KEEP ALL AS OF LATEST input datatable
# for now, keep only the following columns
# gene_pair = gene_pair[['LR pair', 'Ligand', 'Ligand.HGNC', 'Receptor', 'Receptor.HGNC',
#                        'perplexity link', 'PMID', 'binding location', 
#                        'bind in trans?', 'bidirectional signalling?',
#                        'interaction type', 'original source']]

gene_pair = gene_pair.dropna(subset=['LR_pair_orig'])

# some PMIDs kick in with "," so replace
gene_pair["PMID"] = [value.replace(",", "") for value in gene_pair["PMID"]]
gene_pair = gene_pair.dropna(subset=['PMID'])

### NO NEED FOR MAPPING AS OF LATEST input datatable
# Mapping for replacements
# mapping = dict(zip(fetchGSheet.src_info['original source'], fetchGSheet.src_info['shortname']))
# # Replace values in the column based on the mapping
# gene_pair['original source'] = gene_pair['original source'].replace(mapping)

gene_pair.columns = gene_pair.columns.str.strip()
gene_pair[['Ligand', 'Receptor']] = gene_pair['LR Pair Card'].str.split(' ', n=1, expand=True)

## add Ligand/Receptor Location
# def dedup_locations(loc_str):
#     # Split, strip, deduplicate, and sort
#     parts = [loc.strip() for loc in loc_str.split(',') if loc.strip()]
#     unique_sorted = sorted(set(parts), key=str.lower)  # case-insensitive sort
#     return unique_sorted

# def generate_LocToolTip(row, geneloc, loc_col):
#     ligand = row[loc_col]
#     original_locations = [loc.strip() for loc in row["location"].split(',')]
#     original_sources = [src.strip() for src in row["source"].split(',')]

#     # Get deduplicated locations
#     unique_locations = dedup_locations(row["location"])

#     if len(unique_locations) == 1:
#         # Single tooltip case
#         location = unique_locations[0]
#         matching_rows = geneloc[(geneloc[loc_col] == ligand) & (geneloc["location"].str.contains(location))]
#         all_sources = matching_rows["source"].unique()
#         sources_str = ", ".join(sorted(set(all_sources)))
#         return f'<span title="based on {sources_str}">{location}</span>'
#     else:
#         # Multiple tooltips — find each (ligand, location) match in original df
#         spans = []
#         for loc in unique_locations:
#             matching_rows = geneloc[
#                 (geneloc[loc_col] == ligand) &
#                 (geneloc["location"].str.contains(loc))
#             ]
#             all_sources = matching_rows["source"].unique()
#             sources_str = ", ".join(sorted(set(all_sources)))
#             spans.append(f'<span title="based on {sources_str}">{loc}</span>')
#         return ", ".join(spans)

# remove duplicate val per cell and prioritize "secreted"
def dedup_locations(loc_str):
    # Split and strip
    parts = [loc.strip() for loc in loc_str.split(',') if loc.strip()]
    # Deduplicate while preserving order
    seen = dict.fromkeys(parts)
    unique_ordered = list(seen.keys())
    # If "secreted" exists, move it to the front
    if "secreted" in seen:
        unique_ordered = ["secreted"] + [loc for loc in unique_ordered if loc != "secreted"]
    return ", ".join(unique_ordered)


# Group the original loc_info by Ligand
ligand_loc = fetchGSheet.ligand_loc.dropna(axis=1, how='all')
# make location short
ligand_loc.loc[ligand_loc["location"].str.contains("secreted", na=False), "location"] = "secreted"

grouped = ligand_loc.groupby("Ligand").agg({
    "location": lambda x: ', '.join(x),
    "source": lambda x: ', '.join(x)
}).reset_index()

# Generate tooltips
grouped = ligand_loc.groupby("Ligand").agg({
    "location": lambda x: dedup_locations(', '.join(x))
}).reset_index()

# create dict
mapping_loc = dict(zip(grouped['Ligand'], grouped['location'])) 
gene_pair['Ligand location'] = gene_pair['Ligand'].replace(mapping_loc)


# Group the original loc_info by Receptor
receptor_loc = fetchGSheet.receptor_loc.dropna(axis=1, how='all')
# make location short
receptor_loc.loc[receptor_loc["location"].str.contains("secreted", na=False), "location"] = "secreted"
grouped = receptor_loc.groupby("Receptor").agg({
    "location": lambda x: ', '.join(x),
    "source": lambda x: ', '.join(x)
}).reset_index()

# Generate tooltips
# grouped["Receptor location"] = grouped.apply(lambda row: generate_LocToolTip(row, receptor_loc,loc_col="Receptor"), axis=1)
grouped = receptor_loc.groupby("Receptor").agg({
    "location": lambda x: dedup_locations(', '.join(x))
}).reset_index()

# create dict
mapping_loc = dict(zip(grouped['Receptor'], grouped['location'])) 
gene_pair['Receptor location'] = gene_pair['Receptor'].replace(mapping_loc)


# Set missing mappings to 'unknown'
gene_pair.loc[gene_pair['Ligand location'] == gene_pair['Ligand'], 'Ligand location'] = 'unknown'
gene_pair.loc[gene_pair['Receptor location'] == gene_pair['Receptor'], 'Receptor location'] = 'unknown'
# Set "n/a" to unknown
gene_pair['Ligand location'] = [value.replace("n/a", "unknown") for value in gene_pair['Ligand location']]
gene_pair['Receptor location'] = [value.replace("n/a", "unknown") for value in gene_pair['Receptor location']]

# Fetch HGNC IDs from the dataset
hgnc_id = [col for col in gene_pair.columns if "HGNC" in col]
hgnc_id = pd.concat([gene_pair[col] for col in hgnc_id]).unique()

gene_pair['Human LR Pair'] = np.where(
    gene_pair['Human evidence'] == "not conserved", 
    "no human ortholog",                                  
    gene_pair['Homo sapiens_ligand'] + " " + gene_pair['Homo sapiens_receptor'] 
)


# Rename columns for better clarity
gene_pair = gene_pair.rename(columns={
    "LR_pair_orig": "LR Pair",
    "HGNC ligand": "Ligand HGNC ID",
    "HGNC receptor": "Receptor HGNC ID",
    "ENSEMBL ligand": "Ligand ENSEMBL ID",
    "ENSEMBL receptor": "Receptor ENSEMBL ID",
    # "perplexity link": "Perplexity", # will be replaced with actual link later
    # "original source": "Database Source",
    "Ligand location": "Ligand Location",
    "Receptor location": "Receptor Location",
    # "binding location": "Binding Location",
    # "bind in trans?" : "Trans-binding", 
    # "bidirectional signalling?": "Bidirectional Signalling",
    # "interaction type" : "Interaction Type"
    #"PMID": "PMID support" # was PMID support
})
gene_pair = gene_pair.drop(columns=["Ligand", "Receptor"])
gene_pair = gene_pair.rename(columns={"Homo sapiens_ligand": "Ligand", "Homo sapiens_receptor": "Receptor"})
# Merge gene_pair with pop_up_info_lim for Ligand(L)
gene_pair = gene_pair.merge(pop_up_info_lim, how='left', left_on='Ligand HGNC ID', right_on='HGNC ID')

gene_pair = gene_pair.rename(columns={"Approved name": "Ligand Name", 
                                     #"MGI ID": "Ligand MGI ID", # NOT APPLIED YET BUT should be taken from Ensembl BioMart
                                     #"RGD ID": "Ligand RGD ID", # NOT APPLIED YET BUT should be taken from Ensembl BioMart
                                      "Alias symbol": "Ligand Aliases",
                                      "Previous symbol": "Ligand Old symbol",
                                     },
                            )
gene_pair = gene_pair.drop(columns=["HGNC ID", "Approved symbol"])


gene_pair = gene_pair.merge(pop_up_info_lim, how='left', left_on='Receptor HGNC ID', right_on='HGNC ID')

gene_pair = gene_pair.rename(columns={"Approved name": "Receptor Name",
             #"MGI ID": "Receptor MGI ID",
             #"RGD ID": "Receptor RGD ID",
                                      "Alias symbol": "Receptor Aliases",
                                      "Previous symbol": "Receptor Old symbol",}
                            )

### AS OF LATEST DB skip pathways for now (code in addPathwayDiseaseAnnot_temp.py)

# Add new columns where all Ligand Symbol & Aliases and Receptor Symbol & Aliases merged in one column
def format_symbol_aliases(symbol, old_symbol, aliases):
    """
    Formats symbol, old symbols, and aliases.
    If the final formatted string would be empty after considering N/A values
    and empty inputs, it returns "mouse-specific".
    Otherwise, it formats based on the presence of old_symbol and aliases,
    removing unnecessary parentheses or commas, following the structure:
    "Symbol (Old Symbol, Aliases)" if both exist.
    """
    # Normalize inputs to empty strings if they are None/NaN or just whitespace
    symbol_str = str(symbol).strip()
    old_symbol_str = str(old_symbol).strip()
    aliases_str = str(aliases).strip()

    # Filter out values that are empty strings or "N/A" for old_symbol and aliases
    parts_for_join = []
    if old_symbol_str and old_symbol_str != "N/A":
        parts_for_join.append(old_symbol_str)
    if aliases_str and aliases_str != "N/A":
        parts_for_join.append(aliases_str)

    # Construct the preliminary result based on your original logic:
    # "symbol (old_symbol, aliases)" if parts_for_join is not empty, else "symbol"
    if parts_for_join:
        prelim_result = f"{symbol_str} ({', '.join(parts_for_join)})"
    else:
        prelim_result = symbol_str # Just the symbol if no old_symbol or aliases

    return prelim_result

# This is crucial for consistent handling by the function before processing "N/A".
gene_pair['Ligand'] = gene_pair['Ligand'].fillna('')
gene_pair['Ligand Old symbol'] = gene_pair['Ligand Old symbol'].fillna('')
gene_pair['Ligand Aliases'] = gene_pair['Ligand Aliases'].fillna('')


# to later check which ligand-receptor pairs are non-human
def is_mouse_specific(name):
    if not isinstance(name, str):
        return False
    name = name.strip()  # remove leading/trailing spaces
    return any(c.islower() for c in name[1:])

gene_pair = gene_pair.drop(columns=["HGNC ID"])

gene_pair['Ligand Symbols'] = gene_pair.apply(
    lambda row: "no human ortholog" if is_mouse_specific(row['Ligand']) 
               else format_symbol_aliases(row['Ligand'], row['Ligand Old symbol'], row['Ligand Aliases']),
    axis=1
)


# This is crucial for consistent handling by the function before processing "N/A".
gene_pair['Receptor'] = gene_pair['Receptor'].fillna('')
gene_pair['Receptor Old symbol'] = gene_pair['Receptor Old symbol'].fillna('')
gene_pair['Receptor Aliases'] = gene_pair['Receptor Aliases'].fillna('')

gene_pair['Receptor Symbols'] = gene_pair.apply(
    lambda row: "no human ortholog" if is_mouse_specific(row['Receptor']) 
                else format_symbol_aliases(row['Receptor'], row['Receptor Old symbol'], row['Receptor Aliases']),
    axis=1
)

### tooltips 
gene_pair["Ligand Symbols"] = [
    f'<span title="{aliases}">{aliases}</span>'
    for aliases in gene_pair["Ligand Symbols"]
]
gene_pair["Receptor Symbols"] = [
    f'<span title="{aliases}">{aliases}</span>'
    for aliases in gene_pair["Receptor Symbols"]
]

# might be used later just save info for now (for mouse cards)
grab_mouse_info = gene_pair["LR Pair Card"][gene_pair["Human evidence"].isin(["absent in human", "not conserved"])]
grab_mouse_info = grab_mouse_info.unique()
grab_mouse_info

# Add an empty A.I. summary column filled with None (just to save it's order
gene_pair['A.I. summary'] = None
#gene_pair = gene_pair.drop(columns=["Approved symbol_x", "Approved symbol_y"])

### For latest DB, skip (code saved as addOrth_temp.py)

# Add
first_columns=['LR Pair Card', 'Human LR Pair', 'Ligand', 'Receptor', 'Ligand Symbols', 'Receptor Symbols', 'Ligand Location', 'Receptor Location',	'Ligand HGNC ID', 'Receptor HGNC ID', 'A.I. summary', 'Human evidence', 'Ligand ENSEMBL ID', 'Receptor ENSEMBL ID'] # 'Database Source'

end_columns=['PMID', 'Pair_species', 'lig_species', 'rec_species', 'ligand_orig', 'receptor_orig']
gene_pair = gene_pair[first_columns + [col for col in gene_pair.columns if col not in first_columns + end_columns] + end_columns]
# gene_pair = gene_pair[first_columns + [col for col in gene_pair.columns if col not in first_columns]]

# number of unique vars (Human and Mouse both counted)

lrPairsCount = len(gene_pair["LR Pair Card"].unique())
ligand = gene_pair["LR Pair Card"].str.split(' ', expand=True)[0]
receptor = gene_pair["LR Pair Card"].str.split(' ', expand=True)[1]
ligandCount =len(ligand.unique())

receptorCount = len(receptor.unique())


### Remove from here for latest DB
# # Mouse Orthologue
# MouseLigandCount = len(gene_pair["Ligand MGI ID"].unique())

# MouseReceptorCount = len(gene_pair["Receptor MGI ID"].unique())

# # Rat Orthologue
# RatLigandCount = len(gene_pair["Ligand RGD ID"].unique())

# RatReceptorCount = len(gene_pair["Receptor RGD ID"].unique())

# gene_pair["PMID"] = [value.replace(" ", "") for value in gene_pair["PMID"]] # was'PMID support'


source = np.array(gene_pair["PMID"].unique())
source = source.astype(str)
source = ",".join(sorted(set(filter(lambda x: x.lower() != 'nan', source))))
# Split the string into individual elements, filter out empty strings, and get unique values
source = sorted(
    set(filter(lambda x: x.strip() and x.strip().lower() != 'nan', source.split(',')))
)
source = [value.replace(" ", "") for value in source]

# Function to join unique sorted values
agg_func = lambda x: ', '.join(sorted(set(map(str, x))))

# Group and aggregate all columns except 'LR Pair Card'
gene_pair = gene_pair.groupby('LR Pair Card').agg(agg_func).reset_index()
gene_pair = gene_pair[gene_pair['LR Pair Card'] != '']
# Identify rows where BOTH 'Ligand HGNC ID' and 'Receptor HGNC ID' are empty
# We check if the stripped string is empty, as fillna('') converts None/NaN to empty strings
has_hgnc_id = (gene_pair['Ligand HGNC ID'].astype(str).str.strip() != '') | \
              (gene_pair['Receptor HGNC ID'].astype(str).str.strip() != '')

# Separate the DataFrame into two parts (human-based cards and mouse based)
human_rows = gene_pair[~(gene_pair["Human evidence"].isin(["absent in human", "not conserved"]))]
mouse_rows = gene_pair[gene_pair["Human evidence"].isin(["absent in human", "not conserved"])]

# Concatenate the DataFrames: rows with IDs first, then rows without IDs
gene_pair = pd.concat([human_rows, mouse_rows]).reset_index(drop=True)
DBlength = len(gene_pair)

## Create interaction ids (SKIP -- use fixed one)
#gene_pair["Interaction ID"] = [f"CDB{str(i).zfill(5)}" for i in range(1, DBlength + 1)]

# Create mapping for Interaction ID
interaction_id_mapping = dict(zip(fetchGSheet.InteractionIDs['LR Pair Card'], fetchGSheet.InteractionIDs['Interaction ID']))

# Replace values in the column based on the mapping
gene_pair["Interaction ID"] = gene_pair["LR Pair Card"].replace(interaction_id_mapping)

# Sort by interaction ID
gene_pair = gene_pair.sort_values(
    by='Interaction ID',
    key=lambda col: col.str.split(':').str[1].astype(int)
)



# for creating PMIDs
gene_pair00 = gene_pair[['LR Pair Card', 'PMID']] # was "PMID support"

# create Perplexity link
def create_url_basic(perplexity_col):
    query = f"What is the primary evidence that {perplexity_col} bind-each-other-as-a-ligand-and-receptor-pair. Exclude reviews, uniprot, wiki, genecards, PIPS, iuphar as sources."
    encoded_query = query.replace(" ", "%20")
    return f"https://www.perplexity.ai/search?q={encoded_query}"
    
# Option 2 -- new query all together

# def generate_perplexity_link_pmid(row): 
#     query = f"What-is-the-biological-relevance-of-the-ligand-and-receptor-pair-{row['Human LR Pair']}-based-on-Pubmed-ID-{row['PMID']}"
#     return (
#          f'<a href="https://www.perplexity.ai/search?q={query}" target="_blank">'
#         f'<img src="https://img.icons8.com/?size=30&id=0NbBuNOxUwps&format=png&color=000000" alt="Perplexity AI" /></a>'
#     )

# cannot use perplexity logo
def generate_perplexity_link_pmid(row): 
    query = f"What-is-the-biological-relevance-of-the-ligand-and-receptor-pair-{row['Human LR Pair']}-based-on-Pubmed-ID-{row['PMID']}"
    return (
         f'<a href="https://www.perplexity.ai/search?q={query}" target="_blank" style="text-decoration: none;">&#128172;</a>'
    )

# Apply function to the DataFrame
gene_pair["A.I. summary"] = gene_pair.apply(generate_perplexity_link_pmid, axis=1)

# create URLs for the HGNC IDs

# create URLs for the Ensembl


gene_pair["Ligand HGNC ID"] = [
    '<a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{}" target="_blank">{}</a>'.format(ligand, ligand)
    for ligand in gene_pair["Ligand HGNC ID"]
]

# receptor
gene_pair["Receptor HGNC ID"] = [
    '<a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{}" target="_blank">{}</a>'.format(receptor, receptor)
    for receptor in gene_pair["Receptor HGNC ID"]
]


# ligand
gene_pair["Ligand ENSEMBL ID"] = [
    '<a href="http://www.ensembl.org/id/{}" target="_blank">{}</a>'.format(ligand, ligand)
    for ligand in gene_pair["Ligand ENSEMBL ID"]
]

# receptor
gene_pair["Receptor ENSEMBL ID"] = [
    '<a href="http://www.ensembl.org/id/{}" target="_blank">{}</a>'.format(receptor, receptor)
    for receptor in gene_pair["Receptor ENSEMBL ID"]
]



# Function to generate hyperlinks for the "PMID support" column
def generate_links_with_doi(df, gene_column, pmid_column, id_column):
    def create_link(gene, id_col, sources):
        # Replace spaces with "——" in the gene name for the link
        gene_name = gene.replace(" ", "—")
        
        if len(sources) == 1:
            source = sources[0]
            if source.startswith("https://www.biorxiv.org/content/"):
                # If the value starts with "https://doi.org/", use it as the hyperlink
                return f'<a href="{source}" target="_blank">BioRxiv</a>'
            else:
                # If it's a single PMID, hyperlink the PMID text
                return f'<a href="{site_url}cards/{gene}.html">{source}</a>'
        else:
            # If multiple PMIDs, show the count and hyperlink to the page
            return f'<a href="{site_url}cards/{gene}.html" target="_blank">{len(sources)} PMIDs</a>'

    # Process each row to generate the "PMID" column # was "PMID support"
    df["PMID"] = [
        create_link(
            gene=row[gene_column], 
            id_col = row[id_column],
            sources=[s.strip() for s in row[pmid_column].split(',') if s.strip()]
        )
        for _, row in df.iterrows()
    ]
    return df


# Generate the links for the "PMID" column # was "PMID support"
gene_pair = generate_links_with_doi(gene_pair, gene_column="Human LR Pair", 
                                    pmid_column="PMID", id_column= "Interaction ID")

# for disease type, cancer-related and top pathways, when missing say "ask Perplexity"


def generate_perplexity_kegglinks(
    df,
    pathway_col="KEGG Pathway",
    default_query_template="What-biological or other functional-pathways-is-the-ligand-receptor-{pair}-associated-with"
):
    def create_link(row):
        value = row.get(pathway_col, "")
        
        if pd.isna(value) or str(value).strip().lower() in ["nan", "none", "", "unknown"]:
            pair = row["Human LR Pair"]
            label = "ask Perplexity"
            query = default_query_template.format(pair=pair)
            encoded_query = urllib.parse.quote(query)
            return f'<a href="https://www.perplexity.ai/search?q={encoded_query}" target="_blank">{label}</a>'
        else:
            return value

    df[pathway_col] = df.apply(create_link, axis=1)
    return df

gene_pair = generate_perplexity_kegglinks(gene_pair, pathway_col="KEGG Pathway")
    
def generate_perplexity_links(df, pathway_col, default_query_template):
    def create_link(row):
        pathway_value = str(row[pathway_col]).strip().lower()
        pair = row["Human LR Pair"]
        
        if pd.isna(row[pathway_col]) or pathway_value in ["nan", "none", "", "unknown"]:
            label = "ask Perplexity"
            query = default_query_template.format(pair=pair)
            output =  f'<a href="https://www.perplexity.ai/search?q={query}" target="_blank">{label}</a>'
        else:
            label = row[pathway_col]
            query = f"What-is-the-role-of-the-ligand-and-receptor-pair-{pair}-in-{label}"
            output = f'{label} (see <a href="https://www.perplexity.ai/search?q={query}" target="_blank">evidence in Perplexity</a>)'
        
        return output
    
    df[pathway_col] = df.apply(create_link, axis=1)
    return df

### SKIP for latest DB

# gene_pair = generate_perplexity_links(
#     gene_pair,
#     pathway_col="PROGENy Pathway",
#     default_query_template="What-major signalling pathways-is-the-ligand-receptor-pair-{pair}-associated-with"
# )

# gene_pair = generate_perplexity_links(
#     gene_pair,
#     pathway_col="Disease Type",
#     default_query_template="What-disease types-is-the-ligand-receptor-pair-{pair}-associated-with"
# )
# # if it is a yes or no question
# def generate_perplexity_links_yesno(
#     df,
#     pathway_col="Cancer-related",
#     default_query_template="Is-the-{pair}-associated-with-cancer-and-if-so-which-ones"
# ):
#     def create_link(row):
#         pathway_value = str(row[pathway_col]).strip().lower()
#         pair = row["Human LR Pair"]
        
#         if pd.isna(row[pathway_col]) or pathway_value in ["nan", "none", "", "unknown"]:
#             label = "ask Perplexity"
#             query = default_query_template.format(pair=pair)
#             output = f'<a href="https://www.perplexity.ai/search?q={query}" target="_blank">{label}</a>'
#         else:
#             label = row[pathway_col]
#             query = f"Provide evidence to support the statement-{pair}-is-related-to-cancer-answer-is-{label}"
#             output = f'{label} (see <a href="https://www.perplexity.ai/search?q={query}" target="_blank">evidence in Perplexity</a>)'
#         return output
    
#     df[pathway_col] = df.apply(create_link, axis=1)
#     return df

# gene_pair = generate_perplexity_links_yesno(
#     gene_pair,
#     pathway_col="Cancer-related"
# )

# Add tooltip to name

def add_geneToolTip(species):
    def tooltip_html(symbol, name):
        return (
            f'<span class="tooltip">{symbol}'
            f'<span class="tooltiptext">{name}</span></span>'
        )

    gene_pair[species + " Ligand"] = [
        tooltip_html(ligand_symbol, ligand_name)
        for ligand_name, ligand_symbol in zip(gene_pair[species + " Ligand Name"], gene_pair[species + " Ligand"])
    ]
    gene_pair[species + " Receptor"] = [
        tooltip_html(receptor_symbol, receptor_name)
        for receptor_name, receptor_symbol in zip(gene_pair[species + " Receptor Name"], gene_pair[species + " Receptor"])
    ]

## Make the Human evidence consistent
gene_pair["Human evidence"] = np.where(
    gene_pair["Human evidence"].str.contains("DIRECT", na=False),
    "Direct",
    np.where(
        gene_pair["Human evidence"] == "CONSERVATION",
        "Conservation",
        gene_pair["Human evidence"]
    )
)

### Remove tooltip for name for each species for now as only zebrafish has the proper names ###     
# speciesPrime_list = ["Zebrafish"]
# # Loop through each species and update gene_pair
# for species in speciesPrime_list:
#    gene_pair = add_geneToolTip(species)

mouse_columns = ['Mouse Ligand', 'Mouse Receptor','Ligand MGI ID','Receptor MGI ID'] 
rat_columns = ['Rat Ligand','Rat Receptor','Ligand RGD ID','Receptor RGD ID']
zebrafish_columns = ['Zebrafish Ligand','Zebrafish Receptor','Ligand ZFIN ID','Receptor ZFIN ID']

# List of prefixes
prefixes = ("Chimpanzee", "Chicken", "Pig", "Cow", "Dog", "Horse", "Sheep", "Marmoset", "Macaque", "Frog")

# Get column names that start with any of the given prefixes
selected_columns = [col for col in gene_pair.columns if col.startswith(prefixes)]
# was "PMID support"
gene_pair0 = gene_pair[["Interaction ID"]+ first_columns+["PMID"]+["Ligand Name", "Receptor Name"]]
# gene_pair0 = gene_pair[["Interaction ID", "Human LR Pair", "Ligand", "Receptor",
#                        "Ligand Symbols", "Receptor Symbols", 
#                        "Ligand Location", "Receptor Location",
#                        "Ligand HGNC ID", "Receptor HGNC ID",
#                        "Perplexity", "PMID", 'Ligand Name','Receptor Name', 'KEGG Pathway', 'Cancer-related', 'Disease Type', 'Binding Location', 'Trans-binding', 'Bidirectional Signalling', 'Interaction Type', "PROGENy Pathway"] + mouse_columns + rat_columns]

gene_pair = gene_pair[["Interaction ID"]+ first_columns+["Ligand Name", "Receptor Name"]]

# gene_pair = gene_pair[["Interaction ID", "Human LR Pair", "Ligand", "Receptor",
#                        "Ligand Symbols", "Receptor Symbols", 
#                        "Ligand Location", "Receptor Location",
#                        "Ligand HGNC ID", "Receptor HGNC ID",
#                        "Perplexity", "PMID", 
#                        "Database Source", "Binding Location",
#                        "Trans-binding", "Bidirectional Signalling",
#                        "Interaction Type",'Ligand Name','Receptor Name'] + mouse_columns + rat_columns + zebrafish_columns + selected_columns]
# rm  "KEGG Pathway", "PROGENy Pathway", "Cancer-related", "Disease Type" for now

# Quick check if there is mouse-specific
gene_pair['Ligand'] = gene_pair.apply(
    lambda row: "no human ortholog" if is_mouse_specific(row['Ligand']) else row['Ligand'],
    axis=1
)
gene_pair['Receptor'] = gene_pair.apply(
    lambda row: "no human ortholog" if is_mouse_specific(row['Receptor']) else row['Receptor'],
    axis=1
)

# gene symbol
gene_pair["Ligand"] = [
    f'<span title="{ligand_name}">{ligand_symbol}</span>'
    for ligand_name, ligand_symbol in zip(gene_pair["Ligand Name"], 
                                              gene_pair["Ligand"])
]
# gene symbol
gene_pair["Receptor"] = [
    f'<span title="{receptor_name}">{receptor_symbol}</span>'
    for receptor_name, receptor_symbol in zip(gene_pair["Receptor Name"], 
                                              gene_pair["Receptor"])
]

### SKIP for latest DB
# tweak Database Source so when multiple should show multiple with tool tip
# gene_pair["Database Source"] = [
#     f'<span title="{orig_dbSource}">{("multiple (CDB2025 included)" if "connectomeDB2025" in orig_dbSource and "," in orig_dbSource else "multiple" if "," in orig_dbSource else orig_dbSource)}</span>'
#     for orig_dbSource in gene_pair["Database Source"]
# ]


# FOR NOW just make it as simple as em-dash/arrow
# def replace_spaces(row):
#     if 'secreted' in row['Ligand Location'].lower():
#         return row['Human LR Pair'].replace(" ", " <span style='font-size: 14px;'>○</span> <span style='font-size: 24px;'>⤚</span> ")
#     elif row['Ligand Location'] == 'unknown':
#         return row['Human LR Pair'].replace(" ", " <span style='font-size: 14px;'>○</span> <span style='font-size: 24px;'>⤚</span> ")
#     elif 'membrane' in row['Ligand Location'].lower():
#         return row['Human LR Pair'].replace(" ", " <span style='font-size: 24px;'>⤙</span> <span style='font-size: 24px;'>⤚</span> ")
#     else:
#         return row['Human LR Pair'].replace(" ", " \u2192 ")

# # Apply the function to the 'LR Pair' column
# gene_pair['Human LR Pair'] = gene_pair.apply(replace_spaces, axis=1)
# gene_pair["Human LR Pair"] = gene_pair["Human LR Pair"].str.replace(" ", "-")
gene_pair = gene_pair.drop(columns=["Ligand Name", "Receptor Name"])


# Create the links to the HTML cards
gene_pair["LR Pair Card"] = [
    f'<a href="{site_url}cards/{ "mouse" if evidence == "not conserved" else "human" }/{lrPairOrig.replace(" ","-")}.html" target="_blank">{lrPair}</a>'
    for lrPairOrig, lrPair, evidence in zip(gene_pair0["LR Pair Card"], gene_pair["LR Pair Card"], gene_pair["Human evidence"])
]


# Add tooltips to the column headers
gene_pair.columns = [
    f'<span title="Unique ConnectomeDB ID for each ligand–receptor pair">{col}</span>' if col == "Interaction ID" else
    f'<span title="Ligand-receptor Pair">{col}</span>' if col == "Human LR Pair" else
    f'<span title="HGNC gene symbol for the ligand">{col}</span>' if col == "Ligand" else
    f'<span title="HGNC gene symbol for the receptor">{col}</span>' if col == "Receptor" else
     f'<span title="Official gene symbol (aliases, old names)">{col}</span>' if col in ["Ligand Symbols", "Receptor Symbols"] else
    f'<span title="Click the icon below to run Perplexity on the Human LR pair">{col}</span>' if col == "A.I. summary" else
    f'<span title="Official Gene Symbol; Hover on symbols below to show gene names">{col}</span>' if col in ["Ligand", "Receptor"] else
    f'<span title="HGNC gene ID for the ligand (link to HGNC)">{col}</span>' if col == "Ligand HGNC ID" else
    
    f'<span title="HGNC gene ID for the receptor (link to HGNC)">{col}</span>' if col == "Receptor HGNC ID" else
    f'<span title="ENSEMBL gene ID for the ligand (link to ENSEMBL)">{col}</span>' if col  == "Ligand ENSEMBL ID" else
    f'<span title="ENSEMBL gene ID for the receptor (link to ENSEMBL)">{col}</span>' if col == "Receptor ENSEMBL ID" else
    f'<span title=" PubMed IDs (PMID) with Literature Evidence for LR Interaction. Click on the link for more details">{col}</span>' if col == "PMID" else
    f'<span title="Location based on the predicted subcellular localization of the human proteome">{col}</span>' if col in ["Ligand Location", "Receptor Location"] else
    f'<span title="Direct: experimentally verified; Conservation: inferred from orthology">{col}</span>' if col == "Human evidence" else
    f'<span title="Double-click header of {col} to reverse sort">{col}</span>'
    for col in gene_pair.columns
]

gene_pair = gene_pair.reset_index(drop=True)  # Remove the index

#######################################################################
# Identify the column(s) that contain '(PMID)' and temporarily remove for presubmission
pmid_cols = [col for col in gene_pair.columns if '(PMID)' in col]
gene_pair = gene_pair.drop(columns=pmid_cols)
#######################################################################

gene_pair000 = gene_pair.copy()

keywords_to_modify = ["Ligand", "Receptor"]
exclude_keywords = ["HGNC ID", "Location", "Human"]  # Columns containing this will not be modified

# Copy the original columns so we can modify only the first 10
new_columns = gene_pair000.columns.tolist()

# Modify only the first 10 columns
new_columns[:10] = [
    f'{col.split(">")[0]}">Human {col.split(">")[1]}</span>'
    if any(keyword in col for keyword in keywords_to_modify) and not any(exclude in col for exclude in exclude_keywords)
    else col
    for col in new_columns[:10]
]

# Assign the modified column names back to the DataFrame
gene_pair000.columns = new_columns
#######################################################################
### For latest DB, no need to limit columns
human_columns = [col for col in gene_pair000.columns]
#######################################################################
### For latest DB, no need to limit rows
#human_gene_pair = gene_pair.iloc[:, :-36]
# remove mouse specific ones from the datatable
evidence_cols = [col for col in gene_pair.columns if 'Human evidence' in col]
human_gene_pair = gene_pair[~(gene_pair[evidence_cols[0]] == "not conserved")]

### remove LR Pair Card and Just use Interaction ID
# Replace visible text with Interaction ID
def replace_link_text(interaction_id, link_html):
    if pd.isna(link_html):
        return interaction_id  # fallback if no link
    return re.sub(r'>(.*?)<', f'>{interaction_id}<', str(link_html), count=1)
interactionID_cols = [col for col in human_gene_pair.columns if 'Interaction ID' in col][0]
paircard_cols = [col for col in human_gene_pair.columns if 'LR Pair Card' in col][0]
# Apply the transformation
human_gene_pair[interactionID_cols] = human_gene_pair.apply(
    lambda row: replace_link_text(row[interactionID_cols], row[paircard_cols]),
    axis=1
)
# Drop the old LR Pair Card column
human_gene_pair = human_gene_pair.drop(columns=[paircard_cols])

# add number of mouse pair cards
numOfMouseOrth = len(gene_pair[evidence_cols][(gene_pair[evidence_cols[0]] == "not conserved")])