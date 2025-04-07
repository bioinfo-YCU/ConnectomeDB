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

# Suppress SettingWithCopyWarning
warnings.simplefilter("ignore", category=UserWarning)


# Other vertebrates
species_list = [
    "ptroglodytes", "ggallus", "sscrofa", "btaurus", 
    "clfamiliaris", "ecaballus", "oarambouillet"
]

# Select only the relevant columns from pop_up_info

pop_up_info = fetchGSheet.pop_up_info.rename(columns={"Mouse genome informatics (MGI) ID": "MGI ID", 
                                                      "Rat genome database (RGD) ID": "RGD ID"})

pop_up_info_lim = pop_up_info[["Approved symbol", "Approved name", "MGI ID", "RGD ID"]]
pop_up_info_lim = pop_up_info_lim.drop_duplicates(subset="Approved symbol", keep="first")

# Drop columns where all values are NA in gene_pair
gene_pair = fetchGSheet.gene_pair.dropna(axis=1, how='all')

# Mapping for replacements
mapping = {
    'Ramilowski_2015_Literature_supported': "connectomeDB2015",
    'Noël et al. 2020 (ICELLNET)': "ICELLNET",
    'Hou et al. 2020 (connectomeDB2020)': "connectomeDB2020",
    'Efremova et al. 2020 (CellphoneDB)': "CellphoneDB",
    'Cabello-Aguilar et al. 2020 (SingleCellSignalR)': "SingleCellSignalR",
    'Baccin et al. 2020 (RNA-Magnet)': "RNA-Magnet",
    'ConnectomeDB2025 (this publication)': "connectomeDB2025 🆕"
}

# Replace values in the column based on the mapping
gene_pair['Source'] = gene_pair['Source'].replace(mapping)

# Fetch species IDs from the dataset
hgnc_id = [col for col in gene_pair.columns if "HGNC ID" in col]
hgnc_id = pd.concat([gene_pair[col] for col in hgnc_id]).unique()

# Rename columns for better clarity
gene_pair = gene_pair.rename(columns={
    "Ligand receptor pair": "Human LR Pair",
    "Ligand gene symbol": "Ligand",
    "Receptor gene symbol": "Receptor",
    "Perplexity link": "Perplexity",
    "Source": "Database Source"
})

# Merge gene_pair with pop_up_info_lim for Ligand(L)
gene_pair = gene_pair.merge(pop_up_info_lim, how='left', left_on='Ligand', right_on='Approved symbol')

gene_pair = gene_pair.rename(columns={"Approved name": "Ligand name", 
                                     "MGI ID": "Ligand MGI ID",
                                     "RGD ID": "Ligand RGD ID"},
                            )

# Add top pathway per pair
LR_pairs = gene_pair["Human LR Pair"].unique()
df= pd.read_csv("data/pathway_annotations_per_pair.csv")
df = df[df["interaction"].isin(LR_pairs)]
# Sort by absolute value of 'weight', descending (larger abs(weight) first)
df_sorted = df.reindex(df['weight'].abs().sort_values(ascending=False).index)
# Keep only the first occurrence for each unique 'interaction'
df_unique = df_sorted.drop_duplicates(subset='interaction', keep='first')
df = df_unique.reset_index(drop=True)
top_pathway_df = df[["interaction", "source"]]
top_pathway_df = top_pathway_df.rename(columns={
                                      "source": "Top Pathway"
})
gene_pair = gene_pair.merge(top_pathway_df, how='left', left_on='Human LR Pair', right_on='interaction')

# Add Disease Category per pair
df= pd.read_csv("data/diseaseType_per_pair.csv")
disease_df = df[df["interaction_x"].isin(LR_pairs)]

gene_pair = gene_pair.merge(disease_df, how='left', left_on='Human LR Pair', right_on='interaction_x')

# Add MGI annotation
MGI_info = pd.read_csv("data/MGI_ID_biomart.csv")
gene_pair = gene_pair.merge(MGI_info, how='left', left_on='Ligand MGI ID', right_on='MGI ID')

# Add RGD annotation
RGD_info = pd.read_csv("data/RGD_ID_biomart.csv")
RGD_info['RGD ID'] = "RGD:" + RGD_info['RGD ID'].astype(str)
gene_pair = gene_pair.merge(RGD_info, how='left', left_on='Ligand RGD ID', right_on='RGD ID')

# Add ZFIN id and symbol
ZFIN_info = pd.read_csv("data/ZFIN_ID_human_orthos.txt", sep="\t", skiprows=1)
ZFIN_info = ZFIN_info[['ZFIN ID', 'ZFIN Symbol', 'ZFIN Name', 'HGNC ID']]

ZFIN_info = ZFIN_info.dropna(subset=['HGNC ID'])
ZFIN_info = ZFIN_info.drop_duplicates(subset=['HGNC ID'])
ZFIN_info['HGNC ID'] = ZFIN_info['HGNC ID'].apply(lambda x: f'HGNC:{int(x)}')
gene_pair = gene_pair.merge(ZFIN_info, how='left', left_on='Ligand HGNC ID', right_on='HGNC ID')

gene_pair = gene_pair.drop(columns=["RGD ID", "MGI ID", "HGNC ID", "interaction", "interaction_x"])

gene_pair = gene_pair.rename(columns={
                                     "MGI name": "Mouse Ligand", 
                                     "RGD name": "Rat Ligand",
                                     "ZFIN ID": "Ligand ZFIN ID",
                                     "ZFIN Symbol": "Zebrafish Ligand",
                                     "ZFIN Name": "Zebrafish Ligand name"}
                            )

gene_pair = gene_pair.merge(pop_up_info_lim, how='left', left_on='Receptor', right_on='Approved symbol')
gene_pair = gene_pair.rename(columns={"Approved name": "Receptor name",
                                      "MGI ID": "Receptor MGI ID",
                                      "RGD ID": "Receptor RGD ID"}
                            )

# Add MGI name
gene_pair = gene_pair.merge(MGI_info, how='left', left_on='Receptor MGI ID', right_on='MGI ID')
gene_pair = gene_pair.merge(RGD_info, how='left', left_on='Receptor RGD ID', right_on='RGD ID')
gene_pair = gene_pair.merge(ZFIN_info, how='left', left_on='Receptor HGNC ID', right_on='HGNC ID')
gene_pair = gene_pair.drop(columns=["RGD ID", "MGI ID", "HGNC ID"])

gene_pair = gene_pair.rename(columns={
                                     "MGI name": "Mouse Receptor", 
                                     "RGD name": "Rat Receptor",
                                     "ZFIN ID": "Receptor ZFIN ID",
                                     "ZFIN Symbol": "Zebrafish Receptor",
                                     "ZFIN Name": "Zebrafish Receptor name"}
                            )

gene_pair = gene_pair.drop(columns=["Approved symbol_x", "Approved symbol_y"])

# Function to add species-specific species Enseml ID and symbol for all other species except for mouse, rat, and zebrafish
def appendOtherSpeciesInfo(species, origDF):
    species_name = {
    "ptroglodytes": "Chimpanzee",
    "ggallus": "Chicken",
    "sscrofa": "Pig",
    "btaurus": "Cow",
    "clfamiliaris": "Dog",
    "ecaballus": "Horse",
    "oarambouillet": "Sheep",
    }.get(species, "Unknown species")
    
    # Load species-specific data
    species_info = pd.read_csv(f"data/{species}_ID_biomart.csv")

    # Keep relevant columns
    species_info = species_info[[f"{species}_homolog_ensembl_gene", 
                                 f"{species}_homolog_associated_gene_name", 
                                 'hgnc_id']]

    # Remove rows where 'hgnc_id' is NaN and drop duplicates
    species_info = species_info.dropna(subset=['hgnc_id'])
    species_info = species_info.drop_duplicates(subset=['hgnc_id'])

    # Merge with ligand data
    origDF = origDF.merge(species_info, how='left', 
                           left_on='Ligand HGNC ID', right_on='hgnc_id')
    
    # Rename columns for ligand info
    origDF = origDF.rename(columns={
        f"{species}_homolog_associated_gene_name": f"{species_name} Ligand", 
        f"{species}_homolog_ensembl_gene": f"{species_name} Ligand Ensembl ID"
    })

    # Drop duplicate 'hgnc_id' column
    origDF = origDF.drop(columns=['hgnc_id'])

    # Merge with receptor data
    origDF = origDF.merge(species_info, how='left', 
                           left_on='Receptor HGNC ID', right_on='hgnc_id')

    # Rename columns for receptor info
    origDF = origDF.rename(columns={
        f"{species}_homolog_associated_gene_name": f"{species_name} Receptor", 
        f"{species}_homolog_ensembl_gene": f"{species_name} Receptor Ensembl ID"
    })

        # Drop duplicate 'hgnc_id' column
    origDF = origDF.drop(columns=['hgnc_id'])

    # Drop columns where all values are NaN
    origDF = origDF.dropna(axis=1, how='all')

    return origDF


# Loop through each species and update gene_pair
for species in species_list:
    gene_pair = appendOtherSpeciesInfo(species, gene_pair)

# Drop columns where all values are NA in gene_pair
gene_pair = gene_pair.dropna(axis=1, how='all')

gene_pair = gene_pair.fillna(" ")
gene_pair = gene_pair[gene_pair['Human LR Pair'] != ' ']

if "PMID link" in gene_pair.columns:
    gene_pair = gene_pair.drop(columns=["PMID link"])

# Add
first_columns=['Human LR Pair', 'Ligand', 'Receptor', 'Database Source']

end_columns=['HGNC L R', 'sanity check', 'curator', 'secondary source?']
gene_pair = gene_pair[first_columns + [col for col in gene_pair.columns if col not in first_columns + end_columns] + end_columns]


# number of unique vars

lrPairsCount = len(gene_pair["Human LR Pair"].unique())

ligandCount = len(gene_pair["Ligand"].unique())

receptorCount = len(gene_pair["Receptor"].unique())

# Mouse Orthologue
MouseLigandCount = len(gene_pair["Ligand MGI ID"].unique())

MouseReceptorCount = len(gene_pair["Receptor MGI ID"].unique())

# Rat Orthologue
RatLigandCount = len(gene_pair["Ligand RGD ID"].unique())

RatReceptorCount = len(gene_pair["Receptor RGD ID"].unique())

gene_pair["PMID support"] = [value.replace(" ", "") for value in gene_pair["PMID support"]]

source = np.array(gene_pair["PMID support"].unique())
source = source.astype(str)
source = ",".join(sorted(set(filter(lambda x: x.lower() != 'nan', source))))

# Split the string into individual elements, filter out empty strings, and get unique values
source = sorted(
    set(filter(lambda x: x.strip() and x.strip().lower() != 'nan', source.split(',')))
)
source = [value.replace(" ", "") for value in source]
sourceCount = len(source)

# for creating PMIDs
gene_pair00 = gene_pair[['Human LR Pair', 'PMID support']]

# create URLs for the HGNC IDs

# ligand
gene_pair["Ligand HGNC ID"] = [
    '<a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{}" target="_blank">{}</a>'.format(ligand, ligand)
    for ligand in gene_pair["Ligand HGNC ID"]
]

# receptor
gene_pair["Receptor HGNC ID"] = [
    '<a href="https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{}" target="_blank">{}</a>'.format(receptor, receptor)
    for receptor in gene_pair["Receptor HGNC ID"]
]

# Perplexity
gene_pair["Perplexity"] = [
    '<a href="{}" target="_blank"> <img src="https://img.icons8.com/?size=30&id=0NbBuNOxUwps&format=png&color=000000" alt="Perplexity AI" /></a>'.format(url)
    for url in gene_pair["Perplexity"]
]

# Function to generate hyperlinks for the "PMID support" column
# Function to generate hyperlinks for the "PMID support" column
def generate_links_with_doi(df, gene_column, pmid_column):
    def create_link(gene, sources):
        # Replace spaces with "——" in the gene name for the link
        gene_name = gene.replace(" ", "——")
        
        if len(sources) == 1:
            source = sources[0]
            if source.startswith("https://www.biorxiv.org/content/"):
                # If the value starts with "https://doi.org/", use it as the hyperlink
                return f'<a href="{source}" target="_blank">BioRxiv preprint</a>'
            else:
                # If it's a single PMID, hyperlink the PMID text
                return f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/pubmed/{gene_name}_pmid_details.html">{source}</a>'
        else:
            # If multiple PMIDs, show the count and hyperlink to the page
            return f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/pubmed/{gene_name}_pmid_details.html" target="_blank">{len(sources)} PMIDs</a>'

    # Process each row to generate the "PMID support" column
    df["PMID support"] = [
        create_link(
            gene=row[gene_column], 
            sources=[s.strip() for s in row[pmid_column].split(',') if s.strip()]
        )
        for _, row in df.iterrows()
    ]
    return df


# Generate the links for the "PMID support" column
gene_pair = generate_links_with_doi(gene_pair, gene_column="Human LR Pair", pmid_column="PMID support")

# for disease type, cancer-related and top pathways, explicitly say "not available"
gene_pair["Top Pathway"] = gene_pair["Top Pathway"].apply(
    lambda x: "unknown" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x
)
gene_pair["Disease Type"] = gene_pair["Disease Type"].apply(
    lambda x: "unknown" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x
)

gene_pair["Cancer-related"] = gene_pair["Cancer-related"].apply(
    lambda x: "unknown" if pd.isna(x) or str(x).strip().lower() in ["nan", "none", ""] else x
)

gene_pair["Ligand MGI ID"] = [
        f'<a href="https://www.informatics.jax.org/marker/{mouseOrth}" target="_blank">{mouseOrth}</a>' 
        if pd.notna(mouseOrth) and mouseOrth.strip() else "" 
        for mouseOrth in gene_pair["Ligand MGI ID"]
    ]

gene_pair["Receptor MGI ID"] = [
        f'<a href="https://www.informatics.jax.org/marker/{mouseOrth}" target="_blank">{mouseOrth}</a>' 
        if pd.notna(mouseOrth) and mouseOrth.strip() else "" 
        for mouseOrth in gene_pair["Receptor MGI ID"]
    ]

gene_pair["Ligand RGD ID"] = [
        f'<a href="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id={ratOrth.replace("RGD:", "")}" target="_blank">{ratOrth}</a>' 
        if pd.notna(ratOrth) and ratOrth.strip() else "" 
        for ratOrth in gene_pair["Ligand RGD ID"]
    ]

gene_pair["Receptor RGD ID"] = [
        f'<a href="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id={ratOrth.replace("RGD:", "")}" target="_blank">{ratOrth}</a>' 
        if pd.notna(ratOrth) and ratOrth.strip() else "" 
        for ratOrth in gene_pair["Receptor RGD ID"]
    ]
# Add tooltip to name

def add_geneToolTip(species):
    gene_pair[species+" Ligand"] = [
        f'<span title="{ligand_name}">{ligand_symbol}</span>'
        for ligand_name, ligand_symbol in zip(gene_pair[species+" Ligand name"], gene_pair[species+" Ligand"])
    ]
    gene_pair[species+" Receptor"] = [
        f'<span title="{receptor_name}">{receptor_symbol}</span>'
        for receptor_name, receptor_symbol in zip(gene_pair[species+" Receptor name"], gene_pair[species+" Receptor"])
    ]

### Remove tooltip for name for each species for now as only zebrafish has the proper names ###     
# speciesPrime_list = ["Zebrafish"]
# # Loop through each species and update gene_pair
# for species in speciesPrime_list:
#    gene_pair = add_geneToolTip(species)

mouse_columns = ['Mouse Ligand', 'Mouse Receptor','Ligand MGI ID','Receptor MGI ID'] 
rat_columns = ['Rat Ligand','Rat Receptor','Ligand RGD ID','Receptor RGD ID']
zebrafish_columns = ['Zebrafish Ligand','Zebrafish Receptor','Ligand ZFIN ID','Receptor ZFIN ID']

# List of prefixes
prefixes = ("Chimpanzee", "Chicken", "Pig", "Cow", "Dog", "Horse", "Sheep")

# Get column names that start with any of the given prefixes
selected_columns = [col for col in gene_pair.columns if col.startswith(prefixes)]

gene_pair0 = gene_pair[['Human LR Pair', 'Ligand', 'Receptor', 'Perplexity', 'PMID support',
       'Ligand HGNC ID', 'Ligand location', 'Receptor HGNC ID',
       'Receptor location', 'Ligand name', 'Receptor name', 'Top Pathway', 'Cancer-related', 'Disease Type'] + mouse_columns + rat_columns]

gene_pair = gene_pair[['Human LR Pair', 'Database Source', 'Ligand', 'Receptor', 'Perplexity', 'PMID support',
        'Ligand HGNC ID', 'Receptor HGNC ID', 'Ligand location', 'Receptor location',
        'Ligand name', 'Receptor name', 'Top Pathway', 'Cancer-related', 'Disease Type'] + mouse_columns + rat_columns + zebrafish_columns + end_columns + selected_columns]


# gene symbol
gene_pair["Ligand"] = [
    f'<span title="{ligand_name}">{ligand_symbol}</span>'
    for ligand_name, ligand_symbol in zip(gene_pair["Ligand name"], 
                                          gene_pair["Ligand"])
]

# gene symbol
gene_pair["Receptor"] = [
    f'<span title="{receptor_name}">{receptor_symbol}</span>'
    for receptor_name, receptor_symbol in zip(gene_pair["Receptor name"], 
                                              gene_pair["Receptor"])
]

def replace_spaces(row):
    if row['Ligand location'] == 'secreted':
        return row['Human LR Pair'].replace(" ", " <span style='font-size: 15px;'>○</span> <span style='font-size: 24px;'>⤚</span> ")
    elif row['Ligand location'] == '':
        return row['Human LR Pair'].replace(" ", " <span style='font-size: 15px;'>○</span> <span style='font-size: 24px;'>⤚</span> ")
    elif row['Ligand location'] == 'plasma membrane':
        return row['Human LR Pair'].replace(" ", " <span style='font-size: 24px;'>⤙</span> <span style='font-size: 24px;'>⤚</span> ")
    else:
        return row['Human LR Pair'].replace(" ", " \u2192 ")

# Apply the function to the 'LR Pair' column
gene_pair['Human LR Pair'] = gene_pair.apply(replace_spaces, axis=1)

gene_pair = gene_pair.drop(columns=["Ligand name", "Receptor name"])


# Create the links to the HTML cards
gene_pair["Human LR Pair"] = [
    f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/{lrPairOrig}.html">{lrPair}</a>'
    for lrPairOrig, lrPair in zip(gene_pair0["Human LR Pair"], gene_pair["Human LR Pair"])
]



# Add tooltips to the column headers
gene_pair.columns = [
    f'<span title=" Ligand-Receptor Interacting Pair, as described in Liu et al. (PMID: XXXXXX)">{col}</span>' if col == "Human LR Pair" else
    f'<span title="Click the logo below to run Perplexity on the Human LR pair">{col}&nbsp;</span>' if col == "Perplexity" else
    f'<span title=" Official Gene Symbol; Hover on symbols below to show gene names">{col}&nbsp;&nbsp;&nbsp;</span>' if col in ["Ligand", "Receptor"] else
    f'<span title="HUGO Gene Nomenclature Committee (HGNC) ID. Click on the link for more details">{col}&nbsp;&nbsp;</span>' if col in ["Ligand HGNC ID", "Receptor HGNC ID"] else
    f'<span title=" PubMed IDs (PMID) with Literature Evidence for LR Interaction. Click on the link for more details">{col}</span>' if col == "PMID support" else
    f'<span title="Rat Genome Database (RGD) ID. Click on the link for more details">{col}</span>' if col in ["Ligand RGD ID", "Receptor RGD ID"] else
    f'<span title="Mouse Genome Informatics (MGI) ID. Click on the link for more details">{col}</span>' if col in ["Ligand MGI ID", "Receptor MGI ID"]else
    f'<span title="Zebrafish Information Network (ZFIN) ID. Click on the link for more details">{col}</span>' if col in ["Ligand ZFIN ID", "Receptor ZFIN ID"] else
    f'<span title="Location based on the predicted subcellular localization of the human proteome, as described in Ramilowski et al. (PMID: 26198319)">{col}</span>' if col in ["Ligand location", "Receptor location"] else
    f'<span title="Double-click header of {col} to ensure all values are shown">{col}&nbsp;</span>'
    for col in gene_pair.columns
]

gene_pair = gene_pair.reset_index(drop=True)  # Remove the index
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
human_columns = [col for col in gene_pair000.columns][:10]