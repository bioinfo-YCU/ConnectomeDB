## Function to prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the database qmds
import sys, os
sys.path.append(os.path.abspath("src"))  
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


# Select only the relevant columns from pop_up_info

pop_up_info = fetchGSheet.pop_up_info.rename(columns={"Mouse genome informatics (MGI) ID": "MGI ID", "Rat genome database (RGD) ID": "RGD ID"})

pop_up_info_lim = pop_up_info[["Approved symbol", "Approved name", "MGI ID", "RGD ID"]]
pop_up_info_lim = pop_up_info_lim.drop_duplicates(subset="Approved symbol", keep="first")

# Drop columns where all values are NA in gene_pair
gene_pair = fetchGSheet.gene_pair.dropna(axis=1, how='all')
# Fetch species IDs from the dataset
hgnc_id = [col for col in gene_pair.columns if "HGNC ID" in col]
hgnc_id = pd.concat([gene_pair[col] for col in hgnc_id]).unique()

# Rename columns for better clarity
gene_pair = gene_pair.rename(columns={
    "Ligand receptor pair": "Human LR Pair",
    "Ligand gene symbol": "Ligand",
    "Receptor gene symbol": "Receptor",
    "Perplexity link": "Perplexity"
})

# Merge gene_pair with pop_up_info_lim for Ligand(L)
gene_pair = gene_pair.merge(pop_up_info_lim, how='left', left_on='Ligand', right_on='Approved symbol')

gene_pair = gene_pair.rename(columns={"Approved name": "Ligand name", 
                                     "MGI ID": "Ligand MGI ID",
                                     "RGD ID": "Ligand RGD ID"},
                            )

# Add MGI name
MGI_info = pd.read_csv("data/MGI_ID_biomart.csv")
gene_pair = gene_pair.merge(MGI_info, how='left', left_on='Ligand MGI ID', right_on='MGI ID')

# Add RGD name
RGD_info = pd.read_csv("data/RGD_ID_biomart.csv")
RGD_info['RGD ID'] = "RGD:" + RGD_info['RGD ID'].astype(str)
gene_pair = gene_pair.merge(RGD_info, how='left', left_on='Ligand RGD ID', right_on='RGD ID')

gene_pair = gene_pair.drop(columns=["RGD ID", "MGI ID"])

gene_pair = gene_pair.rename(columns={
                                     "MGI name": "Mouse Ligand", 
                                     "RGD name": "Rat Ligand"}
                            )

gene_pair = gene_pair.merge(pop_up_info_lim, how='left', left_on='Receptor', right_on='Approved symbol')

gene_pair = gene_pair.rename(columns={"Approved name": "Receptor name",
                                      "MGI ID": "Receptor MGI ID",
                                      "RGD ID": "Receptor RGD ID"}
                            )

# Add MGI name
gene_pair = gene_pair.merge(MGI_info, how='left', left_on='Receptor MGI ID', right_on='MGI ID')
gene_pair = gene_pair.merge(RGD_info, how='left', left_on='Receptor RGD ID', right_on='RGD ID')
gene_pair = gene_pair.drop(columns=["RGD ID", "MGI ID"])

gene_pair = gene_pair.rename(columns={
                                     "MGI name": "Mouse Receptor", 
                                     "RGD name": "Rat Receptor"}
                            )
gene_pair = gene_pair.drop(columns=["Approved symbol_x", "Approved symbol_y"])

# Drop columns where all values are NA in gene_pair
gene_pair = gene_pair.dropna(axis=1, how='all')

gene_pair = gene_pair.fillna(" ")
gene_pair = gene_pair[gene_pair['Human LR Pair'] != ' ']

if "PMID link" in gene_pair.columns:
    gene_pair = gene_pair.drop(columns=["PMID link"])

# Add
first_columns=['Human LR Pair', 'Ligand', 'Receptor', 'Source']

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

mouse_columns = [col for col in gene_pair.columns if "MGI" in col or "Mouse" in col]
rat_columns = [col for col in gene_pair.columns if "RGD" in col or "Rat" in col]

gene_pair0 = gene_pair[['Human LR Pair', 'Ligand', 'Receptor', 'Perplexity', 'PMID support',
       'Ligand HGNC ID', 'Ligand location', 'Receptor HGNC ID',
       'Receptor location', 'Ligand name', 'Receptor name'] + mouse_columns + rat_columns]

gene_pair = gene_pair[['Human LR Pair', 'Ligand', 'Receptor', 'Source', 'Perplexity', 'PMID support',
       'Ligand location', 'Receptor location', 'Ligand HGNC ID', 'Receptor HGNC ID',
        'Ligand name', 'Receptor name'] + mouse_columns + rat_columns + end_columns]
# gene symbol
gene_pair["Ligand"] = [
    f'<span title="{ligand_name}">{ligand_symbol}</span>'
    for ligand_name, ligand_symbol in zip(gene_pair["Ligand name"], gene_pair["Ligand"])
]

# gene symbol
gene_pair["Receptor"] = [
    f'<span title="{receptor_name}">{receptor_symbol}</span>'
    for receptor_name, receptor_symbol in zip(gene_pair["Receptor name"], gene_pair["Receptor"])
]

def replace_spaces(row):
    if row['Ligand location'] == 'secreted':
        return row['Human LR Pair'].replace(" ", " ○ <span style='font-size: 30px;'>⤚</span> ")
    elif row['Ligand location'] == 'plasma membrane':
        return row['Human LR Pair'].replace(" ", " <span style='font-size: 30px;'>⤙</span> <span style='font-size: 30px;'>⤚</span> ")
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
    f'<span title="Ligand Receptor Pair">{col}</span>' if col == "Human LR Pair" else
    f'<span title="Click the logo below to run Perplexity on the Human LR pair">{col}&nbsp;</span>' if col == "Perplexity" else
    f'<span title="Hover on symbols below to show gene names">{col}&nbsp;&nbsp;&nbsp;</span>' if col in ["Ligand", "Receptor"] else
    f'<span title="Click on HGNC IDs below for more details">{col}&nbsp;&nbsp;</span>' if col in ["Ligand HGNC ID", "Receptor HGNC ID"] else
    f'<span title="Click on the Pubmed IDs (PMID) below for more details">{col}</span>' if col == "PMID support" else
    f'<span title="Click on the Rat Genome Database(RGD) IDs below for more details">{col}</span>' if col in ["Ligand RGD ID", "Receptor RGD ID"] else
    f'<span title="Click on the Mouse Genome Informatics(MGI) IDs below for more details">{col}</span>' if col in ["Ligand MGI ID", "Receptor MGI ID"] else
    f'<span title="Location is defined as predicted subcellular location (HUMAN)">{col}</span>' if col in ["Ligand location", "Receptor location"] else
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

# Find columns with "Mouse" in the name
mouse_columns = [col for col in gene_pair.columns if "MGI" in col or "Mouse" in col]

# Filter rows where all "Mouse" columns are not " "
mouse_gene_pair = gene_pair000[(gene_pair000[mouse_columns].map(str.strip) != "").all(axis=1)]
# Dynamically identify columns containing "Ligand" and "Receptor" in their names 
# since it is now in span format

new_columns = mouse_gene_pair.columns.tolist()

new_columns = [
    col.replace("Mouse ", "").strip()
    if "Mouse Ligand" in col or "Mouse Receptor" in col
    else col
    for col in new_columns
]
mouse_gene_pair.columns = new_columns
        
ligand_col = [col for col in mouse_gene_pair.columns if "Ligand&nbsp;" in col][1]
receptor_col = [col for col in mouse_gene_pair.columns if "Receptor&nbsp;" in col][1]
ligand_location = [col for col in mouse_gene_pair.columns if "Ligand location" in col][0]
receptor_location = [col for col in mouse_gene_pair.columns if "Receptor location" in col][0]


# Combine columns into "Mouse LR Pair" with appropriate replacements
def format_lr_pair(row):
    if row[ligand_location] == 'secreted':
        return f"{row[ligand_col]} ○ <span style='font-size: 30px;'>⤚</span> {row[receptor_col]}"
    elif row[receptor_location] == 'plasma membrane':
        return f"{row[ligand_col]} <span style='font-size: 30px;'>⤙</span> <span style='font-size: 30px;'>⤚</span> {row[receptor_col]}"
    else:
        return f"{row[ligand_col]} \u2192 {row[receptor_col]}"

# Apply the function row-wise and assign to the new column using .loc
mouse_gene_pair1 = mouse_gene_pair.copy() 
mouse_gene_pair1.loc[:, "Mouse LR Pair"] = mouse_gene_pair1.apply(format_lr_pair, axis=1)
mouse_columns = [col for col in mouse_gene_pair1.columns if "MGI" in col]
# Reorder the DataFrame
new_order = ["Mouse LR Pair", ligand_col, receptor_col] + mouse_columns + human_columns
mouse_gene_pair1 = mouse_gene_pair1[new_order]
mouse_gene_pair1 = mouse_gene_pair1.reset_index(drop=True)  

## Limit to those with either Rat Ligand or Receptor
rat_columns = [col for col in gene_pair.columns if "RGD" in col or "Rat" in col]
# Filter rows where all "Rat" columns are not " "
rat_gene_pair = gene_pair000[(gene_pair000[rat_columns].map(str.strip) != "").all(axis=1)]


new_columns = rat_gene_pair.columns.tolist()

new_columns = [
    col.replace("Rat ", "").strip()
    if "Ligand" in col or "Receptor" in col
    else col
    for col in new_columns
]
rat_gene_pair.columns = new_columns

# Dynamically identify columns containing "Ligand" and "Receptor" in their names 
# since it is now in span format
ligand_col = [col for col in rat_gene_pair.columns if "Ligand&nbsp;" in col][1]
receptor_col = [col for col in rat_gene_pair.columns if "Receptor&nbsp;" in col][1]
ligand_location = [col for col in rat_gene_pair.columns if "Ligand location" in col][0]
receptor_location = [col for col in rat_gene_pair.columns if "Receptor location" in col][0]

def format_lr_pair(row):
    if row[ligand_location] == 'secreted':
        return f"{row[ligand_col]} ○ <span style='font-size: 30px;'>⤚</span> {row[receptor_col]}"
    elif row[receptor_location] == 'plasma membrane':
        return f"{row[ligand_col]} <span style='font-size: 30px;'>⤙</span> <span style='font-size: 30px;'>⤚</span> {row[receptor_col]}"
    else:
        return f"{row[ligand_col]} \u2192 {row[receptor_col]}"

rat_gene_pair1 = rat_gene_pair.copy() 
rat_gene_pair1.loc[:, "Rat LR Pair"] = rat_gene_pair1.apply(format_lr_pair, axis=1)
rat_columns = [col for col in rat_gene_pair1.columns if "RGD" in col]
# Reorder the DataFrame
new_order = ["Rat LR Pair", ligand_col, receptor_col] + rat_columns + human_columns
rat_gene_pair1 = rat_gene_pair1[new_order]
rat_gene_pair1 = rat_gene_pair1.reset_index(drop=True)  