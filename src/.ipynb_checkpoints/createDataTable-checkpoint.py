from itables import init_notebook_mode
import pandas as pd
from itables import show
from itables import options
from IPython.display import HTML, display
import numpy as np
import fetchGSheet 

# Select only the relevant columns from pop_up_info

pop_up_info = fetchGSheet.pop_up_info.rename(columns={"Mouse genome informatics (MGI) ID": "Mouse (MGI) ID","Rat genome database (RGD) ID": "Rat (RGD) ID"})

pop_up_info_lim = pop_up_info[["Approved symbol", "Approved name", "Mouse (MGI) ID", "Rat (RGD) ID"]]
pop_up_info_lim = pop_up_info_lim.drop_duplicates(subset="Approved symbol", keep="first")

# Drop columns where all values are NA in gene_pair
gene_pair = fetchGSheet.gene_pair.dropna(axis=1, how='all')

# Rename columns for better clarity
gene_pair = gene_pair.rename(columns={
    "Ligand receptor pair": "LR Pair",
    "Ligand gene symbol": "Ligand",
    "Receptor gene symbol": "Receptor",
    "Perplexity link": "Perplexity"
})

# Merge gene_pair with pop_up_info_lim for Ligand(L)
gene_pair = gene_pair.merge(pop_up_info_lim, how='left', left_on='Ligand', right_on='Approved symbol')

gene_pair = gene_pair.rename(columns={"Approved name": "Ligand name", 
                                     "Mouse (MGI) ID": "Ligand Mouse (MGI) ID",
                                     "Rat (RGD) ID": "Ligand Rat (RGD) ID"},
                            )

gene_pair = gene_pair.merge(pop_up_info_lim, how='left', left_on='Receptor', right_on='Approved symbol')

gene_pair = gene_pair.rename(columns={"Approved name": "Receptor name",
                                      "Mouse (MGI) ID": "Receptor Mouse (MGI) ID",
                                      "Rat (RGD) ID": "Receptor Rat (RGD) ID"}
                            )

gene_pair = gene_pair.drop(columns=["Approved symbol_x", "Approved symbol_y"])

# Drop columns where all values are NA in gene_pair
gene_pair = gene_pair.dropna(axis=1, how='all')

gene_pair = gene_pair.fillna(" ")
gene_pair = gene_pair[gene_pair['LR Pair'] != ' ']

if "PMID link" in gene_pair.columns:
    gene_pair = gene_pair.drop(columns=["PMID link"])

# Add
first_columns=['LR Pair','Source', 'Ligand','Receptor', 'Perplexity']
end_columns=['HGNC L R','sanity check', 'curator','secondary source?']
gene_pair = gene_pair[first_columns + [col for col in gene_pair.columns if col not in first_columns + end_columns] + end_columns]

# number of unique vars

lrPairsCount = len(gene_pair["LR Pair"].unique())

ligandCount = len(gene_pair["Ligand"].unique())

receptorCount = len(gene_pair["Receptor"].unique())

# Mouse Orthologue
MouseLigandCount = len(gene_pair["Ligand Mouse (MGI) ID"].unique())

MouseReceptorCount = len(gene_pair["Receptor Mouse (MGI) ID"].unique())

# Rat Orthologue
RatLigandCount = len(gene_pair["Ligand Rat (RGD) ID"].unique())

RatReceptorCount = len(gene_pair["Receptor Rat (RGD) ID"].unique())

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
gene_pair00 = gene_pair[['LR Pair', 'PMID support']]

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
    '<a href="{}" target="_blank"> <img src="https://img.icons8.com/?size=35&id=0NbBuNOxUwps&format=png&color=000000" alt="Perplexity AI" /></a>'.format(url)
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
                return f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/pubmed/{gene_name}_pmid_details.html" target="_blank">{source}</a>'
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
gene_pair = generate_links_with_doi(gene_pair, gene_column="LR Pair", pmid_column="PMID support")

gene_pair["Ligand Mouse (MGI) ID"] = [
        f'<a href="https://www.informatics.jax.org/marker/{mouseOrth}" target="_blank">{mouseOrth}</a>' 
        if pd.notna(mouseOrth) and mouseOrth.strip() else "" 
        for mouseOrth in gene_pair["Ligand Mouse (MGI) ID"]
    ]

gene_pair["Receptor Mouse (MGI) ID"] = [
        f'<a href="https://www.informatics.jax.org/marker/{mouseOrth}" target="_blank">{mouseOrth}</a>' 
        if pd.notna(mouseOrth) and mouseOrth.strip() else "" 
        for mouseOrth in gene_pair["Receptor Mouse (MGI) ID"]
    ]

gene_pair["Ligand Rat (RGD) ID"] = [
        f'<a href="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id={ratOrth.replace("RGD:", "")}" target="_blank">{ratOrth}</a>' 
        if pd.notna(ratOrth) and ratOrth.strip() else "" 
        for ratOrth in gene_pair["Ligand Rat (RGD) ID"]
    ]

gene_pair["Receptor Rat (RGD) ID"] = [
        f'<a href="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id={ratOrth.replace("RGD:", "")}" target="_blank">{ratOrth}</a>' 
        if pd.notna(ratOrth) and ratOrth.strip() else "" 
        for ratOrth in gene_pair["Receptor Rat (RGD) ID"]
    ]

gene_pair0 = gene_pair[['LR Pair', 'Ligand', 'Receptor', 'Perplexity', 'PMID support',
       'Ligand HGNC ID', 'Ligand location', 'Receptor HGNC ID', 'Ligand Mouse (MGI) ID','Receptor Mouse (MGI) ID',
       'Receptor location', 'Ligand name', 'Receptor name', 'Ligand Rat (RGD) ID', 'Receptor Rat (RGD) ID']]


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
        return row['LR Pair'].replace(" ", " ○ <span style='font-size: 30px;'>⤚</span> ")
    elif row['Ligand location'] == 'plasma membrane':
        return row['LR Pair'].replace(" ", " <span style='font-size: 30px;'>⤙</span> <span style='font-size: 30px;'>⤚</span> ")
    else:
        return row['LR Pair'].replace(" ", " \u2192 ")

# Apply the function to the 'LR Pair' column
gene_pair['LR Pair'] = gene_pair.apply(replace_spaces, axis=1)

gene_pair["Ligand location"] = [
    '<span title="This is the tooltip for {}">{}</span>'.format(loc, loc)
    for loc in gene_pair["Ligand location"]
]

gene_pair["Receptor location"] = [
    '<span title="This is the tooltip for {}">{}</span>'.format(loc, loc)
    for loc in gene_pair["Receptor location"]
]


gene_pair = gene_pair.drop(columns=["Ligand name", "Receptor name"])


# Create the links to the HTML cards
gene_pair["LR Pair"] = [
    f'<a href="https://comp.med.yokohama-cu.ac.jp/collab/connectomeDB/cards/{lrPairOrig}.html" target="_blank">{lrPair}</a>'
    for lrPairOrig, lrPair in zip(gene_pair0["LR Pair"], gene_pair["LR Pair"])
]


# Add tooltips to the column headers
gene_pair.columns = [
    f'<span title="Ligand Receptor Pair">{col}</span>' if col == "LR Pair" else
    f'<span title="Click the logo below to run Perplexity on the LR pair">{col}&nbsp;</span>' if col == "Perplexity" else
    f'<span title="Hover on symbols below to show gene names">{col}&nbsp;&nbsp;&nbsp;</span>' if col in ["Ligand", "Receptor"] else
    f'<span title="Click on HGNC IDs below for more details">{col}&nbsp;&nbsp;</span>' if col in ["Ligand HGNC ID", "Receptor HGNC ID"] else
    f'<span title="Click on the Pubmed IDs (PMID) below for more details">{col}</span>' if col == "PMID support" else
    f'<span title="Click on the Rat Genome Database(RGD) IDs below for more details">{col}</span>' if col in ["Ligand Rat (RGD) ID", "Receptor Rat (RGD) ID"] else
    f'<span title="Click on the Mouse Genome Informatics(MGI) IDs below for more details">{col}</span>' if col in ["Ligand Mouse (MGI) ID", "Receptor Mouse (MGI) ID"] else
    f'<span title="Location is defined as subcellular location">{col}</span>' if col in ["Ligand location", "Receptor location"] else
    f'<span title="Double-click header of {col} to ensure all values are shown">{col}&nbsp;</span>'
    for col in gene_pair.columns
]

gene_pair = gene_pair.reset_index(drop=True)  # Remove the index

## Limit to those with either Mouse Ligand or Receptor
# Find columns with "Mouse" in the name
mouse_columns = ["<span title=\"Ligand Receptor Pair\">LR Pair</span>"] + [col for col in gene_pair.columns if "Mouse" in col]

# Filter rows where all "Mouse" columns are not " "
mouse_gene_pair = gene_pair[(gene_pair[mouse_columns].map(str.strip) != "").all(axis=1)]
# Reorder the DataFrame
new_order = mouse_columns + [col for col in mouse_gene_pair.columns if col not in mouse_columns]
mouse_gene_pair = mouse_gene_pair[new_order]
mouse_gene_pair = mouse_gene_pair.reset_index(drop=True)  

# Find columns with "Rat" in the name
rat_columns = ["<span title=\"Ligand Receptor Pair\">LR Pair</span>"] + [col for col in gene_pair.columns if "Rat" in col]
# Filter rows where any "Mouse" columns have a value
rat_gene_pair = gene_pair[(gene_pair[rat_columns].map(str.strip) != "").all(axis=1)]
# Reorder the DataFrame
new_order = rat_columns + [col for col in rat_gene_pair.columns if col not in rat_columns]
rat_gene_pair = rat_gene_pair[new_order]
rat_gene_pair = rat_gene_pair.reset_index(drop=True)  
