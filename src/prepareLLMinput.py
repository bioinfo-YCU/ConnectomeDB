## Function to prepare input for LLM where all abstract for each LR pair are put together in one place
import pandas as pd
import os, sys
import json
# Add the src directory to the path for importing modules
sys.path.append(os.path.abspath("src"))
from createPMIDpages import gene_pair00

# Load the files
file1 = pd.read_csv("data/pubmed_results.csv") 
file2 = gene_pair00

# Convert the PMIDs column in file2 to lists for easy comparison
file1['PMID'] = file1['PMID'].astype(str)
file2['PMID_List'] = file2['PMID support'].apply(lambda x: x.split(','))

# Create a dictionary for quick PMID to Abstract mapping
pmid_to_abstract = dict(zip(file1['PMID'], file1['Abstract']))

# Function to get all abstracts for a list of PMIDs
def get_abstracts(pmids):
    return [pmid_to_abstract[pmid] for pmid in pmids if pmid in pmid_to_abstract]

# Map abstracts to LR pairs
file2['Abstracts'] = file2['PMID_List'].apply(get_abstracts)

# Convert to a list of dictionaries
data_for_llm = file2[['Human LR Pair', 'Abstracts']].to_dict(orient='records')

# Save as JSON
with open("data/data_for_llm.json", "w") as f:
    json.dump(data_for_llm, f, indent=4)