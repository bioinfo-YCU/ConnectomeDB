import pandas as pd
import os, sys
import json
import numpy as np

# Add the src directory to the path for importing modules
sys.path.append(os.path.abspath("src"))
from createPMIDpages import gene_pair00

# Load the files
file1 = pd.read_csv("data/pubmed_results.csv") 
file2 = gene_pair00

# Convert PMIDs to string
file1['PMID'] = file1['PMID'].astype(str)

# Split PMIDs in gene_pair00
file2['PMID_List'] = file2['PMID support'].apply(lambda x: x.split(',') if isinstance(x, str) else [])

# Create PMID → Abstract map
pmid_to_abstract = dict(zip(file1['PMID'], file1['Abstract']))

# Get list of abstracts for each LR pair
def get_abstracts(pmids):
    return [pmid_to_abstract[pmid] for pmid in pmids if pmid in pmid_to_abstract]

file2['Abstracts'] = file2['PMID_List'].apply(get_abstracts)

# Clean up abstracts: make sure all entries are string lists
file2['Abstracts'] = file2['Abstracts'].apply(lambda x: [str(a) for a in x if isinstance(a, str)])
file2 = file2.replace({np.nan: None})

# Extract clean data
data_for_llm = file2[['Human LR Pair', 'Abstracts']].to_dict(orient='records')

# Save to JSON safely
with open("data/data_for_llm.json", "w", encoding="utf-8") as f:
    json.dump(data_for_llm, f, indent=4, ensure_ascii=False)
