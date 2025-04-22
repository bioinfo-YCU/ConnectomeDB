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

import numpy as np

# Clean data
file2['Abstracts'] = file2['Abstracts'].apply(lambda x: [str(a) for a in x if isinstance(a, str)])
file2 = file2.replace({np.nan: None})

# Convert to list of dicts
data_for_llm = file2[['Human LR Pair', 'Abstracts']].to_dict(orient='records')

# Dump safely
with open("data/data_for_llm.json", "w", encoding="utf-8") as f:
    json.dump(data_for_llm, f, indent=4, ensure_ascii=False)
