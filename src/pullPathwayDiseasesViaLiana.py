## Function to acquire pathway and disease annotations from Liana+

import liana as li
import omnipath as op
import decoupler as dc
import pandas as pd

import sys
import os
sys.path.append(os.path.abspath("src"))  # Add src directory to path
from createDataTable import gene_pair0

### Parameters
topN = 10000 #Number of top pathways to be included
pathway_output_file="data/pathway_annotations_per_pair.csv"


### Pathway Annotations

# load PROGENy pathways, we use decoupler as a proxy as it formats the data in a more convenient way
progeny = dc.get_progeny(top=topN)
# import connectomeDB database ligands and receptors
lr_pairs = gene_pair0[["Ligand", "Receptor"]]
lr_pairs.columns = lr_pairs.columns.str.lower()

# generate ligand-receptor geneset
lr_progeny = li.rs.generate_lr_geneset(lr_pairs, progeny, lr_sep="^")
# some of the pairs are missing
len(lr_progeny["interaction"].unique())

lr_progeny.to_csv(pathway_output_file, index=False)

### Disease Annotations

diseases = op.requests.Annotations.get(
    resources = ['DisGeNet']
    )

diseases = diseases[['genesymbol', 'label', 'value']]
diseases = diseases.pivot_table(index='genesymbol',
                                columns='label', values='value',
                                aggfunc=lambda x: '; '.join(x)).reset_index()
diseases = diseases[['genesymbol', 'disease']]
diseases['disease'] = diseases['disease'].str.split('; ')
diseases = diseases.explode('disease')
lr_diseases = li.rs.generate_lr_geneset(lr_pairs, diseases, source='disease', target='genesymbol', weight=None, lr_sep="^")
lr_diseases.sort_values("interaction")

# some of the pairs are missing
len(lr_diseases["interaction"].unique())

output_file="data/disease_annotations_per_pair.csv"
lr_diseases.to_csv(output_file, index=False)