import requests
import pandas as pd
import sys
import os

sys.path.append(os.path.abspath("src"))  # Add src directory to path
from createDataTable import gene_pair0

# Get all unique genes
ligand_list = gene_pair0["Ligand"].tolist()
receptor_list = gene_pair0["Receptor"].tolist()
unique_genes = list(set(ligand_list + receptor_list))  # Combine and remove duplicates
import re
import pandas as pd
import gzip

# Input FASTA file
fasta_file = "data/uniprotkb_proteome_UP000005640_AND_revi_2025_03_27.fasta.gz"

# Regex pattern to extract details from header
header_pattern = re.compile(r"^>sp\|(?P<uniprot_id>[A-Z0-9]+(?:-\d+)?)\|(?P<protein_name>.+?) OS=Homo sapiens OX=9606 GN=(?P<gene_name>[A-Za-z0-9-]+)")

# Store extracted data
records = []

# Read and parse the file
with gzip.open(fasta_file, "rt") as f:
    header = None
    sequence = []
    
    for line in f:
        line = line.strip()
        
        if line.startswith(">"):
            # Store previous sequence if exists
            if header and sequence:
                isoform_type = "Canonical" if "-" not in header["uniprot_id"] else "Alternative Isoform"
                records.append([header["uniprot_id"], header["gene_name"], header["protein_name"], isoform_type, "".join(sequence)])
            
            # Match new header
            match = header_pattern.match(line)
            if match:
                header = match.groupdict()
                sequence = []
            else:
                header = None
        
        elif header:
            sequence.append(line)

    # Add the last record
    if header and sequence:
        isoform_type = "Canonical" if "-" not in header["uniprot_id"] else "Alternative Isoform"
        records.append([header["uniprot_id"], header["gene_name"], header["protein_name"], isoform_type, "".join(sequence)])

# Convert to pandas DataFrame
df = pd.DataFrame(records, columns=["UniProt ID", "Gene Symbol", "Protein Name", "Isoform Type", "FASTA Sequence"])

# Save as TSV
df.to_csv("data/human_uniprot_isoforms.tsv", sep="\t", index=False)

print(f"✅ Extracted {len(df)} Homo sapiens protein sequences and saved to 'human_uniprot_isoforms.tsv'.")

### Merging with LR pairs
df = df[['UniProt ID', 'Gene Symbol', 'Isoform Type', 'FASTA Sequence']]
lim_df = gene_pair0[["Human LR Pair", "Ligand", "Receptor"]]
lim_df = lim_df.merge(df, how='left', left_on='Ligand', right_on='Gene Symbol')
lim_df = lim_df.drop(columns=["Gene Symbol"])
lim_df = lim_df.rename(columns={"FASTA Sequence": "Ligand Sequence",
                                "UniProt ID": "Ligand Isoform Uniprot ID",
                                "Isoform Type": "Ligand Isoform Type"})
lim_df

lim_df = lim_df.merge(df, how='left', left_on='Receptor', right_on='Gene Symbol')
lim_df = lim_df.drop(columns=["Gene Symbol"])
lim_df = lim_df.rename(columns={"FASTA Sequence": "Receptor Sequence",
                                "UniProt ID": "Receptor Isoform Uniprot ID",
                                "Isoform Type": "Receptor Isoform Type"})
lim_df.to_csv("data/LRpair_uniprot_sequences.tsv", sep="\t", index=False)