import requests
import pandas as pd
import gzip
import re
import sys
import os

sys.path.append(os.path.abspath("src"))  # Add src directory to path
from createDataTable import gene_pair0

# Step 1: Extract Gene Symbol Mapping from GTF
gtf_file = "data/gencode.v47.annotation.gtf.gz"
gene_map = {}

# Read GTF file and extract gene_id -> gene_name mapping
with gzip.open(gtf_file, "rt") as f:
    for line in f:
        if line.startswith("#"):  # Skip comments
            continue
        
        fields = line.strip().split("\t")
        if fields[2] == "gene":  # Only extract gene entries
            info = {key.strip(): value.strip('"') for key, value in re.findall(r'(\S+) "([^"]+)"', fields[8])}
            if "gene_id" in info and "gene_name" in info:
                gene_map[info["gene_id"]] = info["gene_name"]

print(f"✅ Extracted {len(gene_map)} gene mappings from GTF.")

# Step 2: Parse GENCODE Protein FASTA and Add Gene Symbols
fasta_file = "data/gencode.v47.pc_translations.fa.gz"

# Store extracted data
records = []

# Open the GENCODE FASTA file and parse sequences
with gzip.open(fasta_file, "rt") as f:
    header = None
    sequence = []
    
    for line in f:
        line = line.strip()
        
        if line.startswith(">"):
            # Store previous sequence if exists
            if header and sequence:
                # Extract the Gene Symbol using the Gene ID
                gene_symbol = gene_map.get(header["gene_id"], "Unknown")
                isoform_type = "Canonical" if "-1" in header["protein_id"] else "Alternative Isoform"
                
                # Append the parsed data to records
                records.append([header["protein_id"], header["transcript_id"], header["gene_id"], gene_symbol, isoform_type, "".join(sequence)])
            
            # Split header by '|' and extract necessary fields
            fields = line[1:].split("|")  # Skip the '>' symbol and split by '|'
            if len(fields) >= 6:
                header = {
                    "protein_id": fields[0], 
                    "transcript_id": fields[1], 
                    "gene_id": fields[2]  
                }
                sequence = []
            else:
                header = None
        
        elif header:
            sequence.append(line)

    # Add the last record if needed
    if header and sequence:
        gene_symbol = gene_map.get(header["gene_id"], "Unknown")
        isoform_type = "Canonical" if "-1" in header["protein_id"] else "Alternative Isoform"
        records.append([header["protein_id"], header["transcript_id"], header["gene_id"], gene_symbol, isoform_type, "".join(sequence)])

# Step 3: Convert to pandas DataFrame and Save to TSV
df = pd.DataFrame(records, columns=["Ensembl Protein ID", "Ensembl Transcript ID", "Ensembl Gene ID", "Gene Symbol", "Isoform Type", "FASTA Sequence"])

# Save to TSV
df.to_csv("data/human_gencode_isoforms.tsv", sep="\t", index=False)

# Print completion message
print(f"✅ Extracted {len(df)} protein sequences with Gene Symbols and saved to 'gencode_protein_isoforms_with_symbols.tsv'.")

### Merge with LR pairs

lim_df = gene_pair0[["Human LR Pair", "Ligand", "Receptor"]]
lim_df = lim_df.merge(df, how='left', left_on='Ligand', right_on='Gene Symbol')
lim_df = lim_df.drop(columns=["Gene Symbol"])
lim_df = lim_df.rename(columns={"FASTA Sequence": "Ligand Sequence",
                                "Ensembl Protein ID": "Ligand Ensembl Protein ID",
                                "Ensembl Transcript ID": "Ligand Ensembl Transcript ID",
                                "Ensembl Gene ID": "Ligand Ensembl Gene ID",
                                "Isoform Type": "Ligand Isoform Type"})

lim_df = lim_df.merge(df, how='left', left_on='Receptor', right_on='Gene Symbol')
lim_df = lim_df.drop(columns=["Gene Symbol"])
lim_df = lim_df.rename(columns={"FASTA Sequence": "Receptor Sequence",
                                "Ensembl Protein ID": "Receptor Ensembl Protein ID",
                                "Ensembl Transcript ID": "Receptor Ensembl Transcript ID",
                                "Ensembl Gene ID": "Receptor Ensembl Gene ID",
                                "Isoform Type": "Receptor Isoform Type"})

lim_df.to_csv("data/LRpair_gencode_sequences.tsv", sep="\t", index=False)

