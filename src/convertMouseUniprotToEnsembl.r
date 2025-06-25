# Purpose: Download HGNC complete txt file and save in data folder
# remember to use "R" conda env
# Libraries
library(biomaRt)
library(tidyverse)

print(paste0("Retrieving on: ", format(Sys.Date(), "%d-%b-%Y")))

# Function to retrieve orthologs for different species
library(biomaRt)
library(readr)

# Connect to Ensembl BioMart
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Retrieve mapping of Ensembl Gene ID to UniProtKB/Swiss-Prot
mouse_map <- getBM(
    attributes = c("ensembl_gene_id", "uniprotswissprot"),
    mart = ensembl
)

# Optional: filter out empty UniProt entries
mouse_map <- subset(mouse_map, uniprotswissprot != "" & !is.na(uniprotswissprot))

# Optional: remove duplicates
mouse_map <- unique(mouse_map)

# Save to file
write.table(mouse_map, file = "data/mouse_uniprot_to_ensembl.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
