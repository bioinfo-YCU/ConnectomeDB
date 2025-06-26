# Purpose: Download Ensembl to UniProt mapping for any species
# Usage: Change `species_name` below to one of the supported options

library(biomaRt)
library(tidyverse)

# Set species here
species_name <- "horse"  # options: mouse, rat, cow, etc.

# Lookup table: species → Ensembl code
species_lookup <- list(
  mouse       = "mmusculus_gene_ensembl",
  rat         = "rnorvegicus_gene_ensembl",
  zebrafish   = "drerio_gene_ensembl",
  chimpanzee  = "ptroglodytes_gene_ensembl",
  chicken     = "ggallus_gene_ensembl",
  pig         = "sscrofa_gene_ensembl",
  cow         = "btaurus_gene_ensembl",
  dog         = "clfamiliaris_gene_ensembl",
  horse       = "ecaballus_gene_ensembl",
  sheep       = "oarambouillet_gene_ensembl",
  marmoset    = "cjacchus_gene_ensembl",
  macaque     = "mmulatta_gene_ensembl"
)

# Resolve dataset name
if (!(species_name %in% names(species_lookup))) {
  stop(paste("Species not supported:", species_name))
}

dataset_name <- species_lookup[[species_name]]

# Log retrieval time
cat("Retrieving for", species_name, "on", format(Sys.Date(), "%d-%b-%Y"), "\n")

# Connect to BioMart
ensembl <- useMart("ensembl", dataset = dataset_name)

# Retrieve mapping
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "uniprotswissprot"),
  mart = ensembl
)

# Filter and deduplicate
gene_map <- gene_map %>%
  filter(uniprotswissprot != "", !is.na(uniprotswissprot)) %>%
  distinct()

# Save to file
out_file <- paste0("data/", species_name, "_uniprot_to_ensembl.tsv")
write.table(gene_map, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Saved to", out_file, "\n")
