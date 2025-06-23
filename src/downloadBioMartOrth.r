# Purpose: Download HGNC complete txt file and save in data folder
# remember to use "R" conda env
# Libraries
library(biomaRt)
library(tidyverse)

print(paste0("Retrieving on: ", format(Sys.Date(), "%d-%b-%Y")))

# Function to retrieve orthologs for different species
library(biomaRt)
library(readr)

get_species_orthologs <- function(species_name) {
  # Connect to Ensembl BioMart (human)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # Retrieve all human genes with their Ensembl IDs and HGNC symbols
  human_genes <- getBM(
        attributes = c("hgnc_id", "hgnc_symbol", "ensembl_gene_id"),
        mart = ensembl
    )
  # Define species-specific ortholog attribute names
  species_column <- paste0(species_name, "_homolog_ensembl_gene")
  species_gene_name <- paste0(species_name, "_homolog_associated_gene_name")
  species_gene_GOC <- paste0(species_name, "_homolog_goc_score")
  species_gene_WGA <- paste0(species_name, "_homolog_wga_coverage")
  species_gene_identToQuery <- paste0(species_name, "_homolog_perc_id") 
  species_gene_identToTarget <- paste0(species_name, "_homolog_perc_id_r1") 
  species_gene_confidence <- paste0(species_name, "_homolog_orthology_confidence") 

  # Ensure attributes are from a single species' ortholog page
  attributes_list <- c(
    "ensembl_gene_id",
    species_column, 
    species_gene_name,
    species_gene_GOC,
    species_gene_WGA,
    species_gene_identToQuery,
    species_gene_identToTarget,
    species_gene_confidence
  )
  
  # Retrieve orthologs from human dataset
  orth_genes <- getBM(
    attributes = attributes_list,
    mart = ensembl
  )
  result <- merge(human_genes, orth_genes, by = "ensembl_gene_id", all.x = TRUE)
  # Connect to ortholog species dataset
  ensembl_orthologs <- useMart("ensembl", dataset = paste0(species_name, "_gene_ensembl"))
  
  # Set species-specific gene ID and description fields
  species_id <- switch(species_name,
                       "mmusculus" = "mgi_id",
                       "rnorvegicus" = "rgd_id",
                       "drerio" = "zfin_id_id",
                       "ensembl_gene_id")  # default

  species_gene_name_long <- if (species_name == "mmusculus") {
    "mgi_description"
  } else {
    "description"
  }

  # Get gene info from ortholog species dataset
  orthologs <- getBM(
    attributes = c("ensembl_gene_id", species_id, species_gene_name_long),
    mart = ensembl_orthologs
  )

  # Merge based on the ortholog ensembl gene ID (species_column)
  colnames(orthologs)[1] <- species_column  # rename for correct merge
  final_result <- merge(result, orthologs, by = species_column, all.x = TRUE)
  final_result <- final_result[final_result$hgnc_id != "", ]
  final_result <- final_result[!is.na(final_result$hgnc_id), ]
  # Write result
  write_csv(final_result, paste0("data/", species_name, "_ID_biomart.csv"))

  return(final_result)
}

# searchAttributes(ensembl, pattern) useful for looking for column names

# Main

#  Rhesus macaque (Macaca mulatta)
get_species_orthologs("mmulatta")
#  Marmoset (Callithrix jacchus)
get_species_orthologs("cjacchus")
# Chimp (Pan troglodytes)
get_species_orthologs("ptroglodytes")
# Chicken (Gallus gallus)
get_species_orthologs("ggallus")
# Pig (Sus scrofa)
get_species_orthologs("sscrofa")
# Cow (Bos taurus)
get_species_orthologs("btaurus")
# Dog (Canis lupus familiaris)
get_species_orthologs("clfamiliaris")
# Horse (Equus caballus)
get_species_orthologs("ecaballus")
# Sheep (Ovis aries rambouillet)
get_species_orthologs("oarambouillet")
# Mouse (Rattus norvegicus) # test
get_species_orthologs("rnorvegicus")
# Rat (Mus musculus) # test
get_species_orthologs("mmusculus")
# Zebrafish (Danio rerio) # test
get_species_orthologs("drerio")