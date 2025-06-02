# Purpose: Download HGNC complete txt file and save in data folder
# remember to use "R" conda env
# Libraries
library(biomaRt)
library(tidyverse)

print(paste0("Retrieving on: ", format(Sys.Date(), "%d-%b-%Y")))

# Function to retrieve orthologs for different species
get_species_orthologs <- function(species_name) {
  # Connect to Ensembl BioMart
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  # Retrieve all human genes with their Ensembl IDs and HGNC symbols
    human_genes <- getBM(
        attributes = c("hgnc_id", "hgnc_symbol", "ensembl_gene_id"),
        mart = ensembl
    )
  
  # Define species-specific attributes dynamically based on species name
  species_column <- paste0(species_name, "_homolog_ensembl_gene")
  species_gene_name <- paste0(species_name, "_homolog_associated_gene_name")
  
  # Get orthologs for the specified species
  orthologs <- getBM(
    attributes = c("ensembl_gene_id", 
                   species_column, 
                   species_gene_name),
    mart = ensembl
  )
  final_result <- merge(human_genes, orthologs, by = "ensembl_gene_id", all.x = TRUE)
  readr::write_csv(final_result, paste0("data/",species_name, "_ID_biomart.csv"))
  return(final_result)
}

# Function to list available attributes for a given Ensembl species dataset
list_species_attributes <- function(species_dataset) {
  ensembl <- useMart("ensembl", dataset = species_dataset)
  attributes <- listAttributes(ensembl)
  return(attributes)
}

# Main

#  Rhesus macaque (Macaca mulatta)
get_species_orthologs("mmulatta")
#  Marmoset (Callithrix jacchus)
get_species_orthologs("cjacchus")
# Chimp (Pan troglodytes)
get_species_orthologs("ptroglodytes")
# Chicken (Gallus gallus)
chicken_orthologs <- get_species_orthologs("ggallus")
# Pig (Sus scrofa)
pig_orthologs <- get_species_orthologs("sscrofa")
# Cow (Bos taurus)
cow_orthologs <- get_species_orthologs("btaurus")
# Dog (Canis lupus familiaris)
dog_orthologs <- get_species_orthologs("clfamiliaris")
# Horse (Equus caballus)
horse_orthologs <- get_species_orthologs("ecaballus")
# Sheep (Ovis aries rambouillet)
sheep_orthologs <- get_species_orthologs("oarambouillet")
# Mouse (Mus musculus) # test
mouse_orthologs <- get_species_orthologs("mmusculus")
# Zebrafish (Danio rerio) # test
zebra_orthologs <- get_species_orthologs("drerio")