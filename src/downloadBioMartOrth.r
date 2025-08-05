# Purpose: Download Biomart Orth info from original species of interest
# remember to use "R" conda env
# Libraries
library(biomaRt)
library(tidyverse)
library(readr)
print(paste0("Retrieving on: ", format(Sys.Date(), "%d-%b-%Y")))

get_species_orthologs <- function(orig_species, species_name) {
  # Set species-specific gene ID and symbol fields for original species
  orig_species_id <- switch(orig_species,
                            "mmusculus" = "mgi_id",
                            "rnorvegicus" = "rgd_id",
                            "drerio" = "zfin_id_id",
                            "hsapiens" = "hgnc_id",
                            "ensembl_gene_id")  # default
  
  orig_species_symbol <- switch(orig_species,
                                "mmusculus" = "mgi_symbol",
                                "rnorvegicus" = "rgd_symbol",
                                "drerio" = "external_gene_name",
                                "hsapiens" = "hgnc_symbol",
                                "external_gene_name")  # default
  
  # Set gene ID field for target species
  species_id <- switch(species_name,
                       "mmusculus" = "mgi_id",
                       "rnorvegicus" = "rgd_id",
                       "drerio" = "zfin_id_id",
                       "hsapiens" = "hgnc_id",
                       "ensembl_gene_id")  # default
  
  # Connect to Ensembl BioMart (original species)
  ensembl <- biomaRt::useMart("ensembl", dataset = paste0(orig_species, "_gene_ensembl"))
  
  # Retrieve genes from original species
  orig_attributes <- c(
      orig_species_id, 
      orig_species_symbol, 
      "external_synonym",  # common for many species
      "ensembl_gene_id"
    )
    
  orig_genes <- biomaRt::getBM(
    attributes = c(orig_attributes),
    mart = ensembl
  )
  
  # Define ortholog attribute names
  species_column <- paste0(species_name, "_homolog_ensembl_gene")
  species_gene_name <- paste0(species_name, "_homolog_associated_gene_name")
  species_gene_GOC <- paste0(species_name, "_homolog_goc_score")
  species_gene_WGA <- paste0(species_name, "_homolog_wga_coverage")
  species_gene_identToQuery <- paste0(species_name, "_homolog_perc_id")
  species_gene_identToTarget <- paste0(species_name, "_homolog_perc_id_r1")
  species_gene_confidence <- paste0(species_name, "_homolog_orthology_confidence")
  
  # Get orthologs from original species mart
    available_attributes <- biomaRt::listAttributes(ensembl)[["name"]]
    
    wga_attribute <- paste0(species_name, "_homolog_wga_coverage")
    
    include_wga <- wga_attribute %in% available_attributes
    
    attributes_list <- c(
      "ensembl_gene_id",
      species_column,
      species_gene_name,
      species_gene_GOC,
      species_gene_identToQuery,
      species_gene_identToTarget,
      species_gene_confidence
    )
    
    if (include_wga) {
      attributes_list <- append(attributes_list, wga_attribute, after = 4)
    }
    
      orth_genes <- biomaRt::getBM(
        attributes = attributes_list,
        mart = ensembl
      )
      
  result <- merge(orig_genes, orth_genes, by = "ensembl_gene_id", all.x = TRUE)
  
  # Connect to ortholog species dataset
  ensembl_orthologs <- biomaRt::useMart("ensembl", dataset = paste0(species_name, "_gene_ensembl"))
  
  species_gene_name_long <- if (species_name == "mmusculus") {
    "mgi_description"
  } else {
    "description"
  }
  
  # Get additional info from ortholog species
  orthologs <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", species_id, species_gene_name_long),
    mart = ensembl_orthologs
  )
  
  # Merge to enrich with species-specific ID/symbol
  colnames(orthologs)[1] <- species_column  # match column for merge
  final_result <- merge(result, orthologs, by = species_column, all.x = TRUE)
  
  # Remove rows with empty or NA original species ID
  final_result <- final_result[final_result[[orig_species_id]] != "", ]
  final_result <- final_result[!is.na(final_result[[orig_species_id]]), ]
  
  # Write result
  readr::write_csv(final_result, paste0("data/", species_name, "_ID_biomart_", orig_species, "_centric.csv"))
  
  return(final_result)
}

#  mouse to human
get_species_orthologs("mmusculus", "hsapiens")

#  mouse to rat
get_species_orthologs("mmusculus", "btaurus")

#  mouse to cow
get_species_orthologs("mmusculus", "btaurus")

#  Rhesus macaque (Macaca mulatta) to human
get_species_orthologs("mmulatta", "hsapiens")

#  Marmoset (Callithrix jacchus) to human
get_species_orthologs("cjacchus", "hsapiens")

# Chimp (Pan troglodytes) to human
get_species_orthologs("ptroglodytes", "hsapiens")

# Chicken (Gallus gallus) to human
get_species_orthologs("ggallus", "hsapiens")

# Pig (Sus scrofa) to human
get_species_orthologs("sscrofa", "hsapiens")

# Cow (Bos taurus) to human
get_species_orthologs("btaurus", "hsapiens")

# Dog (Canis lupus familiaris) to human
get_species_orthologs("clfamiliaris", "hsapiens")

# Horse (Equus caballus) to human
get_species_orthologs("ecaballus", "hsapiens")

# Sheep (Ovis aries rambouillet) to human
get_species_orthologs("oarambouillet", "hsapiens")

# Rat (Rattus norvegicus) to human
get_species_orthologs("rnorvegicus", "hsapiens")

# Rat (Rattus norvegicus) to mouse
get_species_orthologs("rnorvegicus", "mmusculus")

# Zebrafish (Danio rerio) to human
get_species_orthologs("drerio", "hsapiens")

# Guinea pig to Human
get_species_orthologs("cporcellus", "hsapiens")

# Rabbit to Human
get_species_orthologs("ocuniculus", "hsapiens")

# green-spotted puffer to human
get_species_orthologs("tnigroviridis", "hsapiens")

# Frog has no "hsapiens_homolog_wga_coverage"
get_species_orthologs("xtropicalis", "hsapiens")

# searchAttributes(ensembl, pattern) useful for looking for column names

# Main

#  human to Rhesus macaque (Macaca mulatta)
get_species_orthologs("hsapiens", "mmulatta")
# human to Marmoset (Callithrix jacchus)
get_species_orthologs("hsapiens", "cjacchus")
# human to Chimp (Pan troglodytes)
get_species_orthologs("hsapiens", "ptroglodytes")
# human to Chicken (Gallus gallus)
get_species_orthologs("hsapiens", "ggallus")
# human to Pig (Sus scrofa)
get_species_orthologs("hsapiens", "sscrofa")
# human to Cow (Bos taurus)
get_species_orthologs("hsapiens", "btaurus")
# human to Dog (Canis lupus familiaris)
get_species_orthologs("hsapiens", "clfamiliaris")
# human to Horse (Equus caballus)
get_species_orthologs("hsapiens", "ecaballus")
# human to Sheep (Ovis aries rambouillet)
get_species_orthologs("hsapiens", "oarambouillet")
# human to Rat (Rattus norvegicus)  
get_species_orthologs("hsapiens", "rnorvegicus")
# human to Zebrafish (Danio rerio) 
get_species_orthologs("hsapiens", "drerio")
# human to Mouse (Danio rerio) 
get_species_orthologs("hsapiens", "mmusculus")
# Human to Rabbit
get_species_orthologs("hsapiens", "ocuniculus")
# Human to Guinea pig
get_species_orthologs("hsapiens", "cporcellus")
# Human to green-spotted puffer
get_species_orthologs("hsapiens", "tnigroviridis")
# Human has no frog "xtropicalis_homolog_wga_coverage"
get_species_orthologs("hsapiens", "xtropicalis")

