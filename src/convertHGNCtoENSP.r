# Purpose: Download HGNC complete txt file and save in data folder
# remember to use "R" conda env
# Libraries
library(biomaRt)

# Main
print(paste0("Retrieving on: ", format(Sys.Date(), "%d-%b-%Y")))

# Connect to Ensembl BioMart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human_genes <- getBM(
        attributes = c("hgnc_id", "ensembl_peptide_id", "transcript_mane_select"),
        mart = ensembl
    )
human_genes <- unique(human_genes)
human_genes <- human_genes[!human_genes$hgnc_id=="",]
human_genes <- human_genes[!human_genes$ensembl_peptide_id =="", ]
human_genes <- human_genes[!human_genes$transcript_mane_select =="", ]
human_genes <- human_genes[!is.na(human_genes$ensembl_peptide_id), ]
# Save to file
readr::write_csv(human_genes, "data/hgnc_ensp_biomart.csv")