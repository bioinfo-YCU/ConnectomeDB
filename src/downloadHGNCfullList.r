# Purpose: Download HGNC complete txt file and save in data folder
# remember to use "R" conda env
# Libraries
# library(tidyverse)

# Main
print(paste0("Retrieving on: ", format(Sys.Date(), "%d-%b-%Y")))

# Download HGNC complete txt file
download.file(
  url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
  destfile = "hgnc_complete_set.txt",
  method = "auto"
)

# read into R
hgnc_data <- read.delim("hgnc_complete_set.txt", sep = "\t", header = TRUE)
# head(hgnc_data)
# Save to file
write.table(hgnc_data, file = "data/HGNC_gene_info_full.tsv", sep = "\t", row.names = FALSE, quote = FALSE)