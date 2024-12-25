### Quick guide of what each script does

The following scripts are always executed when qmds are being rendered via quarto:

`fetchGSheet.py`: Fetch original data table created by Prof. Forrest from Google Sheet

`createDataTable.py`: Prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the database qmds

***

The following scripts should be remade when there are new pairs or PMIDs (once a day should be enough):

`createCards.py`: Create Ligand-Receptor pair cards based on `cardTemplate.html` from HTML directory

`webScrapePubMed.py`: Scrape data from Pubmed for Title, Abstract, Journal, and Year and saves results in data directory

`createPMIDpages.py`: Create an evidence page per LR pair with each tab per PMID based on `pmidTemplate.html` from HTML directory
