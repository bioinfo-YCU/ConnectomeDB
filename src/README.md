### Quick guide of what each script does

**The following scripts are ALWAYS executed when qmds are being rendered via quarto:**

`fetchGSheet.py`: Fetch original data table created by Prof. Forrest from Google Sheet

`createDataTable.py`: Prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the database qmds

`createTriplicateDT.py`: Prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the database qmds

`createDataTable_perSpecies.py`: Prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the database qmds

`createFunctionalAnnotTable.py`: Prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the annotated database qmds

***

**The following scripts should be rerun in this particular order when there are new PMIDs:**

`webScrapePubMed.py`: Scrape data from Pubmed for Title, Abstract, Journal, and Year and saves results in data directory

`extractJourAbbv.py`: Make sure you get the abbreviation for all the Journal names

`prepareLLMinput.py`: to prepare input for `runLLMforBiologyRelevance.py`

`runLLMforBiologyRelevance.py`: to run LLM on all abstract for each LR pair to come up with biological relevance (keywords)

**AND**

**The following scripts should be rerun in this particular order when there are new information (disease, pathway, location), pairs, PMIDs:**

`createPMIDpages.py`: Create an evidence page per LR pair with each tab per PMID based on `pmidTemplate.html` from HTML directory

`pullPathwayDiseasesViaLiana.py`: Acquire pathway and disease annotations from Liana+

`plotExprAcrossCellType.py`: Create Expression plots per gene to be placed inside the PDF/HTML card

`createCards_scRNA.py`: Create Ligand-Receptor pair cards based on `cardTemplate_scRNA.html` from HTML directory with scRNA files

`createCardsWithPMID.py`: Create Ligand-Receptor pair cards based on `mergedCard_tabs.html` from HTML directory

<!-- `createCards.py`: Create Ligand-Receptor pair cards based on `cardTemplate.html` from HTML directory -->

**The following scripts should be rerun when there are new LR pairs:**

### for HGNC

`downloadHGNCfullList.r`
`downloadBioMartOrth.r`
`convertHGNCtoENSP.r`

### for an accurate number of pairs in the title
`rename_perSpecies.py`



