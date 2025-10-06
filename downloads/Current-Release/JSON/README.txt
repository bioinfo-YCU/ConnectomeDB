Version: Current-Release
Date released: Sep 25 2025

Column Definitions

Column	Name	Description
0	Interaction_ID	Unique ConnectomeDB identifier for each ligand–receptor pair.
1	LR_Pair	Human-readable label of the ligand–receptor pair (e.g., TGFB1–TGFBR1).
2	Evidence	Direct: experimentally verified; Inferred: inferred from orthology supporting the interaction.
3	AI_summary	AI-generated summary describing the interaction; can be expanded using external tools.
4	Ligand_Symbols	Official gene symbol(s) for the ligand, including known aliases or old names.
5	Receptor_Symbols	Official gene symbol(s) for the receptor, including known aliases or old names.
6	Ligand_ID	Species-specific gene/protein ID for the ligand. See species mapping below.
7	Receptor_ID	Species-specific gene/protein ID for the receptor. See species mapping below.
8	Ligand_ENSEMBL_ID	ENSEMBL gene ID for the ligand.
9	Receptor_ENSEMBL_ID	ENSEMBL gene ID for the receptor.
10	Human_Ligand_Symbols	Human ortholog gene symbol(s) for the ligand.
11	Human_Receptor_Symbols	Human ortholog gene symbol(s) for the receptor.
12	Ligand_Location	Predicted subcellular localization of the human ligand protein.
13	Receptor_Location	Predicted subcellular localization of the human receptor protein.

Species-Specific Mapping for Columns 6 and 7

The contents of columns 6 (Ligand_ID) and 7 (Receptor_ID) depend on the species of the dataset:

Species	Ligand ID Column Name	Receptor ID Column Name	ID Type / Database Used
Mouse	Ligand_MGI_ID	Receptor_MGI_ID	MGI: Mouse Genome Informatics
Rat	Ligand_RGD_ID	Receptor_RGD_ID	RGD: Rat Genome Database
Zebrafish	Ligand_ZFIN_ID	Receptor_ZFIN_ID	ZFIN: Zebrafish Information Network
Frog	Ligand_XEN_ID	Receptor_XEN_ID	XEN: Xenbase (frog model)
Human	Ligand_HGNC_ID	Receptor_HGNC_ID	HGNC: Human Gene Nomenclature
Other	Ligand_XX_ID	Receptor_XX_ID	Placeholder for future species/databases


If you prefer one table/file that includes all species, you can download the `all_species.json`.
It contains 15 columns including a `Species` column that identifies the species for each ligand–receptor pair.
The individual IDs— Ligand {HGNC/MGI/RGD/XEN/ZFIN/XX} ID —are combined into the Ligand Species ID column, 
and Receptor {HGNC/MGI/RGD/XEN/ZFIN/XX} ID are combined into the Receptor Species ID column. 
   
For bulk downloads, you can download a zipped folder, `all_species_sep_files.zip`, 
which contains the data tables for all 14 species in JSON.
Note that this does not include the `all_species.json` file.