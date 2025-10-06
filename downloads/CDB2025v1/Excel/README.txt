Version: CDB2025v1
Date released: Sep 25 2025

Column Definitions

Column	Name	Description
0	Interaction ID	Unique ConnectomeDB identifier for each ligand–receptor pair.
1	LR Pair	Human-readable label of the ligand–receptor pair (e.g., TGFB1–TGFBR1).
2	Evidence	Direct: experimentally verified; Inferred: inferred from orthology supporting the interaction.
3	AI summary	AI-generated summary describing the interaction; can be expanded using external tools.
4	Ligand Symbols	Official gene symbol(s) for the ligand, including known aliases or old names.
5	Receptor Symbols	Official gene symbol(s) for the receptor, including known aliases or old names.
6	Ligand ID	Species-specific gene/protein ID for the ligand. See species mapping below.
7	Receptor ID	Species-specific gene/protein ID for the receptor. See species mapping below.
8	Ligand ENSEMBL ID	ENSEMBL gene ID for the ligand.
9	Receptor ENSEMBL ID	ENSEMBL gene ID for the receptor.
10	Human Ligand Symbols	Human ortholog gene symbol(s) for the ligand.
11	Human Receptor Symbols	Human ortholog gene symbol(s) for the receptor.
12	Ligand Location	Predicted subcellular localization of the human ligand protein.
13	Receptor Location	Predicted subcellular localization of the human receptor protein.

Species-Specific Mapping for Columns 6 and 7

The contents of columns 6 (Ligand ID) and 7 (Receptor ID) depend on the species of the dataset:

Species	Ligand ID Column Name	Receptor ID Column Name	ID Type / Database Used
Mouse	Ligand MGI ID	Receptor MGI ID	MGI: Mouse Genome Informatics
Rat	Ligand RGD ID	Receptor RGD ID	RGD: Rat Genome Database
Zebrafish	Ligand ZFIN ID	Receptor ZFIN ID	ZFIN: Zebrafish Information Network
Frog	Ligand XEN ID	Receptor XEN ID	XEN: Xenbase (frog model)
Human	Ligand HGNC ID	Receptor HGNC ID	HGNC: Human Gene Nomenclature
Other	Ligand XX ID	Receptor XX ID	Placeholder for future species/databases


If you prefer one table/file that includes all species, you can download the `all_species.xlsx`.
It contains 15 columns including a `Species` columnthat identifies the species for each ligand–receptor pair. 
The individual IDs— Ligand {HGNC/MGI/RGD/XEN/ZFIN/XX} ID —are combined into the Ligand Species ID column, 
and Receptor {HGNC/MGI/RGD/XEN/ZFIN/XX} ID are combined into the Receptor Species ID column. 
   
For bulk downloads, you can download a zipped folder, `all_species_sep_files.zip`, 
which contains the data tables for all 14 species in XLSX.
Note that this does not include the `all_species.xlsx` file.