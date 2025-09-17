## Function to prepare datatables (cleaning and hyperlinking, adding tool tips, etc) input for the database qmds
import sys, os
from itables import init_notebook_mode
import pandas as pd
from itables import show
from itables import options
from IPython.display import HTML, display
import numpy as np
import re
from bs4 import BeautifulSoup
from createDataTable import gene_pair, gene_pair000, human_columns, lrPairsCount, replace_link_text
import warnings
import fetchGSheet 

# Suppress SettingWithCopyWarning
warnings.simplefilter("ignore", category=UserWarning)

def process_species_gene_pair(species, fetchGSheet, gene_pair):
    species_name = {
        "Mouse": "mmusculus",
        "Rat": "rnorvegicus",
        "Zebrafish":"drerio" ,
        "Chimp":"ptroglodytes",
        "Chicken":"ggallus",
        "Pig":"sscrofa",
        "Cow":"btaurus",
        "Dog":"clfamiliaris",
        "Horse":"ecaballus",
        "Sheep":"oarambouillet",
        "Marmoset": "cjacchus" ,
        "Macaque": "mmulatta",
        "Frog": "xtropicalis"
    }.get(species, "Unknown species")

    if species == "Mouse":
        species_id = "MGI"
        species_info = pd.read_csv(f"data/MRK_Merged_{species_id}_DB.tsv", sep="\t", dtype=str)
    elif species == "Rat":
        species_id = "RGD"
        species_info = pd.read_csv(f"data/GENES_RAT_{species_id}_DB.tsv", sep="\t", dtype=str)
        # combine all known aliases also old ones
        species_info["ALIASES"] = species_info[["MARKER_SYMBOL", "OLD_SYMBOL"]].apply(
            lambda row: ";".join(pd.unique(row.dropna().astype(str))), axis=1
        )
    elif species == "Frog":
        species_id = "XEN"
        species_info = pd.read_csv(f"data/GenePageGeneralInfo_{species_id.capitalize()}base_DB.tsv", sep="\t", dtype=str)
    elif species == "Zebrafish":
        species_id = "ZFIN"
        species_info = pd.read_csv(f"data/Zebrafish_merged_{species_id}_DB.tsv", sep="\t", dtype=str)
        # combine all known aliases also old ones
        species_info["Aliases"] = species_info[["Current Name", "Previous Name"]].apply(
            lambda row: ", ".join(pd.unique(row.dropna().astype(str))), axis=1
        )
    else: 
        species_id = "ENSEMBL"
        species_info = pd.read_csv(f"data/hsapiens_ID_biomart_{species_name}_centric.csv", dtype=str)
    
    
    species_lower = species.lower()
    
    if species == "Mouse":
        gene_pair_species = getattr(fetchGSheet, f"gene_pair_{species_lower}")  
    else:
        gene_pair_species = fetchGSheet.safe_fetch(fetchGSheet.sheet_ID, f"FROZEN_{species_lower}", fetchGSheet.credentials_file)
    
    def extract_visible_text(col):
        """Extract visible text between '>' and '</a>'."""
        match = re.search(r'>([^<]+)</a>', col)
        if match:
            return match.group(1).strip()
        return None
    
    # Grab the Interaction ID and remove unnecessary columns
    keywords = ["Evidence</span>", "A.I. summary", " HGNC ID</span>", " ENSEMBL ID</span>"]

    # Keep only columns that do NOT contain any of the keywords
    gene_pair = gene_pair.loc[:, ~gene_pair.columns.str.contains('|'.join(keywords))]

    # Next, drop columns at index positions 3 and 4 ("Ligand and Receptor" since we already have ligand symbols and receptor symbols)
    gene_pair = gene_pair.drop(gene_pair.columns[[3, 4]], axis=1)
    exclude_keywords = ["HGNC ID", "Location", "Human"]  # Columns containing this will not be modified
    keywords_to_modify = ["Ligand Symbols", "Receptor Symbols"]
    # Copy the original columns so we can modify only the first 10
    new_columns = gene_pair.columns.tolist()
    
    # Modify only the first 10 columns
    new_columns = [
        f'{col.split(">")[0]}">Human {col.split(">")[1]}</span>'
        if any(keyword in col for keyword in keywords_to_modify) and not any(exclude in col for exclude in exclude_keywords)
        else col
        for col in new_columns
    ]
    # Assign the modified column names back to the DataFrame
    gene_pair.columns = new_columns
    
    gene_pair["LR Pair Card"] = gene_pair.iloc[:, 1].apply(extract_visible_text)
    ligand_symbols_col = [col for col in gene_pair.columns if "Ligand Symbols" in col][0]
    receptor_symbols_col = [col for col in gene_pair.columns if "Receptor Symbols" in col][0]
    gene_pair = gene_pair.rename(columns={ligand_symbols_col: "Human Ligand Symbols",
                                          receptor_symbols_col: "Human Receptor Symbols"})
    gene_pair_species = gene_pair_species.rename(columns={f"{species} evidence": "Evidence"})

    if species in ["Mouse", "Rat", "Frog", "Zebrafish"]:
        gene_pair_species = gene_pair_species[[
            "LR Pair Card",
            f"{species}_ligand",
            f"{species}_receptor",
            "Evidence",
            f"{species_id} ligand",
            f"{species_id} receptor",
            "ENSEMBL ligand",
            "ENSEMBL receptor",
            "PMID"
        ]]
        gene_pair_species = gene_pair_species.rename(columns={
            "ENSEMBL ligand": "Ligand ENSEMBL ID",
            "ENSEMBL receptor": "Receptor ENSEMBL ID"
                                         }
                                )
    else:
        gene_pair_species = gene_pair_species[[
            "LR Pair Card",
            f"{species}_ligand",
            f"{species}_receptor",
            "Evidence",
            f"{species_id} ligand",
            f"{species_id} receptor",
            "PMID"
        ]]
    

    
    gene_pair =gene_pair_species.merge(gene_pair,how="left", on="LR Pair Card")
    
    if species == "Mouse":
        spec_id = f"{species_id} Marker Accession ID"
        spec_name = "Marker Name"
        spec_alias = "Aliases"
    elif species == "Rat":
        spec_id = f"GENE_{species_id}_ID"
        spec_name = "NAME"
        spec_alias = "ALIASES"
    elif species == "Zebrafish":
        spec_id = f"{species_id}_ID"
        spec_name = "Current Name"
        spec_alias = "Aliases"
    elif species == "Frog":
        spec_id = "tropicalis gene ID"
        spec_name = "gene name"
        spec_alias = "gene synonyms"
    else:
        spec_id = "ensembl_gene_id"
        spec_name = "external_gene_name"
        spec_alias = "external_synonym"
    
    
    species_info = species_info[[spec_id, spec_name,spec_alias]]
    
    if species in ["Mouse", "Frog"]:
        species_info[spec_alias] = species_info[spec_alias].str.replace("|", ", ", regex=False)
    elif species == "Rat":
        species_info[spec_alias] = species_info[spec_alias].str.replace(";", ", ", regex=False)

    
    gene_pair = gene_pair.merge(species_info,how="left", left_on = f"{species_id} ligand",right_on=spec_id)
    gene_pair = gene_pair.drop(columns=[spec_id])
    gene_pair = gene_pair.rename(columns={
                                          spec_name: "Ligand Name",
                                          spec_alias: "Ligand Symbols"
                                         }
                                )
    gene_pair = gene_pair.merge(
        species_info,
        how="left",
        left_on=f"{species_id} receptor",
        right_on=spec_id
    )
    
    gene_pair = gene_pair.drop(columns=[spec_id])
    gene_pair = gene_pair.rename(columns={
                                          spec_name: "Receptor Name",
                                          spec_alias: "Receptor Symbols"
                                         }
                                )
    
    
    
    gene_pair["LR Pair"] = np.where(
        gene_pair["Evidence"] == "not conserved", 
        f"no {species_lower} ortholog",                                  
        gene_pair[f"{species}_ligand"] + " " + gene_pair[f"{species}_receptor"] 
    )
    
    gene_pair = gene_pair[~(gene_pair["Evidence"] == "not conserved")]
            
    def format_symbol_aliases(symbol, aliases):
        """
        Formats symbol, old symbols, and aliases.
        If the final formatted string would be empty after considering N/A values
        and empty inputs, it returns "species-specific".
        Otherwise, it formats based on the presence of old_symbol and aliases,
        removing unnecessary parentheses or commas, following the structure:
        "Symbol (Old Symbol, Aliases)" if both exist.
        """
        # Normalize inputs to empty strings if they are None/NaN or just whitespace
        symbol_str = str(symbol).strip()
        # old_symbol_str = str(old_symbol).strip()
        aliases_str = str(aliases).strip()
    
        # Filter out values that are empty strings or "N/A" for old_symbol and aliases
        parts_for_join = []
        # if old_symbol_str and old_symbol_str != "N/A":
        #     parts_for_join.append(old_symbol_str)
        if aliases_str and aliases_str != "N/A":
            parts_for_join.append(aliases_str)
    
        # Construct the preliminary result based on your original logic:
        # "symbol (old_symbol, aliases)" if parts_for_join is not empty, else "symbol"
        if parts_for_join:
            prelim_result = f"{symbol_str} ({', '.join(parts_for_join)})"
        else:
            prelim_result = symbol_str # Just the symbol if no old_symbol or aliases
    
        return prelim_result
    
    # This is crucial for consistent handling by the function before processing "N/A".
    gene_pair[f"{species}_ligand"] = gene_pair[f"{species}_ligand"].fillna('')
    gene_pair['Ligand Symbols'] = gene_pair['Ligand Symbols'].fillna('')
    
    gene_pair['Ligand Symbols'] = gene_pair.apply(
        lambda row: format_symbol_aliases(row[f"{species}_ligand"], row['Ligand Symbols']),
        axis=1
    )
    
    gene_pair[f"{species}_receptor"] = gene_pair[f"{species}_receptor"].fillna('')
    gene_pair['Receptor Symbols'] = gene_pair['Receptor Symbols'].fillna('')
    
    gene_pair['Receptor Symbols'] = gene_pair.apply(
        lambda row: format_symbol_aliases(row[f"{species}_receptor"], row['Receptor Symbols']),
        axis=1
    )

    # The list of columns to group by
    interaction_id_col = [col for col in gene_pair.columns if "Interaction ID" in col][0]
    
    grouping_cols = [
        interaction_id_col, "LR Pair" #, "Ligand Symbols", "Receptor Symbols" 
    ]
    
    aggregation_cols = [
        col for col in gene_pair.columns if col not in grouping_cols
    ]
    
    # 3. Create a dictionary mapping each aggregation column to the joining function
    agg_dict = {
        col: lambda x: ', '.join(x.astype(str).unique()) for col in aggregation_cols
    }
    
    # 4. Perform the groupby and aggregation
    gene_pair = gene_pair.groupby(grouping_cols).agg(agg_dict).reset_index()
    
    # make direct, conservation and conservation, direct the same
    gene_pair["Evidence"] = np.where(
    gene_pair["Evidence"].str.contains("DIRECT", na=False),
    "Direct",
    np.where(
        gene_pair["Evidence"] == "CONSERVATION",
        "Inferred",
        gene_pair["Evidence"]
    )
)
    
    
    def generate_perplexity_link_pmid(row, species, species_lower): 
        query = (
            f"What-is-the-biological-relevance-of-the-ligand-and-receptor-pair-"
            f"{row['LR Pair']}-based-on-Pubmed-ID-"
            f"{row['PMID']}-in-{species_lower}"
        )
        return (
             f'<a href="https://www.perplexity.ai/search?q={query}" target="_blank" style="text-decoration: none;">&#128172;</a>'
        )
    
    
    # Apply function to the DataFrame
    gene_pair["A.I. summary"] = gene_pair.apply(
        generate_perplexity_link_pmid, axis=1, args=(species, species_lower)
    )
    
    
    gene_pair = gene_pair.rename(columns={
                                          f"{species_id} ligand": f"Ligand {species_id} ID",
                                          f"{species_id} receptor": f"Receptor {species_id} ID"
                                         }
                                )
    ligand_loc_col = [col for col in gene_pair.columns if "Ligand Location" in col][0]
    receptor_loc_col = [col for col in gene_pair.columns if "Receptor Location" in col][0]
    lr_pair_card = [col for col in gene_pair.columns if ">LR Pair Card" in col][0]
    # gene_pair.columns
    gene_pair["LR Pair Card"] = gene_pair[lr_pair_card]
    if species in ["Mouse", "Rat", "Frog", "Zebrafish"]:
        gene_pair = gene_pair[[interaction_id_col, "LR Pair Card", "LR Pair", "Evidence", "A.I. summary", 'Ligand Symbols', 'Receptor Symbols', f"Ligand {species_id} ID", f"Receptor {species_id} ID", "Ligand ENSEMBL ID", "Receptor ENSEMBL ID",  "Human Ligand Symbols", "Human Receptor Symbols",ligand_loc_col, receptor_loc_col]]
    else:     
        gene_pair = gene_pair[[interaction_id_col, "LR Pair Card", "LR Pair", "Evidence", "A.I. summary", 'Ligand Symbols', 'Receptor Symbols',  f"{species}_ligand", f"{species}_receptor", f"Ligand {species_id} ID", f"Receptor {species_id} ID", "Human Ligand Symbols", "Human Receptor Symbols", ligand_loc_col, receptor_loc_col]]
        
    if species == "Mouse":
        # Linkify multiple species IDs in Ligand column
        gene_pair[f"Ligand {species_id} ID"] = gene_pair[f"Ligand {species_id} ID"].apply(
            lambda cell: ", ".join(
                f'<a href="https://www.informatics.jax.org/marker/{mgi.strip()}" target="_blank">{mgi.strip()}</a>'
                for mgi in str(cell).split(", ")
                if mgi.strip()
            ) if pd.notna(cell) else ""
        )
        
        # Linkify multiple MGI IDs in Receptor column
        gene_pair[f"Receptor {species_id} ID"] = gene_pair[f"Receptor {species_id} ID"].apply(
            lambda cell: ", ".join(
                f'<a href="https://www.informatics.jax.org/marker/{mgi.strip()}" target="_blank">{mgi.strip()}</a>'
                for mgi in str(cell).split(", ")
                if mgi.strip()
            ) if pd.notna(cell) else ""
        )
        
    elif species == "Rat":
        # Linkify multiple RGD IDs in Receptor column
        gene_pair[f"Ligand {species_id} ID"] = gene_pair["Ligand RGD ID"].apply(
            lambda cell: ", ".join(
                f'<a href="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id={rgd.strip()}" target="_blank">{"RGD:"+rgd.strip()}</a>'
                for rgd in str(cell).split(", ")
                if rgd.strip()
            ) if pd.notna(cell) else ""
        )
        
        
        # Linkify multiple RGD IDs in Receptor column
        gene_pair[f"Receptor {species_id} ID"] = gene_pair["Receptor RGD ID"].apply(
            lambda cell: ", ".join(
                f'<a href="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id={rgd.strip()}" target="_blank">{"RGD:"+rgd.strip()}</a>'
                for rgd in str(cell).split(", ")
                if rgd.strip()
            ) if pd.notna(cell) else ""
        )
        
    elif species == "Zebrafish":
        # Linkify multiple Zebrafish Receptor IDs
        gene_pair[f"Ligand {species_id} ID"] = gene_pair['Ligand ZFIN ID'].apply(
            lambda cell: ", ".join(
                f'<a href="https://zfin.org/{zfin.strip()}" target="_blank">{zfin.strip()}</a>'
                for zfin in str(cell).split(", ")
                if zfin.strip()
            ) if pd.notna(cell) else ""
        )
        
        # Linkify multiple Zebrafish Ligand IDs
        gene_pair[f"Receptor {species_id} ID"] = gene_pair['Receptor ZFIN ID'].apply(
            lambda cell: ", ".join(
                f'<a href="https://zfin.org/{zfin.strip()}" target="_blank">{zfin.strip()}</a>'
                for zfin in str(cell).split(", ")
                if zfin.strip()
            ) if pd.notna(cell) else ""
        )
    elif species == "Frog":
        def make_xenbase_link(cell):
            links = []
            for xid in str(cell).split(","):
                xid = xid.strip()
                if xid.startswith("XB-GENE-"):
                    url = f"https://www.xenbase.org/xenbase/gene/showgene.do?method=display&geneId={xid}"
                    links.append(f'<a href="{url}" target="_blank">{xid}</a>')
            return ", ".join(links)

        gene_pair[f"Ligand {species_id} ID"] = gene_pair['Ligand XEN ID'].apply(make_xenbase_link)
        gene_pair[f"Receptor {species_id} ID"] = gene_pair['Receptor XEN ID'].apply(make_xenbase_link)
        
    def make_ens_link(cell):
        links = []
        for eid in str(cell).split(","):
            eid = eid.strip()
            if eid.startswith("ENS"):
                url = f"http://www.ensembl.org/id/{eid}"
                links.append(f'<a href="{url}" target="_blank">{eid}</a>')
        return ", ".join(links)

    gene_pair[f"Ligand ENSEMBL ID"] = gene_pair['Ligand ENSEMBL ID'].apply(make_ens_link)
    gene_pair[f"Receptor ENSEMBL ID"] = gene_pair['Receptor ENSEMBL ID'].apply(make_ens_link)

    interactionID_cols = [col for col in gene_pair.columns if 'Interaction ID' in col][0]
    paircard_cols = [col for col in gene_pair.columns if 'LR Pair Card' in col][0]
    # Apply the transformation
    gene_pair[interactionID_cols] = gene_pair.apply(
        lambda row: replace_link_text(row[interactionID_cols], row[paircard_cols]),
        axis=1
    )
    # Drop the old LR Pair Card column
    gene_pair = gene_pair.drop(columns=[paircard_cols])
    ### tooltips 
    gene_pair["Ligand Symbols"] = [
        f'<span title="{aliases}">{aliases}</span>'
        for aliases in gene_pair["Ligand Symbols"]
    ]
    gene_pair["Receptor Symbols"] = [
        f'<span title="{aliases}">{aliases}</span>'
        for aliases in gene_pair["Receptor Symbols"]
    ]
        # Rename Species Ligand/Receptor
    gene_pair = gene_pair.rename(columns={
        f"{species}_ligand": "Ligand reserved",
        f"{species}_receptor": "Receptor reserved"
                                         }
                                )
    # Change the tooltips
    gene_pair.columns = [
    f'<span title="Ligand-receptor Pair">{col}</span>' if col == "LR Pair" else
    f'<span title="HGNC gene symbol for the ligand">{col}</span>' if col == "Ligand reserved" else
    f'<span title="HGNC gene symbol for the receptor">{col}</span>' if col == "Receptor reserved" else
     f'<span title="Official gene symbol (aliases, old names)">{col}</span>' if col in ["Ligand Symbols", "Receptor Symbols"] else
    f'<span title="HGNC gene symbols (aliases, old names)">{col}</span>' if col in ["Ligand Symbols", "Receptor Symbols"] else
    f'<span title="Click the icon below to run Perplexity on the LR pair">{col}</span>' if col == "A.I. summary" else
    f'<span title="Official Gene Symbol; Hover on symbols below to show gene names">{col}</span>' if col in ["Ligand", "Receptor"] else
    f'<span title="HGNC gene ID for the ligand (link to HGNC)">{col}</span>' if col == "Ligand HGNC ID" else
    
    f'<span title="HGNC gene ID for the receptor (link to HGNC)">{col}</span>' if col == "Receptor HGNC ID" else
    f'<span title="ENSEMBL gene ID for the ligand (link to ENSEMBL)">{col}</span>' if col == "Ligand ENSEMBL ID" else
    f'<span title="ENSEMBL gene ID for the receptor (link to ENSEMBL)">{col}</span>' if col == "Receptor ENSEMBL ID" else
    f'<span title=" PubMed IDs (PMID) with Literature Evidence for LR Interaction. Click on the link for more details">{col}</span>' if col == "PMID" else
    f'<span title="Xenbase ID (link to XEN). Click on the link for more details">{col}</span>' if col in ["Ligand XEN ID", "Receptor XEN ID"] else
    f'<span title="Rat Genome Database ID (link to RGD). Click on the link for more details">{col}</span>' if col in ["Ligand RGD ID", "Receptor RGD ID"] else
    f'<span title="Mouse Genome Informatics ID (link to MGI)">{col}</span>' if col in ["Ligand MGI ID", "Receptor MGI ID"]else
    f'<span title="Zebrafish Information Network ID (link to ZFIN)">{col}</span>' if col in ["Ligand ZFIN ID", "Receptor ZFIN ID"] else
    f'<span title="Direct: experimentally verified; Conservation: inferred from orthology">{col}</span>' if col == "Human evidence" else
    col
    for col in gene_pair.columns
]
        
    return gene_pair


mouse_gene_pair1 = process_species_gene_pair("Mouse", fetchGSheet, gene_pair)
rat_gene_pair1 = process_species_gene_pair("Rat", fetchGSheet, gene_pair)
zebrafish_gene_pair1 = process_species_gene_pair("Zebrafish", fetchGSheet, gene_pair)
frog_gene_pair1 = process_species_gene_pair("Frog", fetchGSheet, gene_pair)
chicken_gene_pair1 = process_species_gene_pair("Chicken", fetchGSheet, gene_pair)
macaque_gene_pair1 = process_species_gene_pair("Macaque", fetchGSheet, gene_pair)
pig_gene_pair1 = process_species_gene_pair("Pig", fetchGSheet, gene_pair)
dog_gene_pair1 = process_species_gene_pair("Dog", fetchGSheet, gene_pair)
cow_gene_pair1 = process_species_gene_pair("Cow", fetchGSheet, gene_pair)
chimp_gene_pair1 = process_species_gene_pair("Chimp", fetchGSheet, gene_pair)
horse_gene_pair1 = process_species_gene_pair("Horse", fetchGSheet, gene_pair)
marmoset_gene_pair1 = process_species_gene_pair("Marmoset", fetchGSheet, gene_pair)
sheep_gene_pair1 = process_species_gene_pair("Sheep", fetchGSheet, gene_pair)