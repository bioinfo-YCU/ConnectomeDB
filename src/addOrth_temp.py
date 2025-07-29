# Function to add species-specific species Enseml ID and symbol for all other species except for mouse, rat, and zebrafish
def appendOtherSpeciesInfo(species, origDF):
    species_name = {
        "mmusculus": "Mouse",
        "rnorvegicus": "Rat",
        "drerio": "Zebrafish",
        "ptroglodytes": "Chimpanzee",
        "ggallus": "Chicken",
        "sscrofa": "Pig",
        "btaurus": "Cow",
        "clfamiliaris": "Dog",
        "ecaballus": "Horse",
        "oarambouillet": "Sheep",
        "cjacchus": "Marmoset",
        "mmulatta": "Macaque"    
    }.get(species, "Unknown species")

    # Load species-specific data
    species_info = pd.read_csv(f"data/{species}_ID_biomart.csv")

    # Keep relevant columns - use species code, not species_name
    species_mapping = {
        "mmusculus": "mgi_id",
        "rnorvegicus": "rgd_id", 
        "drerio": "zfin_id_id"
    }
    
    # Fix: Use species (not species_name) for lookup
    species_id = species_mapping.get(species, f"{species}_homolog_ensembl_gene")
    
    # Check what columns actually exist in the CSV
    available_cols = [col for col in [species_id, f"{species}_homolog_associated_gene_name", 'hgnc_id'] 
                     if col in species_info.columns]
    
    species_info = species_info[available_cols]

    # Remove rows where 'hgnc_id' is NaN and handle groupby properly
    if 'hgnc_id' in species_info.columns:
        species_info = species_info.dropna(subset=['hgnc_id'])
        # Fix the groupby - aggregate non-hgnc_id columns properly
        species_info = species_info.groupby('hgnc_id').agg({
            col: lambda x: ', '.join(x.dropna().astype(str).unique()) 
            for col in species_info.columns if col != 'hgnc_id'
        }).reset_index()

    # Merge with ligand data
    origDF = origDF.merge(species_info, how='left', 
                           left_on='Ligand HGNC ID', right_on='hgnc_id',
                           suffixes=('', '_lig'))
    
    # Rename columns for ligand info
    rename_dict = {}
    if f"{species}_homolog_associated_gene_name" in origDF.columns:
        rename_dict[f"{species}_homolog_associated_gene_name"] = f"{species_name} Ligand"
    if f"{species}_homolog_ensembl_gene" in origDF.columns:
        rename_dict[f"{species}_homolog_ensembl_gene"] = f"{species_name} Ligand Ensembl ID"
    if species_id in origDF.columns and species_id != f"{species}_homolog_ensembl_gene":
        rename_dict[species_id] = f"{species_name} Ligand ID"
        
    origDF = origDF.rename(columns=rename_dict)

    # Drop hgnc_id column
    if 'hgnc_id' in origDF.columns:
        origDF = origDF.drop(columns=['hgnc_id'])

    # Merge with receptor data
    origDF = origDF.merge(species_info, how='left', 
                           left_on='Receptor HGNC ID', right_on='hgnc_id',
                           suffixes=('', '_rec'))

    # Rename columns for receptor info
    rename_dict = {}
    if f"{species}_homolog_associated_gene_name" in origDF.columns:
        rename_dict[f"{species}_homolog_associated_gene_name"] = f"{species_name} Receptor"
    elif f"{species}_homolog_associated_gene_name_rec" in origDF.columns:
        rename_dict[f"{species}_homolog_associated_gene_name_rec"] = f"{species_name} Receptor"
        
    if f"{species}_homolog_ensembl_gene" in origDF.columns:
        rename_dict[f"{species}_homolog_ensembl_gene"] = f"{species_name} Receptor Ensembl ID"
    elif f"{species}_homolog_ensembl_gene_rec" in origDF.columns:
        rename_dict[f"{species}_homolog_ensembl_gene_rec"] = f"{species_name} Receptor Ensembl ID"
        
    if species_id in origDF.columns and species_id != f"{species}_homolog_ensembl_gene":
        rename_dict[species_id] = f"{species_name} Receptor ID"
    elif f"{species_id}_rec" in origDF.columns:
        rename_dict[f"{species_id}_rec"] = f"{species_name} Receptor ID"
        
    origDF = origDF.rename(columns=rename_dict)

    # Drop hgnc_id columns
    cols_to_drop = [col for col in origDF.columns if col in ['hgnc_id', 'hgnc_id_rec']]
    if cols_to_drop:
        origDF = origDF.drop(columns=cols_to_drop)

    # Drop columns where all values are NaN
    origDF = origDF.dropna(axis=1, how='all')

    return origDF


# # Loop through each species and update gene_pair
# for species in species_list:
#     gene_pair = appendOtherSpeciesInfo(species, gene_pair)

# Drop columns where all values are NA in gene_pair
gene_pair = gene_pair.dropna(axis=1, how='all')

gene_pair = gene_pair.fillna(" ")
gene_pair = gene_pair[gene_pair['Human LR Pair'] != ' ']
# Step 1: Identify rows where 'Ligand HGNC ID' is missing or empty,
# and set 'Mouse Ligand' to match the human 'Ligand' name.
mask = gene_pair['Ligand HGNC ID'].astype(str).str.strip() == ''
gene_pair.loc[mask, 'Mouse Ligand'] = gene_pair.loc[mask, 'Ligand']

# Step 2: Load MGI info table and keep only relevant, non-null mappings
MGI_info = pd.read_csv("data/mouse_name_mapping.csv")
MGI_info = MGI_info[["MGI symbol", "MGI ID"]].dropna()
# Strip leading/trailing whitespace
gene_pair['Mouse Ligand'] = gene_pair['Mouse Ligand'].astype(str).str.strip()
MGI_info['MGI symbol'] = MGI_info['MGI symbol'].astype(str).str.strip()
# Step 3: Merge to get MGI ID where 'Mouse Ligand' matches the mouse gene name
gene_pair = gene_pair.merge(
    MGI_info,
    left_on='Mouse Ligand',
    right_on='MGI symbol',
    how='left')

gene_pair['Mouse Ligand ID'] = gene_pair.apply(
    lambda row: row['MGI ID'] if is_mouse_specific(row['Ligand']) else row['Mouse Ligand ID'],
    axis=1
)

# Step 5: Drop temporary columns from merge
gene_pair = gene_pair.drop(columns=['MGI ID', 'MGI symbol'])

# and set 'Mouse Receptor' to match the human 'Receptor' name.
mask = gene_pair['Receptor HGNC ID'].astype(str).str.strip() == ''
gene_pair.loc[mask, 'Mouse Receptor'] = gene_pair.loc[mask, 'Receptor']

gene_pair['Mouse Receptor'] = gene_pair['Mouse Receptor'].astype(str).str.strip()
MGI_info['MGI symbol'] = MGI_info['MGI symbol'].astype(str).str.strip()
gene_pair = gene_pair.merge(
    MGI_info,
    left_on='Mouse Receptor',
    right_on='MGI symbol',
    how='left')

gene_pair['Mouse Receptor ID'] = gene_pair.apply(
    lambda row: row['MGI ID'] if is_mouse_specific(row['Receptor']) else row['Mouse Receptor ID'],
    axis=1
)

gene_pair = gene_pair.drop(columns=['MGI ID', 'MGI symbol'])


# Linkify multiple Zebrafish Receptor IDs
gene_pair['Zebrafish Receptor ID'] = gene_pair['Zebrafish Receptor ID'].apply(
    lambda cell: ", ".join(
        f'<a href="https://zfin.org/{zfin.strip()}" target="_blank">{zfin.strip()}</a>'
        for zfin in str(cell).split(", ")
        if zfin.strip()
    ) if pd.notna(cell) else ""
)

# Linkify multiple Zebrafish Ligand IDs
gene_pair['Zebrafish Ligand ID'] = gene_pair['Zebrafish Ligand ID'].apply(
    lambda cell: ", ".join(
        f'<a href="https://zfin.org/{zfin.strip()}" target="_blank">{zfin.strip()}</a>'
        for zfin in str(cell).split(", ")
        if zfin.strip()
    ) if pd.notna(cell) else ""
)

gene_pair = gene_pair.rename(columns={
                                     "Mouse Ligand ID": "Ligand MGI ID", 
                                     "Rat Ligand ID": "Ligand RGD ID",
                                     "Zebrafish Ligand ID": "Ligand ZFIN ID",
                                     "Mouse Receptor ID": "Receptor MGI ID", 
                                     "Rat Receptor ID": "Receptor RGD ID",
                                     "Zebrafish Receptor ID": "Receptor ZFIN ID",}
                            )



# add some missing MGI since they were mouse-specific
MGI_info = dict(zip(MGI_info['MGI symbol'], MGI_info['MGI ID']))
# Function to apply the mapping only if the MGI ID is empty
def map_if_empty(row, gene_col, mgi_col, mapping_dict):
    if pd.isna(row[mgi_col]) or row[mgi_col] == '':
        return mapping_dict.get(row[gene_col], '') # Use .get() to avoid KeyError if symbol not in mapping
    return row[mgi_col]

# Apply the mapping to 'Receptor MGI ID'
gene_pair['Receptor MGI ID'] = gene_pair.apply(
    lambda row: map_if_empty(row, 'Receptor', 'Receptor MGI ID', mapping),
    axis=1
)

# Apply the mapping to 'Ligand MGI ID'
gene_pair['Ligand MGI ID'] = gene_pair.apply(
    lambda row: map_if_empty(row, 'Ligand', 'Ligand MGI ID', mapping),
    axis=1
)

# add some missing RGD since they were mouse-specific
mouse_rat_info = pd.read_csv("data/mouse_to_rat_mapping.csv")
mapping_mouse_to_rat = dict(zip(mouse_rat_info['MGI ID'], mouse_rat_info['RGD ID']))
mapping_mr2 = dict(zip(mouse_rat_info['RGD ID'], mouse_rat_info['RGD symbol']))

def format_rgd_ids(cell):
    if pd.isna(cell):
        return ""
    ids = str(cell).split(",")
    cleaned = []
    for rgd in ids:
        rgd = rgd.strip()
        if not rgd or rgd.lower() == "none":
            continue
        try:
            rgd_clean = f"RGD:{int(float(rgd))}"
        except ValueError:
            rgd_clean = f"RGD:{rgd}"
        cleaned.append(rgd_clean)
    return ", ".join(cleaned)

    
gene_pair["Ligand RGD ID"] = gene_pair["Ligand RGD ID"].apply(format_rgd_ids)
gene_pair["Receptor RGD ID"] = gene_pair["Receptor RGD ID"].apply(format_rgd_ids)
# Apply the mapping to 'Ligand MGI ID'
gene_pair['Ligand RGD ID'] = gene_pair.apply(
    lambda row: map_if_empty(row, 'Ligand MGI ID', 'Ligand RGD ID', mapping_mouse_to_rat),
    axis=1
) 

# Apply the mapping to 'Receptor MGI ID'
gene_pair['Receptor RGD ID'] = gene_pair.apply(
    lambda row: map_if_empty(row, 'Receptor MGI ID', 'Receptor RGD ID', mapping_mouse_to_rat),
    axis=1
)

def map_if_empty(row, source_col, target_col, mapping_dict):
    current_val = row[target_col]
    source_val = row[source_col]

    # Don't overwrite if target already has a value
    if pd.notna(current_val) and str(current_val).strip():
        return current_val

    # Only try mapping if source is non-empty
    if pd.isna(source_val) or not str(source_val).strip():
        return ""

    # Clean source value
    key = str(source_val).strip()
    
    # ✅ If multiple IDs (comma-separated), map only the first one
    key = key.split(",")[0].strip()
    
    # ✅ Ensure it starts with "RGD:"
    if not key.startswith("RGD:"):
        try:
            key = f"RGD:{int(float(key))}"
        except ValueError:
            pass

    return mapping_dict.get(key, "")

    
gene_pair['Rat Receptor'] = gene_pair.apply(
    lambda row: map_if_empty(row, 'Receptor RGD ID', 'Rat Receptor', mapping_mr2),
    axis=1
)

gene_pair['Rat Ligand'] = gene_pair.apply(
    lambda row: map_if_empty(row, 'Ligand RGD ID', 'Rat Ligand', mapping_mr2),
    axis=1
)


# extract the mgi that would need to be converted from mouse to rat
mouse_specific_mgi_ids_ligand = gene_pair[gene_pair['Ligand'].apply(is_mouse_specific)]['Ligand MGI ID'].tolist()
mouse_specific_mgi_ids_receptor = gene_pair[gene_pair['Receptor'].apply(is_mouse_specific)]['Receptor MGI ID'].tolist()
mouse_specific_mgi_ids = mouse_specific_mgi_ids_ligand + mouse_specific_mgi_ids_receptor
mouse_specific_mgi_ids = pd.unique(pd.Series(mouse_specific_mgi_ids))


# Linkify multiple RGD IDs in Ligand column

def clean_rgd_id(rgd):
    rgd = str(rgd).strip()

# Linkify multiple MGI IDs in Ligand column
gene_pair["Ligand MGI ID"] = gene_pair["Ligand MGI ID"].apply(
    lambda cell: ", ".join(
        f'<a href="https://www.informatics.jax.org/marker/{mgi.strip()}" target="_blank">{mgi.strip()}</a>'
        for mgi in str(cell).split(", ")
        if mgi.strip()
    ) if pd.notna(cell) else ""
)

# Linkify multiple MGI IDs in Receptor column
gene_pair["Receptor MGI ID"] = gene_pair["Receptor MGI ID"].apply(
    lambda cell: ", ".join(
        f'<a href="https://www.informatics.jax.org/marker/{mgi.strip()}" target="_blank">{mgi.strip()}</a>'
        for mgi in str(cell).split(", ")
        if mgi.strip()
    ) if pd.notna(cell) else ""
)



# Linkify multiple RGD IDs in Receptor column
gene_pair["Ligand RGD ID"] = gene_pair["Ligand RGD ID"].apply(
    lambda cell: ", ".join(
        f'<a href="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id={rgd.strip()}" target="_blank">{rgd.strip()}</a>'
        for rgd in str(cell).split(", ")
        if rgd.strip()
    ) if pd.notna(cell) else ""
)


# Linkify multiple RGD IDs in Receptor column
gene_pair["Receptor RGD ID"] = gene_pair["Receptor RGD ID"].apply(
    lambda cell: ", ".join(
        f'<a href="https://rgd.mcw.edu/rgdweb/report/gene/main.html?id={rgd.strip()}" target="_blank">{rgd.strip()}</a>'
        for rgd in str(cell).split(", ")
        if rgd.strip()
    ) if pd.notna(cell) else ""
)
