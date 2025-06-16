## Function to convert Mouse (ENSEMBL) to RGD ID via biomart
from pybiomart import Server
import pandas as pd
sys.path.append(os.path.abspath("src"))
from createDataTable import mouse_specific_mgi_ids
# Full list of unique MGI IDs
mgi_ids = mouse_specific_mgi_ids

server = Server(host='http://www.ensembl.org')

mouse_dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['mmusculus_gene_ensembl']
rat_dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['rnorvegicus_gene_ensembl']

# Get mouse Ensembl IDs and rat ortholog Ensembl IDs
mouse_genes = mouse_dataset.query(
    attributes=[
        'ensembl_gene_id',
        'mgi_id',
        'mgi_symbol'
    ]
)

mouse_genes = mouse_genes[mouse_genes["MGI ID"].isin(mgi_ids)]
# Step 2: Get mouse Ensembl genes with rat ortholog Ensembl gene IDs (no mouse MGI here)
mouse_rat_orthologs = mouse_dataset.query(
    attributes=[
        'ensembl_gene_id',
        'rnorvegicus_homolog_ensembl_gene',
        'rnorvegicus_homolog_associated_gene_name'
    ]
)

# Step 3: Merge mouse gene info with orthologs on ensembl_gene_id
mouse_to_rat = pd.merge(mouse_genes, mouse_rat_orthologs, on='Gene stable ID')
# Filter out entries without rat ortholog
mouse_to_rat = mouse_to_rat[mouse_to_rat['Norway rat - BN/NHsdMcwi gene stable ID'].notna() & (mouse_to_rat['Norway rat - BN/NHsdMcwi gene stable ID'] != '')]

# Step 4: Retrieve RGD IDs for rat ortholog Ensembl gene IDs
rat_ensembl_ids = mouse_to_rat['Norway rat - BN/NHsdMcwi gene stable ID'].unique().tolist()
rat_ensembl_to_rgd = rat_dataset.query(
    attributes=['ensembl_gene_id', 'rgd_id']
)
rat_ensembl_to_rgd = rat_ensembl_to_rgd[rat_ensembl_to_rgd["Gene stable ID"].isin(rat_ensembl_ids)]
# Step 5: Merge rat RGD IDs back
final_mapping = pd.merge(
    mouse_to_rat,
    rat_ensembl_to_rgd,
    left_on='Norway rat - BN/NHsdMcwi gene stable ID',
    right_on='Gene stable ID',
    how='left'
)
# Clean column names for clarity
final_mapping.rename(columns={
    'Gene stable ID_x': 'mouse_ensembl_gene_id',
    'Gene stable ID_y': 'rat_ensembl_gene_id',
    'Norway rat - BN/NHsdMcwi gene name': 'RGD symbol',
    'rgd_id': 'RGD ID'
}, inplace=True)

final_mapping = final_mapping.drop(columns=["Norway rat - BN/NHsdMcwi gene stable ID"])
final_mapping['RGD ID'] = final_mapping['RGD ID'].apply(lambda x: f"RGD:{int(x)}" if pd.notnull(x) else None)
# Save
final_mapping.to_csv("data/mouse_to_rat_mapping.csv", index=False)