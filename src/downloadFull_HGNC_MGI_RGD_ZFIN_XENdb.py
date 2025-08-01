# download the data from Mouse Genome Informatics (MGI), Rat Genome Database (RGD), ZFIN databases

import pandas as pd
import requests
from io import StringIO
from itertools import product
import csv
import os
from datetime import datetime

# today's date for version control
today = datetime.now().strftime("%Y%m%d")  # e.g., '20250723'


### HGNC database

# --- 1. Download HGNC complete txt file ---
# Define the URL and the destination file path
url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
temp_destfile = "hgnc_complete_set.txt"
final_destfile = "data/HGNC_gene_info_full.tsv"
final_destfile_today = f"data/HGNC_gene_info_full_{today}.tsv"
print(f"Downloading file from {url}...")

# Send a GET request to the URL
response = requests.get(url)

# Check if the request was successful (status code 200)
if response.status_code == 200:
    # Open the destination file in write-binary mode and write the content
    with open(temp_destfile, 'wb') as f:
        f.write(response.content)
    print(f"File successfully downloaded to {temp_destfile}")

    # --- 2. Read the file into a pandas DataFrame ---
    # The 'read_csv' function can handle tab-separated files by specifying the separator
    hgnc_data = pd.read_csv(temp_destfile, sep='\t', header=0, low_memory=False)
    # print("First 5 rows of the data:")
    # print(hgnc_data.head())

    # --- 3. Save the DataFrame to a new file ---
    # Create the 'data' directory if it doesn't exist
    os.makedirs(os.path.dirname(final_destfile), exist_ok=True)
    os.makedirs(os.path.dirname(final_destfile_today), exist_ok=True)
    # Write the DataFrame to a new tab-separated file
    # index=False is equivalent to R's row.names=FALSE
    # quote=False is the default behavior in pandas' to_csv when not specified for non-numeric types
    hgnc_data.to_csv(final_destfile, sep='\t', index=False, quotechar='"', quoting=1) # quoting=1 is csv.QUOTE_MINIMAL
    
    print(f"Data successfully saved to {final_destfile}")

    # Optional: Clean up the temporary downloaded file
    os.remove(temp_destfile)
    print(f"Removed temporary file: {temp_destfile}")

else:
    print(f"Failed to download file. Status code: {response.status_code}")


### MGI database (https://www.informatics.jax.org/downloads/reports/index.html)

import pandas as pd
import requests
from datetime import datetime
import os

# Create data folder if it doesn't exist
os.makedirs("data", exist_ok=True)

# Today's date for versioning
today = datetime.now().strftime("%Y%m%d")

# URLs and filenames
homology_url = "https://www.informatics.jax.org/downloads/reports/HGNC_AllianceHomology.rpt"
mrk_url = "https://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt"
mrk_list_url = "https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"

homology_file = f"data/HGNC_AllianceHomology_{today}.tsv"
mrk_file = f"data/MRK_Sequence_{today}.tsv"
mrk_list_file = f"data/MRK_List2_{today}.tsv"
merged_file = f"data/MRK_Merged_{today}_MGI_DB.tsv"

# Download HGNC_AllianceHomology.rpt
response = requests.get(homology_url)
response.raise_for_status()
with open(homology_file, "wb") as f:
    f.write(response.content)

# Download MRK_Sequence.rpt
response = requests.get(mrk_url)
response.raise_for_status()
with open(mrk_file, "wb") as f:
    f.write(response.content)

# Download MRK_List2.rpt
response = requests.get(mrk_list_url)
response.raise_for_status()
with open(mrk_list_file, "wb") as f:
    f.write(response.content)

# Load files
df_mrk = pd.read_csv(mrk_file, sep="\t", dtype=str)
df_homology = pd.read_csv(homology_file, sep="\t", dtype=str, index_col = False)

df_mrk_list = pd.read_csv(mrk_list_file, sep="\t", dtype=str)
df_mrk_list= df_mrk_list[["MGI Accession ID","Marker Synonyms (pipe-separated)"]].rename(columns={"MGI Accession ID": "MGI Marker Accession ID",
                                                                                                  "Marker Synonyms (pipe-separated)": "Aliases"})
df_homology = df_homology.rename(columns={"MGI Accession ID": "MGI Marker Accession ID"})
# Keep only relevant columns from HGNC_AllianceHomology
columns_to_keep = [
    "MGI Marker Accession ID",
    "EntrezGene ID",
    "CCDS IDs",
    "HGNC ID"
]
df_homology_subset = df_homology[columns_to_keep]
# Merge: MRK + homology
df_merged = df_mrk.merge(df_homology_subset, how="left", on="MGI Marker Accession ID")

# Merge with synonyms
df_merged = df_merged.merge(df_mrk_list, how="left", on="MGI Marker Accession ID")

# Save result
df_merged.to_csv(merged_file, sep="\t", index=False)
print(f"Merged file saved to: {merged_file}")

### RGD database (https://download.rgd.mcw.edu/data_release/)

# 1. Define URL and dated filename
import csv
url = "https://download.rgd.mcw.edu/data_release/GENES_RAT.txt"
out_file = f"data/GENES_RAT_RGD_DB.tsv"
out_file_today = f"data/GENES_RAT_{today}_RGD_DB.tsv"
# 2. Fetch raw content as text
response = requests.get(url)
response.raise_for_status()
content = response.text

# 3. Split into lines, remove first 102 header lines
lines = content.splitlines()[102:]

# 4. Parse into rows
rows = [line.split("\t") for line in lines]
header = rows[1]
data_rows = rows[2:]


# 5. Write as proper TSV
with open(out_file_today, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header)
    writer.writerows(data_rows)

print(f"✅ Written {len(rows)} rows to {out_file_today}")

with open(out_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header)
    writer.writerows(data_rows)


### ZFIN database (https://zfin.org/downloads)

# Define file URLs and your exact header specifications
files = {
    "aliases": {
        "url": "https://zfin.org/downloads/aliases.txt",
        "headers": ["Current ZFIN ID", "Current Name", "Current Symbol",
                    "Previous Name", "SO ID"],
        "key": "Current ZFIN ID"
    },
    "gene": {
        "url": "https://zfin.org/downloads/gene.txt",
        "headers": ["ZFIN ID", "SO ID", "Symbol", "NCBI Gene ID"],
        "key": "ZFIN ID"
    },
    "orthos": {
        "url": "https://zfin.org/downloads/human_orthos.txt",
        "headers": ["ZFIN ID", "ZFIN Symbol", "ZFIN Name", "Human Symbol",
                    "Human Name", "OMIM ID", "Gene ID", "HGNC ID",
                    "Evidence", "Pub ID", "ZFIN Abbreviation Name",
                    "ECO ID", "ECO Term Name"],
        "key": "ZFIN ID"
    }
}

# 1. Download and load each file with your defined headers
dfs = {}
for name, params in files.items():
    r = requests.get(params["url"])
    r.raise_for_status()
    df = pd.read_csv(
        StringIO(r.text),
        sep="\t",
        comment="#",
        header=None,
        names=params["headers"],
        dtype=str
    )
    # Standardize the key column name
    df = df.rename(columns={params["key"]: "ZFIN_ID"})
    dfs[name] = df

# 2. Merge on ZFIN_ID using outer join
merged = dfs["aliases"]
for name in ["gene", "orthos"]:
    merged = merged.merge(dfs[name], on="ZFIN_ID", how="outer", suffixes=("", f"_{name}"))

# 3. Drop duplicated columns (identical names and contents)
merged = merged.loc[:, ~merged.columns.duplicated()]
# Also drop perfectly identical-content duplicates
columns = merged.columns
to_drop = [
    col2 for i, col1 in enumerate(columns)
    for col2 in columns[i+1:]
    if merged[col1].equals(merged[col2])
]
merged = merged.drop(columns=to_drop)
merged = merged.drop(columns=['SO ID_gene', 'Symbol', 'Evidence','Pub ID', 'ZFIN Symbol', 'ZFIN Name'])
merged = merged.drop_duplicates()

# 4. Save 
output = f"data/Zebrafish_merged_ZFIN_DB.tsv"
output_today = f"data/Zebrafish_merged_{today}_ZFIN_DB.tsv"
merged.to_csv(output, sep="\t", index=False, encoding="utf-8")
merged.to_csv(output_today, sep="\t", index=False, encoding="utf-8")

print(f"✅ Merged file saved: {output_today} — shape: {merged.shape}")

### Xenbase database (https://www.xenbase.org/entry/)
# Define Xenbase URL and output files
### Xenbase database (https://www.xenbase.org/entry/)

# Define Xenbase URL and output files
xenbase_url = "https://download.xenbase.org/xenbase/GenePageReports/GenePageGeneralInfo_ManuallyCurated.txt"
xenbase_file = f"data/GenePageGeneralInfo_Xenbase_DB.tsv"
xenbase_file_today = f"data/GenePageGeneralInfo_{today}_Xenbase_DB.tsv"

# Download the file
response = requests.get(xenbase_url)
response.raise_for_status()

# Manually assign headers since the file has no header row
full_header = [
    "Xenbase genepage ID",      # 0
    "gene symbol",              # 1
    "gene name",                # 2
    "gene function",            # 3
    "gene synonyms",            # 4
    "JGI ID",                   # 5
    # ... The file likely has more columns after this, which we will ignore
]

# Read the file without header and assign column names
# Define Xenbase URL and output files
xenbase_url = "https://download.xenbase.org/xenbase/GenePageReports/GenePageGeneralInfo_ManuallyCurated.txt"
xenbase_file = f"data/GenePageGeneralInfo_Xenbase_DB.tsv"
xenbase_file_today = f"data/GenePageGeneralInfo_{today}_Xenbase_DB.tsv"

# Download the file
response = requests.get(xenbase_url)
response.raise_for_status()

# Manually assign headers since the file has no header row
full_header = [
    "Xenbase genepage ID",      # 0
    "gene symbol",              # 1
    "gene name",                # 2
    "gene function",            # 3
    "gene synonyms",            # 4
    "JGI ID",                   # 5
    # ... The file likely has more columns after this, which we will ignore
]

# Read the file without header and assign column names
df_xenbase = pd.read_csv(
    StringIO(response.text),
    sep="\t",
    header=None,
    names=full_header + [f"extra_{i}" for i in range(6, 20)],  # pad extra dummy headers
    dtype=str
)

# Keep only the desired columns
df_xenbase = df_xenbase[full_header]

xenbase_url_add = "https://download.xenbase.org/xenbase/GenePageReports/XenbaseGenepageToGeneIdMapping.txt"
xenbase_file = f"data/GenePageGeneralInfo_Xenbase_DB.tsv"
xenbase_file_today = f"data/GenePageGeneralInfo_{today}_Xenbase_DB.tsv"

# Download the file
response_add = requests.get(xenbase_url_add)
response_add.raise_for_status()

full_header = [
    "Xenbase genepage ID",     
    "Xenbase genepage symbol",             
    "tropicalis gene ID",              
    "tropicalis gene symbol",          
    "laevis.L gene ID",            
    "laevis.L gene symbol",   
    "laevis.S gene ID",
    "laevis.S gene symbol"
]
# Read the file without header and assign column names
df_xenbase_add = pd.read_csv(
    StringIO(response_add.text),
    sep="\t",
    header=None,
    names=full_header + [f"extra_{i}" for i in range(6, 20)],  # pad extra dummy headers
    dtype=str
)
# Keep only the desired columns
df_xenbase_add = df_xenbase_add[full_header]
df_xenbase_add = df_xenbase_add[['Xenbase genepage ID', 'tropicalis gene ID', 'tropicalis gene symbol']]
df_xenbase=df_xenbase.merge(df_xenbase_add, how="left", on="Xenbase genepage ID")

# Save to file
df_xenbase.to_csv(xenbase_file, sep="\t", index=False, encoding="utf-8")
df_xenbase.to_csv(xenbase_file_today, sep="\t", index=False, encoding="utf-8")

print(f"✅ Xenbase gene info saved: {xenbase_file_today} — shape: {df_xenbase.shape}")
