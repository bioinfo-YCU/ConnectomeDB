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

# 1. URL and date-formatted filename
url = "https://www.informatics.jax.org/downloads/reports/HGNC_AllianceHomology.rpt"
filename = f"data/HGNC_AllianceHomology_{today}_MGI_DB.tsv"

# 2. Fetch the file
response = requests.get(url)
response.raise_for_status()  # stop if download failed

# 3. Save with TSV extension
with open(filename, "wb") as f:
    f.write(response.content)

print(f"✅ Saved as {filename}")


### RGD database (https://download.rgd.mcw.edu/data_release/)

# 1. Define URL and dated filename
import csv
url = "https://download.rgd.mcw.edu/data_release/GENES_RAT.txt"
out_file = f"data/GENES_RAT_{today}_RGD_DB.tsv"

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
with open(out_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header)
    writer.writerows(data_rows)

print(f"✅ Written {len(rows)} rows to {out_file}")


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
output = f"data/Zebrafish_merged_{today}_ZFIN_DB.tsv"
merged.to_csv(output, sep="\t", index=False, encoding="utf-8")

print(f"✅ Merged file saved: {output} — shape: {merged.shape}")
