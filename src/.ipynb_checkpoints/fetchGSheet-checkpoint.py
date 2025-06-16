import gspread
from google.oauth2.service_account import Credentials
import pandas as pd
import time
import random
from gspread.exceptions import APIError, WorksheetNotFound

sheet_ID = "1XP5wBDN_orSlE8RLb2TxSclIopVO1mb_1S3ENf2qYSw"
credentials_file = 'data/connectomedb2025-a9acdf562a84.json'

def fetch_google_sheet_data(sheet_ID, tab_name, credentials_file, max_retries=5):
    """
    Fetch data from a specific tab of a Google Sheet with retry/backoff on quota errors.
    """
    SCOPES = ['https://www.googleapis.com/auth/spreadsheets', 'https://www.googleapis.com/auth/drive']
    credentials = Credentials.from_service_account_file(credentials_file, scopes=SCOPES)
    client = gspread.authorize(credentials)
    
    for attempt in range(max_retries):
        try:
            sheet = client.open_by_key(sheet_ID)
            worksheet = sheet.worksheet(tab_name)
            data = worksheet.get_all_values()
            df = pd.DataFrame(data)
            df.columns = df.iloc[0]
            df = df[1:].reset_index(drop=True)
            return df
        except WorksheetNotFound:
            raise ValueError(f"Tab '{tab_name}' not found in Google Sheet with ID '{sheet_ID}'.")
        except APIError as e:
            if e.response.status_code == 429:
                wait_time = (2 ** attempt) + random.uniform(0, 1)
                print(f"[429 Error] Rate limit hit. Retrying in {wait_time:.2f} seconds...")
                time.sleep(wait_time)
            else:
                raise
    raise RuntimeError(f"Failed to fetch data after {max_retries} retries due to API quota limits.")

# ADD RATE LIMITING BETWEEN CALLS - This is the key fix!
def rate_limited_fetch(sheet_ID, tab_name, credentials_file, delay=1.5):
    """Wrapper that adds delay between API calls"""
    print(f"Fetching '{tab_name}'...")
    result = fetch_google_sheet_data(sheet_ID, tab_name, credentials_file)
    print(f"Successfully fetched {len(result)} rows from '{tab_name}'")
    time.sleep(delay)  # Wait between calls to avoid quota issues
    return result

# Your existing code with rate limiting added:
print("Starting data fetch with rate limiting...")

# Fetching data from Google Sheets
gene_pair = rate_limited_fetch(sheet_ID, "FROZEN LIST HUMAN", credentials_file)
gene_pair_mouse = rate_limited_fetch(sheet_ID, "FROZEN LIST MOUSE", credentials_file)

# For now append mouse and remove in human pair later
gene_pair = pd.concat([gene_pair, gene_pair_mouse])

# Ligand and receptor location
ligand_loc = rate_limited_fetch(sheet_ID, "Ligand_location_HUMAN", credentials_file)
receptor_loc = rate_limited_fetch(sheet_ID, "Receptor_location_HUMAN", credentials_file)

# For now append mouse and remove in human pair later
ligand_loc_mouse = rate_limited_fetch(sheet_ID, "Ligand_location_Mouse", credentials_file)
ligand_loc = pd.concat([ligand_loc, ligand_loc_mouse])

receptor_loc_mouse = rate_limited_fetch(sheet_ID, "Receptor_location_Mouse", credentials_file)
receptor_loc = pd.concat([receptor_loc, receptor_loc_mouse])

# Pathways
kegg_pathway_info = rate_limited_fetch(sheet_ID, "KEGG_metadata_pairs in frozen", credentials_file)

# HGNC gene group
gene_group = rate_limited_fetch(sheet_ID, "HGNC gene group", credentials_file)

src_info = rate_limited_fetch(sheet_ID, "sourceAbbv", credentials_file)

# Your existing processing code stays the same:
human_gene_pair = gene_pair.iloc[:, :-36]
# remove mouse info
human_gene_pair = human_gene_pair.iloc[:-13]