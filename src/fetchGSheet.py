# fetchGSheet.py - Module that creates variables at module level for import

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
            if len(df) > 0:
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

# Rate limiting function
def safe_fetch(sheet_ID, tab_name, credentials_file, delay=1.5):
    """Fetch with rate limiting to avoid quota issues"""
    try:
        result = fetch_google_sheet_data(sheet_ID, tab_name, credentials_file)
        time.sleep(delay)  # Rate limiting delay
        return result
    except Exception as e:
        # Return empty DataFrame so other code doesn't break
        return pd.DataFrame()

# ============================================================================
# MODULE LEVEL VARIABLES - These are what other files import
# ============================================================================

# Loading Google Sheets data silently

# Fetching and combining human/mouse data
gene_pair_human = safe_fetch(sheet_ID, "FROZEN_human", credentials_file)
gene_pair_mouse = safe_fetch(sheet_ID, "FROZEN_mouse", credentials_file)
numOfMouseOrth = len(gene_pair_mouse[0:])
# Combine human and mouse gene pairs
if not gene_pair_human.empty and not gene_pair_mouse.empty:
    gene_pair = pd.concat([gene_pair_human, gene_pair_mouse])
elif gene_pair_mouse.empty and not gene_pair_human.empty:
    gene_pair = gene_pair_human  # Keep just human data
elif gene_pair_human.empty and not gene_pair_mouse.empty:
    gene_pair = gene_pair_mouse  # Use mouse data if human fails
else:
    gene_pair = pd.DataFrame()  # Both failed

# Ligand location data
ligand_loc_human = safe_fetch(sheet_ID, "Ligand_location_HUMAN", credentials_file)
ligand_loc_mouse = safe_fetch(sheet_ID, "Ligand_location_Mouse", credentials_file)



# Combine ligand location data
if not ligand_loc_human.empty and not ligand_loc_mouse.empty:
    ligand_loc = pd.concat([ligand_loc_human, ligand_loc_mouse])
elif ligand_loc_mouse.empty and not ligand_loc_human.empty:
    ligand_loc = ligand_loc_human
elif ligand_loc_human.empty and not ligand_loc_mouse.empty:
    ligand_loc = ligand_loc_mouse
else:
    ligand_loc = pd.DataFrame()

# Receptor location data
receptor_loc_human = safe_fetch(sheet_ID, "Receptor_location_HUMAN", credentials_file)
receptor_loc_mouse = safe_fetch(sheet_ID, "Receptor_location_Mouse", credentials_file)

# Combine receptor location data
if not receptor_loc_human.empty and not receptor_loc_mouse.empty:
    receptor_loc = pd.concat([receptor_loc_human, receptor_loc_mouse])
elif receptor_loc_mouse.empty and not receptor_loc_human.empty:
    receptor_loc = receptor_loc_human
elif receptor_loc_human.empty and not receptor_loc_mouse.empty:
    receptor_loc = receptor_loc_mouse
else:
    receptor_loc = pd.DataFrame()

# Other data
kegg_pathway_info = safe_fetch(sheet_ID, "KEGG_metadata_pairs in frozen", credentials_file)
gene_group = safe_fetch(sheet_ID, "HGNC gene group", credentials_file)
src_info = safe_fetch(sheet_ID, "sourceAbbv", credentials_file)

# Process human gene pairs
# if not gene_pair.empty and gene_pair.shape[1] > 36:
#     human_gene_pair = gene_pair.iloc[:, :-36]
#     if len(human_gene_pair) > 13:
#         human_gene_pair = human_gene_pair.iloc[:-13]
#     # Removed warning print
# else:
#     # Removed warning print
#     human_gene_pair = pd.DataFrame()

# Data loading complete - removed final print statement

# Make only the final merged variables available for import
__all__ = [
    'gene_pair',           # Combined human + mouse
    'ligand_loc',          # Combined human + mouse  
    'receptor_loc',        # Combined human + mouse
    'kegg_pathway_info',   # Single sheet
    'gene_group',          # Single sheet
    'src_info',            # Single sheet
    'gene_pair_human',      # Processed gene pairs
    'gene_pair_mouse'
]