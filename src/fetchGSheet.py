## Function to fetch data from a Google Sheet tab

import gspread
from google.oauth2.service_account import Credentials
import pandas as pd

sheet_ID = "1XP5wBDN_orSlE8RLb2TxSclIopVO1mb_1S3ENf2qYSw" #"15FfI7cVpJmAcytTBmhVE2Z77tgVQYfEZUEe700_Wzlg"
credentials_file = 'data/connectomedb2025-a9acdf562a84.json'

import time
import random
import pandas as pd
import gspread
from google.auth.transport.requests import Request
from google.oauth2.service_account import Credentials
from gspread.exceptions import APIError, WorksheetNotFound

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


# Fetching data from Google Sheets
gene_pair = fetch_google_sheet_data(sheet_ID, "FROZEN LIST HUMAN", credentials_file)
#loc_info = fetch_google_sheet_data(sheet_ID, "proteome_HPA", credentials_file) 
# Ligand and receptor location # previously based on localization 
ligand_loc = fetch_google_sheet_data(sheet_ID, "Ligand_location_HUMAN", credentials_file) 
receptor_loc = fetch_google_sheet_data(sheet_ID, "Receptor_location_HUMAN", credentials_file) 
# Pathways
kegg_pathway_info = fetch_google_sheet_data(sheet_ID, "KEGG_metadata_pairs in frozen", credentials_file) 
# HGNC gene group
gene_group = fetch_google_sheet_data(sheet_ID, "HGNC gene group", credentials_file)
 
src_info = fetch_google_sheet_data(sheet_ID, "sourceAbbv", credentials_file)
#pop_up_info = fetch_google_sheet_data(sheet_ID, "HGNC_Dec2024", credentials_file)
