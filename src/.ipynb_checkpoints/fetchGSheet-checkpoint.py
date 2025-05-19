## Function to fetch data from a Google Sheet tab

import gspread
from google.oauth2.service_account import Credentials
import pandas as pd

sheet_ID = "1XP5wBDN_orSlE8RLb2TxSclIopVO1mb_1S3ENf2qYSw" #"15FfI7cVpJmAcytTBmhVE2Z77tgVQYfEZUEe700_Wzlg"
credentials_file = 'data/connectomedb2025-a9acdf562a84.json'

def fetch_google_sheet_data(sheet_ID, tab_name, credentials_file):
    """
    Fetch data from a specific tab of a Google Sheet and return it as a pandas DataFrame.
    
    Parameters:
        sheet_ID (str): The ID of the Google Sheet.
        tab_name (str): The name of the tab/worksheet to fetch.
        credentials_file (str): Path to the service account JSON credentials file.
    
    Returns:
        pd.DataFrame: Data from the specified Google Sheet tab.
    """
    # Define the scopes and authenticate
    SCOPES = ['https://www.googleapis.com/auth/spreadsheets', 'https://www.googleapis.com/auth/drive']
    credentials = Credentials.from_service_account_file(credentials_file, scopes=SCOPES)
    client = gspread.authorize(credentials)

    # Open the Google Sheet and get the specified tab
    try:
        sheet = client.open_by_key(sheet_ID)
        worksheet = sheet.worksheet(tab_name)
    except gspread.exceptions.WorksheetNotFound:
        raise ValueError(f"Tab '{tab_name}' not found in Google Sheet with ID '{sheet_ID}'.")

    # Get all values from the tab and convert to DataFrame
    data = worksheet.get_all_values()
    df = pd.DataFrame(data)
    df.columns = df.iloc[0]  # Set the first row as the header
    df = df[1:].reset_index(drop=True)  # Remove the header row from the data
    return df

# Fetching data from Google Sheets
gene_pair = fetch_google_sheet_data(sheet_ID, "FROZEN LIST HUMAN", credentials_file)
loc_info = fetch_google_sheet_data(sheet_ID, "proteome_HPA", credentials_file) 
# Integrate 
ligand_loc = fetch_google_sheet_data(sheet_ID, "Ligand_location_HUMAN", credentials_file) 
receptor_loc = fetch_google_sheet_data(sheet_ID, "Receptor_location_HUMAN", credentials_file) 
kegg_pathway_info = fetch_google_sheet_data(sheet_ID, "KEGG_metadata_pairs in frozen", credentials_file) 
# previously based on localization # another option is proteome_HPA
src_info = fetch_google_sheet_data(sheet_ID, "sourceAbbv", credentials_file)
#pop_up_info = fetch_google_sheet_data(sheet_ID, "HGNC_Dec2024", credentials_file)
