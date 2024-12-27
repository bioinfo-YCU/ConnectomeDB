## Function to scrape data from Pubmed for Title, Abstract, Journal, and Year
### IMPORTANT: TURN OFF VPN and make sure you have the data directory (from Sakura)

import sys
import requests
import pandas as pd
import time
import os
import xml.etree.ElementTree as ET

sys.path.append(os.path.abspath("src"))  
from createDataTable import source

# Read the API key from a file
with open("data/ncbi_api_key.txt", "r") as file:
    ncbi_api_key = file.read().strip()

# File to save the results
output_file = "data/pubmed_results.csv"

# Load your list of PMIDs
pmid_list = source

if os.path.exists(output_file):
    existing_data = pd.read_csv(output_file)
    existing_pmids = set(existing_data["PMID"].astype(str))  # Ensure PMIDs are strings
else:
    existing_pmids = set()
    
# Limit PMIDs to the intersection with an existing file
#pmid_list = list(set(pmid_list).intersection(existing_pmids))

def fetch_pubmed_data(pmid_list):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    results = []

    # Load existing data if output file exists
    if os.path.exists(output_file):
        existing_data = pd.read_csv(output_file)
    else:
        existing_data = pd.DataFrame(columns=["PMID", "Title", "Abstract", "Journal", "Year"])

    # Split PMIDs into batches
    batch_size = 50
    pmid_batches = [pmid_list[i:i + batch_size] for i in range(0, len(pmid_list), batch_size)]

    # Iterate over the batches
    for batch in pmid_batches:
        params = {
            "db": "pubmed",
            "id": ",".join(batch),  # Join PMIDs as comma-separated
            "retmode": "xml",
            "api_key": ncbi_api_key
        }

        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()

            # Parse the XML response
            root = ET.fromstring(response.text)
            for article in root.findall(".//PubmedArticle"):
                # Extract Title and Abstract
                title = article.findtext(".//ArticleTitle", default="N/A")
                abstract = article.findtext(".//AbstractText", default="No abstract available")

                # Extract Journal Title
                journal_tag = article.find(".//Journal/Title")
                journal = journal_tag.text.strip() if journal_tag is not None and journal_tag.text else "N/A"

                # Extract Publication Year
                pub_date = article.find(".//PubDate")
                if pub_date is not None:
                    year_tag = pub_date.find("Year")
                    year = year_tag.text if year_tag is not None else "N/A"

                    # Fallback to MedlineDate if Year is missing
                    if year == "N/A":
                        medline_date_tag = pub_date.find("MedlineDate")
                        year = medline_date_tag.text.split()[0] if medline_date_tag is not None else "N/A"
                else:
                    year = "N/A"  # PubDate is completely missing

                # Append the result
                results.append({
                    "PMID": article.findtext(".//MedlineCitation/PMID"),
                    "Title": title,
                    "Abstract": abstract,
                    "Journal": journal,
                    "Year": year
                })

        except Exception as e:
            print(f"Error fetching batch {batch}: {e}")
            # Optionally save the response for debugging
            with open(f"error_batch_{batch[0]}_{batch[-1]}.xml", "w") as f:
                f.write(response.text)

        # Rate limiting to avoid API overload
        time.sleep(1)  # Increase delay for better API compliance

    # Save results
    new_data = pd.DataFrame(results)
    if not new_data.empty:
        # Merge existing and new data, updating missing values
        updated_data = pd.concat([existing_data, new_data])
        
        # Ensure all PMIDs are strings
        updated_data["PMID"] = updated_data["PMID"].astype(str)
        
        # Drop rows with missing PMIDs
        updated_data = updated_data.dropna(subset=["PMID"])
        
        # Ensure rows are ordered and remove duplicates
        updated_data = (
            updated_data.sort_values(by="PMID")  # Ensure rows are ordered
            .drop_duplicates(subset="PMID", keep="last")  # Keep the latest data
        )
        updated_data["Journal"] = updated_data["Journal"].str.split(" (", n=1, expand=False, regex=False).str[0]
        updated_data.to_csv(output_file, index=False)
    else:
        print("No new data fetched.")

    return results

# Fetch data for the intersected PMIDs
fetch_pubmed_data(pmid_list)