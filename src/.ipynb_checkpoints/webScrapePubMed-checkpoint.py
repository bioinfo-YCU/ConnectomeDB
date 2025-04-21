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
import fetchGSheet

# Read the API key from a file
with open("data/ncbi_api_key.txt", "r") as file:
    ncbi_api_key = file.read().strip()

# File to save the results
output_file = "data/pubmed_results.csv"

# Load your list of PMIDs
pmid_list = source

# Fetch HGNC gene symbols (you should have the `fetchGSheet.pop_up_info` dataframe ready)
def extract_hgnc_symbols(fetchGSheet):
    # Concatenate Approved, Alias, and Previous symbols, then extract unique symbols
    hgnc_symbols = pd.concat([
        fetchGSheet['Approved symbol'],
        fetchGSheet['Alias symbol'],
        fetchGSheet['Previous symbol']
    ], axis=0).dropna().str.upper().unique()  # Remove NaNs and make uppercase for matching
     # Remove any empty strings from the list
    hgnc_symbols = [symbol for symbol in hgnc_symbols if symbol != ""]
    return set(hgnc_symbols)  # Return as a set for fast lookup
hgnc_symbols = extract_hgnc_symbols(fetchGSheet.pop_up_info)

# Official species names and their corresponding terms (scientific names)
species_dict = {
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
    "rat": "Rattus norvegicus",
    "rabbit": "Oryctolagus cuniculus",
    "monkey": "Macaca spp.",
    "dog": "Canis lupus familiaris",
    "pig": "Sus scrofa",
    "zebra fish": "Danio rerio",
    "chicken": "Gallus gallus",
    "horse": "Equus ferus caballus",
    "cat": "Felis catus",
    "sheep": "Ovis aries",
    "cow": "Bos taurus",
    "fruit fly": "Drosophila melanogaster",
    "c. elegans": "Caenorhabditis elegans",
}

def fetch_pubmed_data(pmid_list, hgnc_symbols):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    results = []

    # Load existing data if output file exists
    if os.path.exists(output_file):
        existing_data = pd.read_csv(output_file)
    else:
        existing_data = pd.DataFrame(columns=["PMID", "Title", "Abstract", "Journal", "Year", "Species"])

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

                # Initialize species as N/A
                species = "N/A"

                # Check if the word "patient" is detected in title or abstract (assume human)
                if "patient" in title.lower() or "patient" in abstract.lower():
                    species = "Homo sapiens"
                elif "human" in title.lower() or "human" in abstract.lower():
                    species = "Homo sapiens"
                else:
                    # Look for HGNC gene symbols in title or abstract (assume human if found)
                    for gene in hgnc_symbols:
                        if gene in title or gene in abstract:
                            species = "Homo sapiens"
                            break
                    else:
                        # Look for MeSH terms related to species
                        for mesh_heading in article.findall(".//MeshHeadingList/MeshHeading"):
                            descriptor_name = mesh_heading.findtext("DescriptorName")
                            if descriptor_name:
                                # Match official species names using the species_dict
                                for species_term, scientific_name in species_dict.items():
                                    if species_term in descriptor_name.lower():
                                        species = scientific_name
                                        break  # Stop after finding the first match

                # Append the result
                results.append({
                    "PMID": article.findtext(".//MedlineCitation/PMID"),
                    "Title": title,
                    "Abstract": abstract,
                    "Journal": journal,
                    "Year": year,
                    "Species": species
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

# Fetch PubMed data with your list of PMIDs, output file path, and NCBI API key
fetch_pubmed_data(pmid_list, hgnc_symbols)