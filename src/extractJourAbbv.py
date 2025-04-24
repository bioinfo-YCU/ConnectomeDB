## Function to scrape data from Pubmed for Title, Abstract, Journal, and Year
### IMPORTANT: TURN OFF VPN and make sure you have the data directory (from Sakura)

import sys
import requests
import pandas as pd
import time
import os
import xml.etree.ElementTree as ET

sys.path.append(os.path.abspath("src"))  
import re
from fuzzywuzzy import process, fuzz
import csv

# Output CSV file
output_file = "data/journal_abbv.csv"

pubmed_df = pd.read_csv("data/pubmed_results.csv")
journal_names = pubmed_df["Journal"].unique().tolist()

manual_abbr_dict = {
    "The Journal of biological chemistry": "J. Biol. Chem.",
    "Journal of immunology": "J. Immunol.",
    "Acta physiologica": "Acta Physiol. (Oxf.)",
    "The Biochemical journal": "Biochem. J.",
    "Hepatology": "Hepatology",
    "Chemical reviews": "Chem. Rev.",
    "Molecular endocrinology": "Mol. Endocrinol.",
    "Journal of molecular biology": "J. Mol. Biol.",
    "The Journal of experimental medicine": "J. Exp. Med.",
    "Growth factors": "Growth Factors",
    "Development": "Development",
    "Structure": "Structure",
    "Journal of neurochemistry": "J. Neurochem.",
    "Cancer research": "Cancer Res.",
    "Advanced science": "Adv. Sci.",
    "Arthritis & rheumatology": "Arthritis Rheumatol.",
    "Clinical immunology": "Clin. Immunol.",
    "Hypertension": "Hypertension",
    "Neoplasia": "Neoplasia",
    "The Journal of general physiology": "J. Gen. Physiol.",
    "Acta biochimica Polonica": "Acta Biochim. Pol.",
    "Cell cycle": "Cell Cycle",
    "Human reproduction": "Hum. Reprod.",
    "American journal of reproductive immunology": "Am. J. Reprod. Immunol.",
    "Methods": "Methods",
    "Lung cancer": "Lung Cancer",
    "Gut": "Gut",
    "Archives of surgery": "Arch. Surg.",
    "Lancet": "Lancet"
}


medline_file= "data/J_Medline.txt"
def load_journal_info(filename):
    journals = []
    with open(filename, 'r', encoding='utf-8') as file:
        content = file.read()
        # Split the content by the separator that divides journal entries
        entries = content.split('--------------------------------------------------------')
        
        for entry in entries:
            title_match = re.search(r"JournalTitle:\s*(.*?)\n", entry)
            abbr_match = re.search(r"MedAbbr:\s*(.*?)\n", entry)
            
            if title_match and abbr_match:
                title = title_match.group(1).strip()
                abbr = abbr_match.group(1).strip()
                journals.append((title, abbr))
    
    return journals


def get_abbreviations(journal_names, journal_dict, score_threshold=98):
    journal_keys = [title for title, _ in journal_dict]
    results = []

    for name in journal_names:
        cleaned_name = name.strip().lower()

        # 1. Check for exact match
        abbr = None
        for title, abbreviation in journal_dict:
            if cleaned_name == title.strip().lower():
                abbr = abbreviation
                results.append((name, abbr, "Exact"))
                break

        # 2. Fuzzy match
        if not abbr:
            matches = process.extract(cleaned_name, journal_keys, scorer=fuzz.ratio, limit=5)
            best_match = matches[0] if matches else None
            if best_match:
                matched_name, score = best_match
                if score >= score_threshold:
                    abbr = journal_dict[journal_keys.index(matched_name)][1]
                    results.append((name, abbr, f"Fuzzy (score: {score})"))

        # 3. Manual backup
        if not abbr:
            manual_match = manual_abbr_dict.get(name.strip())
            if manual_match:
                results.append((name, manual_match, "Manual Backup"))
            else:
                results.append((name, "Not Found", "No Match"))

    return results

def save_to_csv(results, output_file):
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Journal Name', 'Abbreviation', 'Match Type']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for result in results:
            writer.writerow({'Journal Name': result[0], 'Abbreviation': result[1], 'Match Type': result[2]})
# Load journal information from the file
journal_list = load_journal_info(medline_file)

#Get the matches
results = get_abbreviations(journal_names, journal_list)

# Save results to CSV
save_to_csv(results, output_file)
