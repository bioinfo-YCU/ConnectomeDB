## Function to curate a redirect table based on the info from all 14 tables including humans
from pathlib import Path
import sys
import pandas as pd
import re
import html
import json
from jinja2 import Environment, FileSystemLoader
from datetime import datetime

# === Import from createDataTable.py
sys.path.append("src")
import createDataTable_perSpecies
import createDataTable

# today's date for version control
today = datetime.now().strftime("%Y%m%d")  # e.g., '20250723'

def build_combined_species_table(species_list, source_obj):
    """
    Builds a combined DataFrame for all species in species_list
    by extracting Interaction ID, Card URL, LR Pair, and redirect URLs.
    
    Parameters:
        species_list (list): list of species names (strings)
        source_obj (object): object containing {species}_gene_pair1 DataFrames
    
    Returns:
        pd.DataFrame: combined table across all species
    """
    all_species_dfs = []

    for species in species_list:
        # Get the species-specific DataFrame
        if species == "human":
            Species_df = getattr(createDataTable, f"{species}_gene_pair")
        else:
            Species_df = getattr(createDataTable_perSpecies, f"{species}_gene_pair1")

        
        # Grab only the first two columns
        Species_df = Species_df.iloc[:, :2].copy()
        
        # Identify the actual column names
        interaction_id_col = [col for col in Species_df.columns if "Interaction ID" in col][0]
        LRpair_col = [col for col in Species_df.columns if "LR Pair" in col][0]

        # Extract Card URL and Interaction ID from first column (HTML)
        Species_df["Card URL"] = Species_df[interaction_id_col].str.extract(r'href="([^"]+)"')
        Species_df["Interaction ID"] = Species_df[interaction_id_col].str.extract(r'>([^<]+)<')

        # Rename 2nd column to "LR Pair"
        Species_df = Species_df.rename(columns={LRpair_col: "LR Pair"})

        # Create "Species redirect URL"
        Species_df["Species redirect URL"] = Species_df["LR Pair"].apply(
            lambda x: f"https://connectomedb.org/cards/{species}/{x.replace(' ', '-')}.html"
        )

        # Create "Direct URL"
        Species_df["Direct URL"] = Species_df["Interaction ID"].apply(
            lambda x: f"https://connectomedb.org/cards/human/{x}.html"
        )

        # Drop original HTML column
        Species_df = Species_df.drop(columns=[interaction_id_col])

        # Reorder columns
        first_cols = ["Interaction ID", "Card URL", "LR Pair", "Species redirect URL", "Direct URL"]
        Species_df = Species_df[first_cols + [c for c in Species_df.columns if c not in first_cols]]

        # Keep track of which species it came from
        Species_df["Species"] = species

        all_species_dfs.append(Species_df)

    # Combine all species into one long DataFrame
    combined_df = pd.concat(all_species_dfs, ignore_index=True)
    return combined_df


# Usage example:
species_list = ["human","mouse", "rat", "zebrafish", "frog", "chicken", "macaque", "pig", "dog", "cow", "chimp", "horse", "marmoset", "sheep"]

combined_df = build_combined_species_table(species_list, createDataTable_perSpecies)

# Save to CSV
combined_df.to_csv("data/redirect_table.csv", index=False)
combined_df.to_csv(f"data/redirect_table_{today}.csv", index=False) # version control
