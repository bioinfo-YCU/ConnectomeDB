import os
import jinja2
import sys
import pandas as pd
import numpy as np
import time
import base64
import re
from jinja2 import Environment, FileSystemLoader  

sys.path.append(os.path.abspath("src"))  
import fetchGSheet
from createTriplicateDT import gene_pair_trip

human_gene_pairTrip = gene_pair_trip.iloc[:, :11]
human_gene_pairTrip.rename(columns={human_gene_pairTrip.columns[0]: "Interaction ID"}, inplace=True)
# Paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TEMPLATE_DIR = "HTML"
TEMPLATE_FILE = "filteredDBTemplate.html"
OUTPUT_DIR = "data/filter/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === LOAD DATA ===
human_gene_pairTrip["BaseID"] = human_gene_pairTrip["Interaction ID"].str[:-1]
#human_gene_pairTrip = human_gene_pairTrip[human_gene_pairTrip["BaseID"] == "CDB00219"]

# === JINJA2 SETUP ===
env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
template = env.get_template(TEMPLATE_FILE)

# === GENERATE HTML FILES ===
for base_id, group_df in human_gene_pairTrip.groupby("BaseID"):
    group_df = group_df.drop(columns=["BaseID"])  # Drop it before rendering
    group_df= group_df.sort_values(by='Year', ascending=False)
    rendered_html = template.render(interaction_id=base_id, table=group_df)
    output_path = os.path.join(OUTPUT_DIR, f"{base_id}.html")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(rendered_html)

print(f"✅ Generated {len(human_gene_pairTrip['BaseID'].unique())} HTML files in '{OUTPUT_DIR}'")
