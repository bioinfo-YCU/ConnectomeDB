## Function to create json files to be used for datatable and prepare qmd (Mouse and other orth)
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



def process_species_table(species, species_addl_search, createDataTable_perSpecies):
    # === Paths (relative to project/)
    output_json = Path(f"JSON/{species}_gene_pair.json")  
    qmd_template = Path(f"database/qmd_template/{species}Orth_template.qmd") 
    qmd_output = Path(f"database/{species}Orth.qmd") if species == "mouse" else Path(f"database/other/{species}Orth.qmd")
    template_dir = "HTML"
    template_name = "datatableOrth_template.html"

    # === Create output directories if needed
    output_json.parent.mkdir(parents=True, exist_ok=True)

    # === Clean column names and generate metadata
    def clean_column_names_and_generate_metadata(df):
        def visible_text(html_string):
            return re.sub(r'<[^>]*>', '', html.unescape(html_string)).strip()

        raw_columns = df.columns.tolist()
        visible_columns = [visible_text(col) for col in raw_columns]

        df_cleaned = df.copy()
        df_cleaned.columns = visible_columns

        column_metadata = [
            {"data": visible, "title": html_col}
            for visible, html_col in zip(visible_columns, raw_columns)
        ]
        return df_cleaned, column_metadata

    # === Get the gene_pair DataFrame
    try:
        gene_pair = getattr(createDataTable_perSpecies, f"{species}_gene_pair1")
    except AttributeError:
        raise ValueError(f"Missing attribute for {species}_gene_pair1 in createDataTable_perSpecies")

    # === Clean and export to JSON
    df_cleaned, columns_metadata = clean_column_names_and_generate_metadata(gene_pair)
    df_cleaned.to_json(output_json, orient="records")
    columns_json = json.dumps(columns_metadata, indent=2)

    # === Jinja2 render
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template(template_name)
    rendered_html = template.render(
        columns_json=columns_json,
        species=species, 
        species_addl_search=species_addl_search,
        table_id=f"{species}-table", 
        json_path=f"../JSON/{species}_gene_pair.json" if species == "mouse" else Path(f"../../JSON/{species}_gene_pair.json")
    )

    # === Replace placeholder in .qmd
    if not qmd_template.exists():
        raise FileNotFoundError(f"Missing template file: {qmd_template}")
    
    contents = qmd_template.read_text()
    # Add today's date to any title that ends with </span>
    today = datetime.now().strftime("%Y-%m-%d")
    contents = re.sub(
        r'(title: ".*?) </span>"',
        rf'\1 ({today}) </span>"',
        contents
    )
    if "{{ table_block }}" not in contents:
        raise ValueError(f"Placeholder '{{{{ table_block }}}}' not found in {qmd_template}")

    qmd_output.parent.mkdir(parents=True, exist_ok=True)
    qmd_output.write_text(contents.replace("{{ table_block }}", rendered_html))

    print(f"✔️ Updated {qmd_output}")
    print(f"✔️ Saved JSON to {output_json}")


# === Run for each species
species_list = ["mouse", "rat", "zebrafish", "frog", "chicken", "macaque", "pig", "dog", "cow", "chimp", "horse", "marmoset", "sheep"]
species_addl_search_dict = {
    "mouse": "Epha3, 5430401F13Rik",
    "rat": "Cck, Vav1",
    "zebrafish": "adam10b, zp3a.2",
    "frog": "ccl5, rantes",
    "chicken": "BDKRB2, LOC395551",
    "macaque": "EPHA3, CADM1",
    "pig": "ACKR4, CCL21",
    "dog": "LOC480600, CTLA4",
    "chimp": "CEACAM5, NRXN3", 
    "cow": "LOC508666, CCR10",
    "horse": "LOC100065387, OSCAR",
    "marmoset": "LOC100065387, OSCAR",
    "sheep": "ENSOARG00020015296, IL10RA",
}

for sp in species_list:
    process_species_table(sp, species_addl_search_dict[sp], createDataTable_perSpecies)
