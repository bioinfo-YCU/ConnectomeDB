from pathlib import Path
import sys
import pandas as pd
import re
import html
import json
from jinja2 import Environment, FileSystemLoader

# === Import from createDataTable.py
sys.path.append("src")
import createDataTable_perSpecies

# === Paths (relative to project/)
species = "mouse"
# additional placeholder for search
species_addl_search = "Epha3, 5430401F13Rik, no human ortholog"
output_json = Path(f"JSON/{species}_gene_pair.json")  
qmd_template = Path(f"database/qmd_template/{species}Orth_template.qmd") 
if species == "mouse":
    qmd_output = Path(f"database/{species}Orth.qmd")  
else:
    qmd_output = Path(f"database/other/{species}Orth.qmd")
template_dir = "HTML"
template_name = "datatableOrth_template.html"

# === Create output directories if needed
output_json.parent.mkdir(parents=True, exist_ok=True)

# === Clean column names and Save DataFrame to JSON
def clean_column_names_and_generate_metadata(df):
    def visible_text(html_string):
        return re.sub(r'<[^>]*>', '', html.unescape(html_string)).strip()
    
    raw_columns = df.columns.tolist()
    visible_columns = [visible_text(col) for col in raw_columns]
    
    df_cleaned = df.copy()
    df_cleaned.columns = visible_columns
    
    column_metadata = [
        {
            "data": visible,
            "title": html_col
        }
        for visible, html_col in zip(visible_columns, raw_columns)
    ]
    return df_cleaned, column_metadata

# === Generate DataTables column definitions
# Use getattr to dynamically access the attribute
gene_pair = getattr(createDataTable_perSpecies, f"{species}_gene_pair1")  
df_cleaned, columns_metadata = clean_column_names_and_generate_metadata(gene_pair)

# Export to JSON with f-string
df_cleaned.to_json(f"JSON/{species}_gene_pair.json", orient="records") 
columns_json = json.dumps(columns_metadata, indent=2)

# === Jinja2 render
env = Environment(loader=FileSystemLoader(template_dir))
template = env.get_template(template_name)
rendered_html = template.render(
    columns_json=columns_json,
    species=species, 
    species_addl_search = species_addl_search,
    table_id=f"{species}-table", 
    json_path=f"../JSON/{species}_gene_pair.json" 
)

# === Inject into species.qmd at placeholder
with open(qmd_template, "r") as f:
    contents = f.read()

if "{{ table_block }}" not in contents:
    raise ValueError(f"Placeholder {{{{ table_block }}}} not found in {qmd_template}")  # ✅ Better error message

with open(qmd_output, "w") as f:
    f.write(contents.replace("{{ table_block }}", rendered_html))

print(f"✔️ Updated {qmd_output}")
print(f"✔️ Saved JSON to {output_json}")