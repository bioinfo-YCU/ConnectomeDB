## Function to create json files to be used for datatable prepare qmd (Human)
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
from createDataTable import human_gene_pair  # must return a DataFrame

# === DEBUG: Print column information
print("Original DataFrame columns:")
for i, col in enumerate(human_gene_pair.columns):
    print(f"  {i}: {repr(col)}")

print(f"\nDataFrame shape: {human_gene_pair.shape}")

# === Paths (relative to project/)
output_json = Path("JSON/human_gene_pair.json")
qmd_template = Path("database/qmd_template/human_template.qmd")
qmd_output = Path("database/human.qmd")
template_dir = "HTML"  # where datatable_template.html lives
template_name = "datatable_template.html"

# === Create output directories if needed
output_json.parent.mkdir(parents=True, exist_ok=True)

# === Clean column names and Save DataFrame to JSON
def clean_column_names_and_generate_metadata(df):
    # Extract visible column names from HTML
    def visible_text(html_string):
        cleaned = re.sub(r'<[^>]*>', '', html.unescape(html_string)).strip()
        return cleaned
    
    # Create safe column names for JSON keys (no spaces, periods, special chars)
    def make_safe_key(text):
        # Replace spaces and periods with underscores, remove other special chars
        safe = re.sub(r'[^a-zA-Z0-9_]', '_', text)
        # Remove multiple underscores
        safe = re.sub(r'_+', '_', safe)
        # Remove leading/trailing underscores
        safe = safe.strip('_')
        return safe

    raw_columns = df.columns.tolist()
    visible_columns = [visible_text(col) for col in raw_columns]
    safe_columns = [make_safe_key(col) for col in visible_columns]

    print(f"\nColumn name mapping:")
    for i, (raw, visible, safe) in enumerate(zip(raw_columns, visible_columns, safe_columns)):
        print(f"  {i}: HTML='{raw}' -> Visible='{visible}' -> Safe='{safe}'")

    # Create a copy of the DataFrame with safe column names
    df_cleaned = df.copy()
    df_cleaned.columns = safe_columns

    print(f"\nDataFrame after renaming columns:")
    for i, col in enumerate(df_cleaned.columns):
        print(f"  {i}: {repr(col)}")

    # Create DataTables metadata using safe column names as 'data' and HTML as 'title'
    column_metadata = [
        {
            "data": safe_col,     # Use safe column name for data mapping
            "title": raw_col      # Use original HTML for column headers
        }
        for safe_col, raw_col in zip(safe_columns, raw_columns)
    ]

    print(f"\nColumn metadata:")
    for i, meta in enumerate(column_metadata):
        print(f"  {i}: data='{meta['data']}', title='{meta['title']}'")

    return df_cleaned, column_metadata

# === Generate DataTables column definitions
df_cleaned, columns_metadata = clean_column_names_and_generate_metadata(human_gene_pair)

# Export the cleaned DataFrame to JSON
df_cleaned.to_json("JSON/human_gene_pair.json", orient="records")

# Check what's actually in the JSON
with open("JSON/human_gene_pair.json", "r") as f:
    json_data = json.load(f)
    if json_data:
        print(f"\nFirst record keys in JSON:")
        for key in json_data[0].keys():
            print(f"  {repr(key)}")

# Generate columns JSON for DataTables
columns_json = json.dumps(columns_metadata, indent=2)

print(f"\nColumns JSON (first 3 entries):")
print(json.dumps(columns_metadata[:3], indent=2))

# === Jinja2 render
env = Environment(loader=FileSystemLoader(template_dir))
template = env.get_template(template_name)

rendered_html = template.render(
    columns_json=columns_json,
    table_id="human-table",
    json_path="../JSON/human_gene_pair.json"  # relative to `database/human.qmd`
)

# === Inject into human.qmd at placeholder
with open(qmd_template, "r") as f:
    contents = f.read()
    # Add today's date to any title that ends with </span>
    today = datetime.now().strftime("%Y-%m-%d")
    contents = re.sub(
        r'(title: ".*?) </span>"',
        rf'\1 ({today}) </span>"',
        contents
    )

if "{{ table_block }}" not in contents:
    raise ValueError("Placeholder {{ table_block }} not found in human.qmd")

with open(qmd_output, "w") as f:
    f.write(contents.replace("{{ table_block }}", rendered_html))

print(f"✅ Updated {qmd_output}")
print(f"✅ Saved JSON to {output_json}")