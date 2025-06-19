from pathlib import Path
import sys
import pandas as pd
import re
import html
import json
from jinja2 import Environment, FileSystemLoader

# === Import from createDataTable.py
sys.path.append("src")
from createDataTable import human_gene_pair  # must return a DataFrame

# === Paths (relative to project/)
output_json = Path("JSON/human_gene_pair.json")
qmd_template = Path("database/qmd_template/human_template.qmd")
qmd_output = Path("database/human.qmd")
template_dir = "HTML"  # where datatable_template.html lives
template_name = "datatable_template.html"

# === Create output directories if needed
output_json.parent.mkdir(parents=True, exist_ok=True)

# === Clean column names and Save DataFrame to JSON
human_gene_pair = human_gene_pair.iloc[:, :-48]
def clean_column_names_and_generate_metadata(df):
    # Extract visible column names from HTML
    def visible_text(html_string):
        return re.sub(r'<[^>]*>', '', html.unescape(html_string)).strip()

    raw_columns = df.columns.tolist()
    visible_columns = [visible_text(col) for col in raw_columns]

    # Rename DataFrame columns to their visible text
    df_cleaned = df.copy()
    df_cleaned.columns = visible_columns

    # Create DataTables metadata with original HTML as `title`
    column_metadata = [
        {
            "data": visible,
            "title": html_col
        }
        for visible, html_col in zip(visible_columns, raw_columns)
    ]

    return df_cleaned, column_metadata
# === Generate DataTables column definitions
df_cleaned, columns_metadata = clean_column_names_and_generate_metadata(human_gene_pair)

# Now use df_cleaned to export to JSON
df_cleaned.to_json("JSON/human_gene_pair.json", orient="records")

# And use columns_metadata as your DataTables columns
columns_json = json.dumps(columns_metadata, indent=2)


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

if "{{ table_block }}" not in contents:
    raise ValueError("Placeholder {{ table_block }} not found in human.qmd")

with open(qmd_output, "w") as f:
    f.write(contents.replace("{{ table_block }}", rendered_html))

print(f"✔️ Updated {qmd_output}")
print(f"✔️ Saved JSON to {output_json}")
