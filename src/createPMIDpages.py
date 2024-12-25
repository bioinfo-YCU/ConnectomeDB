## Function to create an evidence page per LR pair with each tab per PMID

import sys
import os
import pandas as pd

# Paths
TEMPLATE_PATH = "HTML/pmidTemplate.html"
OUTPUT_DIR = "data/pubmed/"

# Add the src directory to the path for importing modules
sys.path.append(os.path.abspath("src"))
from createDataTable import source

# Load PubMed data
pubmed_data = pd.read_csv("data/pubmed_results.csv")
pubmed_data["Year"] = pubmed_data["Year"].astype(str).str.replace(".0", 
                                                                  "", 
                                                                  regex=False).astype(int)

pubmed_data["PMID"] = pubmed_data["PMID"].astype(str)
# Replace spaces in "LR Pair" with a placeholder
gene_pair00["LR Pair"] = gene_pair00["LR Pair"].str.replace(" ", "——")

pubmed_data = pubmed_data.reset_index(drop=True)  # Remove the index

def load_template(template_path):
    """
    Load the HTML template from a file.

    Parameters:
        template_path (str): Path to the HTML template file.

    Returns:
        str: The contents of the template file as a string.
    """
    with open(template_path, "r") as file:
        return file.read()


def create_detailed_pages_with_tabs(df, gene_column, pmid_column, pubmed_data, template):
    """
    Generate detailed HTML pages with tabs and PubMed links for each gene pair.

    Parameters:
        df (pd.DataFrame): DataFrame containing gene pair information.
        gene_column (str): Name of the column containing gene pair names.
        pmid_column (str): Name of the column containing PMIDs.
        pubmed_data (pd.DataFrame): DataFrame containing PubMed details.
        template (str): HTML template to use for generating the pages.
    """
    for idx, row in df.iterrows():
        gene_name = row[gene_column]
        pmids = row[pmid_column]

        # Ensure PMIDs are properly processed
        sources = [pmid.strip() for pmid in pmids.split(',') if pmid.strip()]

        if sources:
            tab_headers = []
            tab_contents = []

            for i, pmid in enumerate(sources):
                pubmed_row = pubmed_data[pubmed_data["PMID"] == pmid]

                if not pubmed_row.empty:
                    title = pubmed_row["Title"].values[0]
                    abstract = pubmed_row["Abstract"].values[0]
                    journal = pubmed_row["Journal"].values[0]
                    year = pubmed_row["Year"].values[0]
                else:
                    title = "No Title Found"
                    abstract = "No Abstract Found"
                    journal = "Journal Unknown"
                    year = "Year Unknown"

                active_class = "active" if i == 0 else ""
                tab_headers.append(f'<button class="tablinks {active_class}" onclick="openTab(event, \'tab{pmid}\')">{pmid}</button>')
                tab_contents.append(f"""
                <div id="tab{pmid}" class="tabcontent {active_class}">
                    <h2>{title}</h2>
                    <p><strong>{journal}, {year}; <a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">For more details, see PubMed</a></strong></p>
                    <p>{abstract}</p>
                    
                </div>
                """)

            # Fill the template with dynamic content
            page_content = template.replace("{{GENE_NAME}}", gene_name)
            page_content = page_content.replace("{{TAB_HEADERS}}", "".join(tab_headers))
            page_content = page_content.replace("{{TAB_CONTENTS}}", "".join(tab_contents))

            # Save the HTML file
            os.makedirs(OUTPUT_DIR, exist_ok=True)
            output_path = os.path.join(OUTPUT_DIR, f"{gene_name}_pmid_details.html")
            with open(output_path, "w") as file:
                file.write(page_content)


if __name__ == "__main__":
    # Load the HTML template
    template = load_template(TEMPLATE_PATH)

    # Generate pages
    create_detailed_pages_with_tabs(
        df=gene_pair00,
        gene_column="LR Pair",
        pmid_column="PMID support",
        pubmed_data=pubmed_data,
        template=template,
    )
