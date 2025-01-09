## Function to run LLM on all abstract for each LR pair to come up with biological relevance (keywords)
import pandas as pd
import re
import ijson
from transformers import AutoTokenizer
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import ENGLISH_STOP_WORDS

# Set output file name
output_file = "data/llm_results.csv"
top_n = 8

# Load tokenizer
tokenizer = AutoTokenizer.from_pretrained("facebook/bart-large-cnn")

# Define your custom stop/include words
custom_stop_words = [
    "make", "significant", "activity", "made", "makes", "significantly",
    "activities", "activity", "attenuated", "induced", "enhanced", "attenuates", 
    "induces", "enhances", "available", "abstract", "rarely", "result", "results",
    "produced", "produce", "important", "prominent", "role", "11", "12", "10",
    "18", "19", "",
]
custom_include_keywords = [
    "cancer", "tumor", "tumour", "signaling", "signalling", "cell communication", "protein-coding",
    "protein coding", "non-protein coding"
    "lncrna", "lnc-rna", "differentiation", "immune response", "non-inflammatory", "inflammatory",
    "hypoxia", "TME", "microenvironment", "cell-cell", "cell-to-cell", "dendritic differentiation",
    "gene regulation", "heterogeneity", "growth factor"
]

# Combine default English stop words and custom stop words
combined_stop_words = list(ENGLISH_STOP_WORDS.union(custom_stop_words))

# Initialize CountVectorizer with the updated stop words
vectorizer = CountVectorizer(max_features=top_n, stop_words=combined_stop_words)


def count_tokens(text):
    """
    Count the number of tokens in a given text using the tokenizer.
    """
    return len(tokenizer.encode(text, truncation=False))


def chunk_abstracts(abstracts, max_tokens=16000):
    """
    Split the abstracts into chunks based on the maximum token limit.
    """
    chunks = []
    current_chunk = []
    current_tokens = 0

    for abstract in abstracts:
        abstract_tokens = count_tokens(abstract)
        if current_tokens + abstract_tokens > max_tokens:
            chunks.append(" ".join(current_chunk))
            current_chunk = [abstract]
            current_tokens = abstract_tokens
        else:
            current_chunk.append(abstract)
            current_tokens += abstract_tokens

    if current_chunk:
        chunks.append(" ".join(current_chunk))

    return chunks

def extract_keywords(text, top_n=top_n):
    """
    Extract top N keywords (including multi-word keywords) from a text using CountVectorizer.
    Exclude combined stop words, ensure inclusion of custom keywords, and filter out numeric keywords.
    Handle empty vocabulary errors.
    """
    if not text.strip():
        return ""  # Return an empty string if the text is empty

    try:
        # Remove purely numeric words from the text
        filtered_text = " ".join(word for word in text.split() if not re.match(r"^\d+$", word))

        # Configure CountVectorizer to include n-grams (e.g., unigrams and bigrams)
        vectorizer = CountVectorizer(
            max_features=top_n,
            stop_words=combined_stop_words,
            ngram_range=(1, 5)  # This includes unigrams (single words), bigrams (two words) till 5 words
        )

        # Fit the vectorizer to the filtered text and extract keywords
        X = vectorizer.fit_transform([filtered_text])
        keywords = set(vectorizer.get_feature_names_out())

        # Check which custom keywords (including multi-word) are present in the original text
        found_custom_keywords = [keyword for keyword in custom_include_keywords if keyword in text.lower()]

        # Add the found custom keywords to the set of extracted keywords
        keywords.update(found_custom_keywords)

        return ", ".join(keywords)
    except ValueError as e:
        # Handle empty vocabulary error
        if "empty vocabulary" in str(e):
            return ""  # Return an empty string if no valid tokens remain
        raise  # Re-raise other unexpected errors


def process_human_lr_pair(human_lr_pair, abstracts):
    """
    Process each ligand-receptor pair and extract relevant biological keywords.
    """
    # Chunk abstracts into smaller parts
    chunks = chunk_abstracts(abstracts)
    all_keywords = []

    for chunk in chunks:
        # Extract keywords from each chunk
        keywords = extract_keywords(chunk)
        all_keywords.append(keywords)

    # Combine all the extracted keywords from each chunk
    return {"Human LR Pair": human_lr_pair, "Relevance Keywords": ", ".join(all_keywords)}


# Stream JSON and process incrementally
with open("data_for_llm.json", "r") as f:
    parser = ijson.items(f, "item")
    results = []

    for entry in parser:
        human_lr_pair = entry["Human LR Pair"]
        abstracts = entry["Abstracts"]
        result = process_human_lr_pair(human_lr_pair, abstracts)

        # Skip entries where keywords are empty or only numerical
        if not result["Relevance Keywords"] or result["Relevance Keywords"].replace(",", "").isdigit():
            continue

        results.append(result)

# Write the final results to CSV (overwrite mode)
pd.DataFrame(results).to_csv(output_file, mode="w", header=True, index=False)

print(f"Results written to {output_file}")
