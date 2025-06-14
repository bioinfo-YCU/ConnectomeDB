import os
from bs4 import BeautifulSoup
from sentence_transformers import SentenceTransformer
import faiss
import numpy as np
import json

# Define relative paths
ROOT_DIR = os.path.dirname(os.path.dirname(__file__))  # goes up one level from /src
DATA_DIR = os.path.join(ROOT_DIR, "_site")
OUT_DIR = os.path.join(ROOT_DIR, "data")
INDEX_FILE = os.path.join(OUT_DIR, "faiss_index.bin")
META_FILE = os.path.join(OUT_DIR, "faiss_metadata.json")

model = SentenceTransformer('all-MiniLM-L6-v2')

def extract_text_from_html(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        soup = BeautifulSoup(f, 'html.parser')
        texts = soup.stripped_strings
        return " ".join(texts)

def build_index():
    texts = []
    metadata = []

    for root, _, files in os.walk(DATA_DIR):  # recursive through _site/
        for filename in files:
            if filename.endswith(".html"):
                full_path = os.path.join(root, filename)
                text = extract_text_from_html(full_path)
                texts.append(text)
                metadata.append({
                    "filename": os.path.relpath(full_path, DATA_DIR),
                    "text": text[:200]
                })

    embeddings = model.encode(texts, convert_to_numpy=True)
    dim = embeddings.shape[1]
    index = faiss.IndexFlatL2(dim)
    index.add(embeddings)

    os.makedirs(OUT_DIR, exist_ok=True)
    faiss.write_index(index, INDEX_FILE)

    with open(META_FILE, "w") as f:
        json.dump(metadata, f)
        print(f"Index built with {len(metadata)} documents.")

if __name__ == "__main__":
    build_index()
    
