from http.server import BaseHTTPRequestHandler, HTTPServer
import os
import json
from sentence_transformers import SentenceTransformer
import faiss
import numpy as np
from ollama import Client

# Paths
ROOT_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(ROOT_DIR, "data")
INDEX_FILE = os.path.join(DATA_DIR, "faiss_index.bin")
META_FILE = os.path.join(DATA_DIR, "faiss_metadata.json")


embed_model = SentenceTransformer('all-MiniLM-L6-v2')
index = faiss.read_index(INDEX_FILE)
with open(META_FILE, "r") as f:
    metadata = json.load(f)

client = Client()

def retrieve_documents(query, top_k=3):
    query_vec = embed_model.encode([query], convert_to_numpy=True)
    distances, indices = index.search(query_vec, top_k)
    docs = []
    for idx in indices[0]:
        docs.append(metadata[idx]['text'])
    return docs

def generate_answer_with_ollama(query, docs):
    context = "\n\n".join(docs)
    messages = [
        {"role": "system", "content": "You are a helpful assistant who answers questions using only the provided context."},
        {"role": "user", "content": f"Context:\n{context}\n\nQuestion: {query}\nAnswer:"}
    ]
    response = client.chat(
        model="llama2",
        messages=messages
    )
    return response.message.content


class SimpleHTTPRequestHandler(BaseHTTPRequestHandler):
    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.end_headers()

    def do_POST(self):
        if self.path == '/query':
            content_length = int(self.headers['Content-Length'])
            post_data = self.rfile.read(content_length)
            data = json.loads(post_data)
            query = data.get('query', '')

            docs = retrieve_documents(query)
            answer = generate_answer_with_ollama(query, docs)

            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.send_header('Access-Control-Allow-Origin', '*')
            self.end_headers()

            response = {'answer': answer}
            self.wfile.write(json.dumps(response).encode())
        else:
            self.send_error(404)

def run(server_class=HTTPServer, handler_class=SimpleHTTPRequestHandler):
    server_address = ('', 8000)
    httpd = server_class(server_address, handler_class)
    print("Running RAG chatbot at http://localhost:8000")
    httpd.serve_forever()

if __name__ == '__main__':
    run()
