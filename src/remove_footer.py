# remove footer from index.html

import sys, os
from bs4 import BeautifulSoup

def remove_footer():
    # Get the path to index.html in _site folder
    script_dir = os.path.dirname(os.path.abspath(__file__))  # src folder
    parent_dir = os.path.dirname(script_dir)  # parent of src
    html_file = os.path.join(parent_dir, "_site", "index.html")
    
    if not os.path.exists(html_file):
        print(f"File not found: {html_file}")
        return
    
    with open(html_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    soup = BeautifulSoup(content, 'html.parser')
    
    # Remove footer element
    footer = soup.find('footer')
    if footer:
        footer.decompose()
        print("Footer found and removed")
    else:
        print("No footer element found")
    
    # Alternative: remove by class if footer has specific class
    # footer = soup.find('div', class_='footer')  # adjust class name
    # if footer:
    #     footer.decompose()
    
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(str(soup))
    
    print(f"Footer removed from {html_file}")

if __name__ == "__main__":
    remove_footer()