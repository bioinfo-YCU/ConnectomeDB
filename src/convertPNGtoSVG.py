# Convert PNG to SVG

import base64
from PIL import Image

script_dir = os.path.dirname(os.path.abspath(__file__))  # src folder
parent_dir = os.path.dirname(script_dir)  # parent of src
# Parameter (Relative Directory and image name)
filename = "images/data4_mobile"
filename = parent_dir + filename

def png_to_svg(png_path, svg_path):
    # Read and encode PNG
    with open(png_path, "rb") as img_file:
        img_data = base64.b64encode(img_file.read()).decode()
    
    # Get dimensions
    img = Image.open(png_path)
    width, height = img.size
    
    # Create SVG
    svg_content = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">
  <image href="data:image/png;base64,{img_data}" width="{width}" height="{height}"/>
</svg>'''
    
    with open(svg_path, "w") as f:
        f.write(svg_content)

# Usage
png_to_svg(filename+".png", filename+".svg")