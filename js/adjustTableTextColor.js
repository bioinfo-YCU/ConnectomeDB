function adjustTableTextColor() {
  const tables = document.querySelectorAll('table');
  
  tables.forEach(table => {
    // Get the background color of the table
    const bgColor = window.getComputedStyle(table).backgroundColor;
    
    // Check if the background color is black or dark
    if (bgColor === 'rgb(0, 0, 0)' || bgColor === 'black') {
      table.style.color = 'white';
    } else if (bgColor === 'rgb(255, 255, 255)' || bgColor === 'white') {
      table.style.color = 'black';  // Default to black if background is white
    }
  });
}

// Call on window load
window.onload = adjustTableTextColor;
