document.addEventListener("DOMContentLoaded", function () {
    // Get the current URL path
    const currentPath = window.location.pathname;
    
    // Documentation pages (in main directory)
    const documentationPages = [
        'policy.html',
        'about.html', 
        'methods.html',
        'tutorials.html',
        'help.html'
    ];
    
    // Check if current page is a documentation page
    const isDocumentationPage = documentationPages.some(page => currentPath.includes(page));
    
    // Check if the path contains "database/" for Ligand-Receptor Browser
    if (currentPath.includes("/database/")) {
        // Find the 'Ligand-Receptor Browser' dropdown and mark it as active
        const dropdownItems = document.querySelectorAll(".navbar-nav .nav-item.dropdown");
        dropdownItems.forEach(item => {
            const dropdownText = item.querySelector('.nav-link')?.textContent?.trim();
            if (dropdownText && dropdownText.includes('Ligand-Receptor Browser')) {
                item.classList.add("active");
            }
        });
    } 
    // Check if current page is a documentation page
    else if (isDocumentationPage) {
        // Find the 'Documentation' dropdown and mark it as active
        const dropdownItems = document.querySelectorAll(".navbar-nav .nav-item.dropdown");
        dropdownItems.forEach(item => {
            const dropdownText = item.querySelector('.nav-link')?.textContent?.trim();
            if (dropdownText && dropdownText.includes('Documentation')) {
                item.classList.add("active");
            }
        });
    }
});