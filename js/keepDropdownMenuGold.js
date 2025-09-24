document.addEventListener("DOMContentLoaded", function () {
    try {
        // Get the current URL path
        const currentPath = window.location.pathname;
        console.log('Current path:', currentPath);
        
        // Documentation pages (in main directory)
        const documentationPages = [
            'policy.html',
            'about.html', 
            'methods.html',
            'tutorials.html',
            'help.html',
            'citeus.html',
            'updates.html'
        ];
        
        // Check if current page is a documentation page
        const isDocumentationPage = documentationPages.some(page => currentPath.includes(page));
        
        // Check if the path contains "database/" for Ligand-Receptor Browser
        if (currentPath.includes("/database/")) {
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
            const dropdownItems = document.querySelectorAll(".navbar-nav .nav-item.dropdown");
            dropdownItems.forEach(item => {
                const dropdownText = item.querySelector('.nav-link')?.textContent?.trim();
                if (dropdownText && dropdownText.includes('Documentation')) {
                    item.classList.add("active");
                }
            });
        }
        
    } catch (error) {
        console.error('Error in navigation script:', error);
    }
});