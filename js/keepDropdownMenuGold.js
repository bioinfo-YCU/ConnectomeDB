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
        
        // Check for downloads pattern (downloads/XXX/XXX/index.html or index.qmd)
        const downloadsPattern = /\/downloads\/[^\/]+\/[^\/]+\/index\.(html|qmd)/;
        const simpleDownloadsCheck = currentPath.includes('/downloads/');
        
        if (downloadsPattern.test(currentPath) || simpleDownloadsCheck) {
            console.log('Downloads path detected - applying gold styling');
            
            // Function to apply gold styling to downloads nav item
            const applyDownloadsGoldStyling = () => {
                // Get all nav-link elements
                const navLinks = document.querySelectorAll('.nav-link');
                
                navLinks.forEach(navLink => {
                    // Look for .menu-text span inside the nav-link
                    const menuTextSpan = navLink.querySelector('.menu-text');
                    const linkText = menuTextSpan ? menuTextSpan.textContent?.trim() : navLink.textContent?.trim();
                    
                    // Check if this is the Downloads link
                    if (linkText && linkText.toLowerCase().includes('download')) {
                        console.log('Found Downloads nav link, applying gold styling');
                        
                        // Add active class to parent nav-item if it exists
                        const parentNavItem = navLink.closest('.nav-item');
                        if (parentNavItem) {
                            parentNavItem.classList.add("active");
                        }
                        
                        // Apply gold styling to the nav-link
                        navLink.style.setProperty('color', 'gold', 'important');
                        navLink.style.setProperty('font-weight', 'bold', 'important');
                        navLink.style.setProperty('background-color', 'rgba(255, 215, 0, 0.1)', 'important');
                        navLink.style.setProperty('text-shadow', '0 0 3px rgba(255, 215, 0, 0.5)', 'important');
                        navLink.style.setProperty('border-radius', '4px', 'important');
                        navLink.style.setProperty('padding', '8px 12px', 'important');
                        
                        // Also style the .menu-text span if it exists
                        if (menuTextSpan) {
                            menuTextSpan.style.setProperty('color', 'gold', 'important');
                            menuTextSpan.style.setProperty('font-weight', 'bold', 'important');
                        }
                        
                        console.log('Gold styling applied to Downloads');
                    }
                });
            };
            
            // Try immediately and with small delays to ensure elements are loaded
            applyDownloadsGoldStyling();
            setTimeout(applyDownloadsGoldStyling, 100);
            setTimeout(applyDownloadsGoldStyling, 500);
        }
        
        // Check if the path contains "database/" for Ligand-Receptor Browser
        else if (currentPath.includes("/database/")) {
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