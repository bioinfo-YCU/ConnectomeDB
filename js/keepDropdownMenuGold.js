document.addEventListener("DOMContentLoaded", function () {
    // Get the current URL path
    const currentPath = window.location.pathname;

    // Check if the path contains "database/"
    if (currentPath.includes("/database/")) {
        // Find the 'ConnectomeDB' dropdown item and mark it as active
        const dropdownItem = document.querySelector(".navbar-nav .nav-item.dropdown");
        if (dropdownItem) {
            dropdownItem.classList.add("active");
        }
    }
});
