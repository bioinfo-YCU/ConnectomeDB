function changeTheme() {
    const selectedTheme = document.getElementById("theme-select").value;
    const head = document.getElementsByTagName("head")[0];
    const link = document.createElement("link");
    link.rel = "stylesheet";
    link.type = "text/css";
    link.href = `https://cdn.jsdelivr.net/npm/bootswatch@4.5.2/dist/${selectedTheme}/bootstrap.min.css`;  // Update to your CDN or local path

    // Remove existing theme link if any
    const existingLink = document.querySelector('link[rel="stylesheet"]');
    if (existingLink) {
      head.removeChild(existingLink);
    }

    // Add the new theme
    head.appendChild(link);
}