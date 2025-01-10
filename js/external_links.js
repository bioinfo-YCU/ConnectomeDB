document.addEventListener("DOMContentLoaded", function () {
  const externalLinks = document.querySelectorAll('a[href^="https"]');
  externalLinks.forEach(link => {
    if (!link.href.includes(window.location.hostname)) {
      link.setAttribute("target", "_blank");
      link.setAttribute("rel", "noopener noreferrer");
    }
  });
});
