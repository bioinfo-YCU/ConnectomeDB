function changeTheme() {
  const selectedTheme = document.getElementById("theme-select").value;
  const head = document.getElementsByTagName("head")[0];
  const rootElement = document.documentElement;

  // Remove any existing theme link
  const existingLink = document.querySelector('link[data-theme="custom-theme"]');
  if (existingLink) {
    head.removeChild(existingLink);
  }

  // Remove any inline styles for CSS variables related to the previous theme
  rootElement.style.removeProperty('--navbar-bg-color');
  rootElement.style.removeProperty('--navbar-text-color');
  rootElement.style.removeProperty('--valuebox-bg-color');
  rootElement.style.removeProperty('--valuebox-text-color');
  rootElement.style.removeProperty('--primary-button-bg-color');
  rootElement.style.removeProperty('--primary-button-text-color');
  rootElement.style.removeProperty('--heading-color');
  rootElement.style.removeProperty('--bslib-color-bg');
  rootElement.style.removeProperty('--bslib-color-fg');
  rootElement.style.removeProperty('--valuebox-primary-bg-color');
  rootElement.style.removeProperty('--valuebox-primary-text-color');

  // Create a new link element for the selected theme
  const link = document.createElement("link");
  link.rel = "stylesheet";
  link.type = "text/css";
  link.href = `https://cdn.jsdelivr.net/npm/bootswatch@4.5.2/dist/${selectedTheme}/bootstrap.min.css?v=${new Date().getTime()}`; // Prevent caching
  link.setAttribute("data-theme", "custom-theme"); // Mark as custom theme

  // Append the new link element to the head
  head.appendChild(link);

  // Define a configuration object for theme properties
  const themeConfig = {
    cosmos: {
      navbarBgColor: '#1d1899',
      navbarTextColor: '#ffffff',
      valueboxBgColor: '#1d1899',
      valueboxTextColor: '#ffffff',
      primaryButtonColor: '#ffffff',
      primaryButtonTextColor: '#1d1899',
      headingColor: '#ffffff',
      valueboxPrimaryBgColor: '#cf0649',  // Background for valuebox.bg-primary
      valueboxPrimaryTextColor: '#ffffff', // Text for valuebox.bg-primary
    },
    flatly: {
      navbarBgColor: '#349442',
      navbarTextColor: '#0f0f0f',
      valueboxBgColor: '#349442',
      valueboxTextColor: '#ffffff',
      primaryButtonColor: '#349442',
      primaryButtonTextColor: '#ffffff',
      headingColor: '#349442',
      valueboxPrimaryBgColor: '#349442',  // Background for valuebox.bg-primary
      valueboxPrimaryTextColor: '#ffffff', // Text for valuebox.bg-primary
    },
    cerulean: {
      navbarBgColor: '#007bff',
      navbarTextColor: '#ffffff',
      valueboxBgColor: '#007bff',
      valueboxTextColor: '#ffffff',
      primaryButtonColor: '#007bff',
      primaryButtonTextColor: '#ffffff',
      valueboxPrimaryBgColor: '#007bff',  // Background for valuebox.bg-primary
      valueboxPrimaryTextColor: '#ffffff', // Text for valuebox.bg-primary
    },
    simplex: {
      navbarBgColor: '#c4353a',
      navbarTextColor: '#ffffff',
      valueboxBgColor: '#c4353a',
      valueboxTextColor: '#ffffff',
      primaryButtonColor: '#c4353a',
      primaryButtonTextColor: '#ffffff',
      valueboxPrimaryBgColor: '#c4353a',  // Background for valuebox.bg-primary
      valueboxPrimaryTextColor: '#ffffff', // Text for valuebox.bg-primary
    },
    default: {
      navbarBgColor: '#63676b', // Default light background
      navbarTextColor: '#000000', // Default dark text color
      valueboxBgColor: '#0f0c02',
      valueboxTextColor: '#000000',
      primaryButtonColor: '#63676b',
      primaryButtonTextColor: '#ffffff',
      valueboxPrimaryBgColor: '#63676b',  // Default background for valuebox.bg-primary
      valueboxPrimaryTextColor: '#ffffff', // Default text color for valuebox.bg-primary
    }
  };

  // Apply styles based on selected theme
  const theme = themeConfig[selectedTheme] || themeConfig.default;

  // Set navbar background and text color based on selected theme
  rootElement.style.setProperty('--navbar-bg-color', theme.navbarBgColor);
  rootElement.style.setProperty('--navbar-text-color', theme.navbarTextColor);
  
  // Set valuebox background and text color based on selected theme
  rootElement.style.setProperty('--valuebox-bg-color', theme.valueboxBgColor);
  rootElement.style.setProperty('--valuebox-text-color', theme.valueboxTextColor);

  // Set primary button background and text color based on selected theme
  rootElement.style.setProperty('--primary-button-bg-color', theme.primaryButtonColor);
  rootElement.style.setProperty('--primary-button-text-color', theme.primaryButtonTextColor);

  // Set valuebox bg-primary colors based on selected theme
  rootElement.style.setProperty('--valuebox-primary-bg-color', theme.valueboxPrimaryBgColor);
  rootElement.style.setProperty('--valuebox-primary-text-color', theme.valueboxPrimaryTextColor);

  // Optionally set heading color (only for themes that define it)
  if (theme.headingColor) {
    rootElement.style.setProperty('--heading-color', theme.headingColor);
  }
}
