document.addEventListener('DOMContentLoaded', () => {
  const checkbox = document.getElementById('popupCheckbox');
  const popup = document.getElementById('popup');
  const closeBtn = document.querySelector('.close-btn');

  // Show popup when checkbox is checked
  checkbox.addEventListener('change', () => {
    popup.style.display = checkbox.checked ? 'flex' : 'none';
  });

  // Close popup when close button is clicked
  closeBtn.addEventListener('click', () => {
    popup.style.display = 'none';
  });

  // Close popup when clicking outside of it
  window.addEventListener('click', (event) => {
    if (event.target === popup) {
      popup.style.display = 'none';
    }
  });
});
