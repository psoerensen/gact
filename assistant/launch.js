'use strict';
for (const button of document.querySelectorAll('[data-task]')) {
  button.addEventListener('click', async () => {
    const status = document.getElementById('copy-status');
    const prompt = document.getElementById('instructions').textContent;
    try {
      await navigator.clipboard.writeText(prompt + '\nInitial task: ' + button.dataset.task);
      status.textContent = 'Complete instructions and selected task copied. Paste into your own AI conversation.';
    } catch (error) {
      document.getElementById('portable').open = true;
      status.textContent = 'Clipboard unavailable. Copy the expanded instructions manually, then add this task: ' + button.dataset.task;
    }
  });
}
