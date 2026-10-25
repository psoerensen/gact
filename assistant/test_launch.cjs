// Exercise real copy handlers, including denied clipboard access, without packages.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
(async () => {
  for (const denied of [false, true]) {
    const handlers = [];
    const tasks = ['ingest', 'annotations', 'analyze', 'extend', 'report-a-bug-or-suggest-a-feature'];
    const status = {textContent: ''};
    const portable = {open: false};
    const copied = [];
    const context = {
      document: {
        querySelectorAll: () => tasks.map(task => ({dataset: {task}, addEventListener: (_, fn) => handlers.push(fn)})),
        getElementById: id => ({'copy-status': status, portable, instructions: {textContent: 'shared instructions\n'}}[id])
      },
      navigator: {clipboard: {writeText: async text => {
        if (denied) throw new Error('denied');
        copied.push(text);
      }}}
    };
    vm.runInNewContext(fs.readFileSync(__dirname + '/launch.js', 'utf8'), context);
    for (let i = 0; i < handlers.length; i++) {
      await handlers[i]();
      if (denied) {
        assert.equal(portable.open, true);
        assert.ok(status.textContent.includes(tasks[i]));
      } else {
        assert.equal(copied[i], 'shared instructions\n\nInitial task: ' + tasks[i]);
        assert.ok(status.textContent.includes('copied'));
      }
    }
    assert.equal(handlers.length, 5);
  }
  console.log('All five copy actions and clipboard-denied fallbacks passed.');
})();
