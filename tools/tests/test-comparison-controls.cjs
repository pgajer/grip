// Run the actual serialized selector hooks without a browser or WebGL.
// Usage: node tools/tests/test-comparison-controls.cjs path/to/vignette.html ...
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');

class Element {
  constructor(tag) {
    this.tagName = tag.toUpperCase(); this.children = []; this.style = {};
    this.events = {}; this.textContent = ''; this.value = '';
  }
  appendChild(child) { this.children.push(child); return child; }
  get options() { return this.children; }
  setAttribute() {}
  addEventListener(name, fn) { this.events[name] = fn; }
  querySelector() { return null; }
  querySelectorAll() { return []; }
  remove() {}
}

let checked = 0, singleObject = 0, multipleObjects = 0;
assert.ok(process.argv.length > 2, 'Supply a rendered comparison vignette');
for (const file of process.argv.slice(2)) {
  const html = fs.readFileSync(file, 'utf8');
  for (const match of html.matchAll(/<script\b([^>]*)>([\s\S]*?)<\/script>/g)) {
    if (!/type="application\/json"/.test(match[1]) || !/data-for=/.test(match[1])) continue;
    const widget = JSON.parse(match[2]);
    const hook = (widget.jsHooks?.render || []).find(h => h.code.includes('grip-layout-selector'));
    if (!hook) continue;
    const ids = hook.data.ids;
    const root = widget.x.objects[widget.x.rootSubscene];
    const visible = new Set([].concat(root.objects));
    for (const values of Object.values(ids)) {
      if (Array.isArray(values)) multipleObjects++; else singleObject++;
      for (const id of [].concat(values)) {
        assert.ok(widget.x.objects[id], `Missing drawing object ${id}`);
        assert.ok(visible.has(id), `Object ${id} is absent from the initial scene`);
      }
    }
    const scene = {
      scene: widget.x,
      addToSubscene(id) { visible.add(id); },
      delFromSubscene(id) { visible.delete(id); },
      drawScene() {}
    };
    const el = new Element('div'); el.rglinstance = scene;
    el.after = panel => { el.panel = panel; };
    const run = vm.runInNewContext('(' + hook.code + ')', {
      document: { createElement: tag => new Element(tag) },
      ResizeObserver: class { observe() {} disconnect() {} }
    });
    run(el, widget.x, hook.data);
    const select = el.panel.children[0].children[0];
    assert.equal(el.panel.style.position, 'static', 'Controls must sit below the canvas');
    const preferred = ['mds', 'SGD', 'SGD: uniform', 'Full SGD'].find(name => name in ids);
    if (preferred) assert.equal(select.value, preferred, 'Metric MDS must be the initial layout');
    function checkVisibility(choice) {
      for (const [name, values] of Object.entries(ids)) {
        const expected = ['Reference', 'reference'].includes(name) ||
          ['all', 'both'].includes(choice) || choice === name;
        for (const id of [].concat(values)) assert.equal(visible.has(id), expected,
          `${file}: ${name} visibility with ${choice}`);
      }
    }
    checkVisibility(select.value);
    for (const choice of select.options) {
      if (choice.disabled) continue;
      select.value = choice.value; select.events.change();
      checkVisibility(select.value);
    }
    checked++;
  }
}
assert.ok(checked > 0, 'No selector hooks tested');
assert.ok(singleObject > 0, 'Include a cloud-only figure to cover scalar IDs');
assert.ok(multipleObjects > 0, 'Include figures with edges or meshes');
console.log(`Passed ${checked} selectors: ${singleObject} single-object and ${multipleObjects} multi-object layouts.`);
