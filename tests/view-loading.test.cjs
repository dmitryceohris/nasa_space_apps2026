const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const html = fs.readFileSync(path.join(__dirname, '../frontend/index.html'), 'utf8');
const code = html.slice(html.indexOf('let reqId = 0, visibleFrameController'), html.indexOf('// ================= VIEW SETTINGS ================='));

function setup() {
  class ImageElement extends EventTarget {
    constructor() { super(); this.style = {}; this.complete = false; this.naturalWidth = 0; }
    removeAttribute() { this.complete = false; this.naturalWidth = 0; }
  }
  const elements = Object.fromEntries(['view-loading', 'image-loading', 'loading-detail', 'wrap', 'view-message',
    'play', 'frame-date', 'date-step', 'flagimg', 'flagbox', 'marker', 'measured', 'warnings', 'target-info',
    'field-center', 'frame-quality', 'method', 'share', 'target', 'ra', 'dec', 'fov', 'det', 'snapshot',
    'info', 'settings-note'].map(id => [id, {style: {}, value: '', textContent: '', setAttribute() {}}]));
  elements.frame = new ImageElement();
  elements.target.value = 'barnard';
  const context = vm.createContext({
    $: id => elements[id], document: {querySelectorAll: () => []},
    AbortController, Promise, Error, clearInterval, location: {protocol: 'https:'},
    syncOverlayControls() {}, plotMessage() {}, syncUrl() {}, syncAladin() {},
    getFrames: async () => [{url: 'frame.png'}], show() {}, englishText: text => text,
    timer: null, frames: [], lastGrid: null, plotId: 0, mode: 'coadd', USE_MOCK: false,
  });
  vm.runInContext(code, context);
  return {elements, context};
}

test('loading stays visible until the displayed image finishes loading', async () => {
  const {elements, context} = setup();
  const pending = context.update();
  await new Promise(setImmediate);
  assert.equal(elements['view-loading'].hidden, false);
  assert.equal(elements['image-loading'].hidden, false);
  elements.frame.complete = true; elements.frame.naturalWidth = 100;
  elements.frame.dispatchEvent(new Event('load'));
  await pending;
  assert.equal(elements['view-loading'].hidden, true);
  assert.equal(elements['image-loading'].hidden, true);
});

test('an older image loading result cannot hide the latest loading state', async () => {
  const {elements, context} = setup();
  const old = context.update();
  await new Promise(setImmediate);
  const current = context.update();
  await new Promise(setImmediate);
  await old;
  assert.equal(elements['view-loading'].hidden, false);
  elements.frame.complete = true; elements.frame.naturalWidth = 100;
  elements.frame.dispatchEvent(new Event('load'));
  await current;
  assert.equal(elements['view-loading'].hidden, true);
});

test('a displayed image error replaces loading with an error message', async () => {
  const {elements, context} = setup();
  const pending = context.update();
  await new Promise(setImmediate);
  elements.frame.dispatchEvent(new Event('error'));
  await pending;
  assert.equal(elements['view-loading'].hidden, true);
  assert.match(elements['view-message'].textContent, /could not be loaded/i);
  assert.equal(elements['view-message'].hidden, false);
});
