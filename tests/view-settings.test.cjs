const {test} = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');

const html = fs.readFileSync(path.join(__dirname, '../frontend/index.html'), 'utf8');
const marker = '// ================= VIEW SETTINGS =================';
const start = html.includes(marker) ? html.indexOf(marker) : html.indexOf('["ra", "dec", "fov", "det"].forEach');
const handlers = html.slice(start, html.indexOf('async function loadPresets()', start));
const manifest = require('../frontend/challenge-demo/snapshots/manifest.json');

function setup(values = {}) {
  const elements = {};
  const defaults = {target:'barnard', ra:'', dec:'', fov:'3', det:'3', snapshot:'2026-11-15'};
  for (const id of ['target','ra','dec','fov','det','snapshot','preset','fov-value','det-band','settings-note']) {
    elements[id] = {value: values[id] ?? defaults[id] ?? '', textContent:'',
      options:[{},{}], selectedOptions:[], setAttribute(){}};
  }
  const calls = [];
  const context = {$: id => elements[id], imageSourceChosen:true, USE_MOCK:false,
    presetList:[{id:'barnard',name:"Barnard's Star",fov:3,det:3}],
    snapshotProfiles:manifest.entries.map(e => e.request), P:new URLSearchParams(),
    update:() => calls.push(Object.fromEntries(Object.entries(elements).map(([id,e])=>[id,e.value])))};
  vm.runInNewContext(handlers, context);
  return {elements,calls,context};
}

test('editing field of view leaves the fixed saved selection and requests the chosen live field', () => {
  const {elements:e,calls} = setup();
  e.fov.value='6'; e.fov.onchange();
  assert.equal(e.snapshot.value,'');
  assert.equal(calls.at(-1).fov,'6');
});

test('editing detector leaves the saved selection and preserves the chosen detector', () => {
  const {elements:e,calls} = setup();
  e.det.value='6'; e.det.onchange();
  assert.equal(e.snapshot.value,'');
  assert.equal(calls.at(-1).det,'6');
});

test('returning to saved observations restores their actual field and detector', () => {
  const {elements:e,calls} = setup({fov:'12',det:'6'});
  e.snapshot.onchange();
  assert.equal(e.fov.value,'3');
  assert.equal(e.det.value,'3');
  assert.equal(calls.at(-1).snapshot,'2026-11-15');
});

test('the field value is visible while dragging and does not submit intermediate requests', () => {
  const {elements:e,calls} = setup();
  e.fov.value='8';
  assert.equal(typeof e.fov.oninput,'function');
  e.fov.oninput();
  assert.match(e['fov-value'].textContent,/8/);
  assert.equal(calls.length,0);
});

test('saved observations recognise the exact WISE coordinate example', () => {
  const {elements:e,calls} = setup({target:'',ra:'08h55m10.831s',dec:'-7d14m42.53s',fov:'12',det:'3'});
  e.snapshot.onchange();
  assert.equal(e.fov.value,'6');
  assert.equal(e.det.value,'6');
  assert.equal(calls.at(-1).ra,'08h55m10.831s');
});

test('changing the field keeps Latest IRSA when a preset is subsequently selected', () => {
  const {elements:e,calls,context} = setup();
  context.imageSourceChosen = false;
  e.fov.value='8'; e.fov.onchange();
  assert.equal(context.imageSourceChosen,true);
  assert.match(e['settings-note'].textContent,/Switched to Latest IRSA/);
  e.preset.selectedOptions=[{value:'barnard',dataset:{fov:'3',det:'3'}}];
  e.preset.onchange();
  assert.equal(e.snapshot.value,'');
  assert.equal(calls.at(-1).snapshot,'');
  assert.equal(e['fov-value'].textContent,'3′');
});

test('an unavailable saved field keeps its coordinates and settings for the API error', () => {
  const {elements:e,calls} = setup({ra:'10',dec:'20',fov:'12',det:'6'});
  e.snapshot.onchange();
  assert.equal(calls.at(-1).fov,'12');
  assert.equal(calls.at(-1).det,'6');
  assert.equal(calls.at(-1).ra,'10');
  assert.equal(calls.at(-1).dec,'20');
  assert.match(e['settings-note'].textContent,/No saved observations match/);
});

test('explicit coordinates take priority over a named target when restoring saved settings', () => {
  const {elements:e,calls} = setup({target:'barnard',ra:'08h55m10.831s',dec:'-7d14m42.53s',fov:'12',det:'3'});
  e.snapshot.onchange();
  assert.equal(calls.at(-1).fov,'6');
  assert.equal(calls.at(-1).det,'6');
});

test('leaving saved mode while its manifest loads does not restore obsolete settings', async () => {
  const {elements:e,calls,context} = setup({fov:'12',det:'6'});
  context.snapshotProfiles = null;
  let finishLoading;
  context.snapshotProfilesReady = new Promise(resolve => { finishLoading=resolve; });
  const pending = e.snapshot.onchange();
  assert.equal(calls.length,0);
  e.fov.value='8'; e.fov.onchange();
  context.snapshotProfiles = manifest.entries.map(e => e.request);
  finishLoading();
  await pending;
  assert.equal(calls.length,1);
  assert.equal(e.snapshot.value,'');
  assert.equal(e.fov.value,'8');
  assert.equal(e.det.value,'6');
});
