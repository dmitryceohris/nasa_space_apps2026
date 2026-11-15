// Run with playwright-cli run-code --filename after opening the local viewer
// with ?target=barnard&fov=3&det=3&snapshot=2026-11-15.
// Saved frames and the sky map are real; live responses are controlled fixtures.
async page => {
  const assert = (condition, message) => { if (!condition) throw new Error(message); };
  const ready = () => page.waitForFunction(() => document.getElementById('info').textContent.startsWith('Frame'));
  const commit = value => page.locator('#fov').evaluate((el, value) => {
    el.value = value;
    el.dispatchEvent(new Event('input', {bubbles:true}));
    el.dispatchEvent(new Event('change', {bubbles:true}));
  }, value);
  const calls = [];
  const passed = [];
  await ready();
  await page.waitForFunction(() => aladin !== null);
  assert(await page.locator('#fov-value').textContent() === '3′', 'URL field size is not visible');
  assert(await page.locator('#frame').evaluate(el => el.naturalWidth > 0), 'saved image did not load');
  passed.push('saved image and initial value from URL');

  let releaseSlow, markSlowSeen;
  const slow = new Promise(resolve => { releaseSlow = resolve; });
  const slowSeen = new Promise(resolve => { markSlowSeen = resolve; });
  const viewRoute = async route => {
    const url = new URL(route.request().url());
    calls.push(url);
    if (url.searchParams.has('snapshot')) return route.continue();
    const fov = Number(url.searchParams.get('fov'));
    if (fov === 8) { markSlowSeen(); await slow; }
    return route.fulfill({json:{status:'ready',query:{ra:269.445979,dec:4.770669,fov},
      grid:{width_px:20,height_px:20,wcs:{CRVAL1:269.445979,CRVAL2:4.770669}},
      frames:[{date:'controlled-field-'+fov,image_url:'/cache/fov-qa.svg',n_exp:1}],warnings:[]}});
  };
  const imageRoute = route => route.fulfill({contentType:'image/svg+xml',
    body:'<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"><rect width="20" height="20" fill="black"/></svg>'});
  await page.route('**/api/view?*', viewRoute);
  await page.route('**/cache/fov-qa.svg', imageRoute);
  try {
    await page.locator('#fov').evaluate(el => { el.value='6'; el.dispatchEvent(new Event('input', {bubbles:true})); });
    assert(await page.locator('#fov-value').textContent() === '6′', 'dragging did not update the value');
    assert(calls.length === 0, 'dragging submitted an intermediate request');
    await page.locator('#fov').dispatchEvent('change');
    await ready();
    assert(await page.locator('#snapshot').inputValue() === '', 'field change did not leave saved mode');
    assert(calls[0].searchParams.get('fov') === '6' && !calls[0].searchParams.has('snapshot'), 'wrong live query');
    assert((await page.locator('#settings-note').textContent()).includes('Switched to Latest IRSA'), 'source change explanation missing');
    const share = new URL(await page.locator('#share-url').inputValue());
    assert(share.searchParams.get('fov') === '6' && !share.searchParams.has('snapshot'), 'share URL has stale settings');
    const mapFov = await page.evaluate(() => aladin.getFov()[0]);
    assert(Math.abs(mapFov - 0.15) < 0.001, 'sky map did not follow the selected field');
    passed.push('drag preview, committed live query, source explanation, share URL and real map update');

    await page.reload(); await ready();
    await page.waitForFunction(() => aladin !== null);
    assert(await page.locator('#fov').inputValue() === '6', 'reload lost the chosen size');
    assert(await page.locator('#fov-value').textContent() === '6′', 'reload lost the numeric label');
    assert(await page.locator('#snapshot').inputValue() === '', 'reload selected saved observations');
    passed.push('live settings survive URL reload');

    await page.locator('#snapshot').selectOption('2026-11-15'); await ready();
    assert(await page.locator('#fov').inputValue() === '3', 'saved size was not restored');
    assert(await page.locator('#det').inputValue() === '3', 'saved detector was not restored');
    assert((await page.locator('#frame').getAttribute('src')).startsWith('/challenge-demo/'), 'saved image was not restored');
    await page.locator('#det').selectOption('6'); await ready();
    assert(await page.locator('#snapshot').inputValue() === '', 'detector change did not leave saved mode');
    await page.locator('#preset').selectOption('epsindib'); await ready();
    assert(await page.locator('#snapshot').inputValue() === '', 'preset unexpectedly returned to saved mode');
    assert(await page.locator('#fov-value').textContent() === '5′', 'preset label was not updated');
    passed.push('return to saved frames, detector change and sticky live source on preset selection');

    await commit('8'); await slowSeen;
    await commit('12'); await ready();
    const lateResponse = page.waitForResponse(r => new URL(r.url()).searchParams.get('fov') === '8');
    releaseSlow(); await lateResponse;
    await page.waitForTimeout(100);
    assert((await page.locator('#info').textContent()).includes('controlled-field-12'), 'late response replaced the chosen field');
    assert(Math.abs(await page.evaluate(() => aladin.getFov()[0]) - 0.3) < 0.001, 'late response changed the map size');
    passed.push('rapid changes ignore late frame responses and retain the map size');

    const slider = await page.locator('#fov').boundingBox();
    await page.mouse.move(slider.x + slider.width * 11 / 19, slider.y + slider.height / 2);
    await page.mouse.down();
    await page.mouse.move(slider.x + slider.width * 15 / 19, slider.y + slider.height / 2, {steps:5});
    await page.mouse.up(); await ready();
    assert(calls.at(-1).searchParams.get('fov') === await page.locator('#fov').inputValue(), 'pointer drag did not load its selected size');
    passed.push('physical pointer drag commits the selected size');
    return {provenance:'real saved frames and Aladin map; controlled live responses',passed};
  } finally {
    releaseSlow();
    await page.unroute('**/api/view?*', viewRoute);
    await page.unroute('**/cache/fov-qa.svg', imageRoute);
  }
}
