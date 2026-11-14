// Run with playwright-cli run-code after opening the local application.
// This is a controlled browser check: API responses and PNGs are synthetic fixtures.
async page => {
  const base = new URL(page.url()).origin;
  const assert = (condition, message) => { if (!condition) throw new Error(message); };
  const calls = [];
  const states = [];
  await page.unroute('**/api/**');
  await page.route('**/api/**', async route => {
    const url = new URL(route.request().url()); calls.push(url);
    if (url.pathname === '/api/targets') return route.fulfill({json: []});
    const ra = url.searchParams.get('ra');
    const target = url.searchParams.get('target');
    if (ra === '10' || target === 'slow') await new Promise(r => setTimeout(r, 350));
    if (ra === '99') return route.fulfill({status:502,json:{error:{code:'upstream',message:'controlled archive failure'}}});
    if (url.pathname !== '/api/view') {
      if (url.pathname === '/api/spectrum') await new Promise(r => setTimeout(r, 350));
      return route.fulfill({json:{status:'ready',warnings:[],wavelength:[1,2],dates:['2026-01-01','2026-01-02'],flux:[1,2],flux_err:[.1,.1]}});
    }
    const bad = ra === '98';
    const date = ra || target;
    const frame = index => ({date:`${date}-${index}`, image_url:bad?'/cache/qa-fail.svg':'/cache/qa.svg',
      flags_url:null, wave_um:2, n_exp:1, predicted_px:[4,5], measured_px:null,
      measured_candidates:[{xy:[6,7],quality:'fair',kind:'central'}],
      measurement:{status:'ambiguous',reason:'controlled candidates'},
      prediction_status:'in_frame',quality:{coverage_fraction:0.8,status:'limited',verified:false,reasons:['incomplete_pixel_coverage']}});
    await route.fulfill({json:{status:'ready',query:{ra:Number(ra)||1,dec:20},
      grid:{width_px:20,height_px:20},frames:[frame(1),frame(2)],warnings:['warning-'+date],
      target_info:null,processing:null}});
  });
  await page.route('**/cache/qa.svg', route => route.fulfill({contentType:'image/svg+xml',body:'<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"><rect width="20" height="20" fill="black"/></svg>'}));
  await page.route('**/cache/qa-fail.svg', route => route.abort());
  const open = query => page.goto(base+'/'+query);
  const ready = () => page.waitForFunction(()=>document.getElementById('info').textContent.startsWith('Frame'));
  await open('');
  assert(calls.every(u=>u.pathname==='/api/targets'), 'empty form requested data');
  states.push('empty form: presets only');
  calls.length = 0;
  await open('?target=barnard&ra=133.795');
  await page.getByRole('button',{name:'Spectrum',exact:true}).click();
  await page.getByRole('button',{name:'Light curve',exact:true}).click();
  assert(calls.every(u=>u.pathname==='/api/targets'), 'partial coordinates fell back to name');
  states.push('partial coordinates: neither graph requested');
  calls.length = 0;
  await open('?ra=08h55m10.831s&dec=-7d14m42.53s&max_per_epoch=1'); await ready();
  const request = calls.find(u=>u.pathname==='/api/view');
  assert(request.searchParams.get('ra')==='08h55m10.831s','RA was converted early');
  assert(request.searchParams.get('max_per_epoch')==='1','max_per_epoch lost');
  assert(!calls.some(u=>/spectrum|lightcurve/.test(u.pathname)), 'graphs loaded automatically');
  await page.getByRole('checkbox',{name:'Show centroid / candidates'}).check();
  const candidate = page.locator('.candidate');
  assert(await candidate.count()===1,'candidate marker missing');
  assert(await candidate.evaluate(e=>e.style.left)==='32.5%','candidate pixel normalization differs');
  const warning = await page.locator('#warnings').textContent();
  await page.getByRole('button',{name:'Next date →',exact:true}).click();
  assert(await page.locator('#warnings').textContent()===warning,'player cleared warnings');
  await page.setViewportSize({width:375,height:760});
  assert(await page.locator('#frame').evaluate(e=>e.getBoundingClientRect().width)<=375,'mobile image overflow');
  states.push('coordinate strings, URL, candidates, persistent warnings and mobile width');
  await open('?ra=10&dec=20');
  await page.getByPlaceholder('RA (deg or hh:mm:ss)').fill('30');
  await page.getByPlaceholder('RA (deg or hh:mm:ss)').press('Tab');
  await ready();
  await page.waitForTimeout(450);
  assert((await page.locator('#info').textContent()).includes('30-1'),'stale frame replaced current frame');
  assert(await page.locator('#warnings').textContent()==='warning-30','stale warnings replaced current warnings');
  await page.getByRole('button',{name:'Play',exact:true}).click();
  await page.getByRole('button',{name:'diff',exact:true}).click(); await ready();
  assert(await page.getByRole('button',{name:'Play',exact:true}).count()===1,'timer not stopped on mode change');
  states.push('late response ignored; mode switch stops player');
  await open('?target=slow');
  await page.getByPlaceholder("Object name, e.g. Barnard's Star").fill('fast');
  await page.getByPlaceholder("Object name, e.g. Barnard's Star").press('Tab');
  await ready(); await page.waitForTimeout(450);
  assert((await page.locator('#info').textContent()).includes('fast-1'),'old target replaced current target');
  states.push('rapid target switch ignores old target');
  await page.getByRole('button',{name:'Spectrum',exact:true}).click();
  await page.getByRole('button',{name:'Light curve',exact:true}).click();
  await page.waitForFunction(()=>document.getElementById('plot').data?.[0]?.x?.[0]==='2026-01-01');
  await page.waitForTimeout(450);
  assert(await page.locator('#plot').evaluate(e=>e.data[0].x[0])==='2026-01-01','old graph replaced current tab');
  states.push('late spectrum cannot replace current light curve');
  await page.evaluate(()=>{delete window.Plotly;});
  await page.getByRole('button',{name:'Spectrum',exact:true}).click();
  await page.waitForFunction(()=>document.getElementById('plot').textContent.includes('Could not load'));
  assert(await page.locator('#frame').evaluate(e=>e.naturalWidth)>0,'missing graph library broke images');
  states.push('missing Plotly does not stop images');
  await open('?ra=98&dec=20');
  await page.waitForFunction(()=>document.getElementById('info').textContent.includes('could not be loaded'));
  assert(await page.locator('#frame').getAttribute('src')===null,'failed PNG left an old image');
  states.push('PNG failure does not report successful display');
  await open('?ra=99&dec=20');
  await page.waitForFunction(()=>document.getElementById('info').textContent.includes('controlled archive failure'));
  assert(await page.locator('#frame').getAttribute('src')===null,'error left an old image');
  assert(await page.locator('#demo-note').isHidden(),'real API failure enabled mock');
  states.push('archive error clears image and does not enable mock');
  await open('?mock=1'); await ready();
  assert(await page.locator('#demo-note').isVisible(),'mock label missing');
  states.push('mock data visibly labelled');
  await page.unroute('**/api/**'); await page.unroute('**/cache/qa.svg'); await page.unroute('**/cache/qa-fail.svg');
  return {provenance:'controlled synthetic browser fixtures',passed:states};
}
