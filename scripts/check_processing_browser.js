// Controlled browser checks for B1 display semantics. Run after opening the local viewer.
async page => {
  const assert = (ok, message) => { if (!ok) throw new Error(message); };
  const origin = new URL(page.url()).origin;
  let association = 'far_from_prediction';
  await page.unroute('**/api/**');
  await page.route('**/api/**', route => {
    if (new URL(route.request().url()).pathname === '/api/targets') return route.fulfill({json: []});
    const frame = i => ({date:'2026-01-0'+(i+1),image_url:'/cache/processing-qa.svg',flags_url:null,
      predicted_px:association==='far_from_prediction'?[4,5]:null,measured_px:null,
      measured_candidates:[{xy:[6,7],kind:'central',quality:'fair'}],
      measurement:{status:'ok',quality:'fair',target_association:association,
        offset_from_prediction_px:association==='far_from_prediction'?16.9:null},
      prediction_status:association==='far_from_prediction'?'in_frame':'outside_field',
      quality:{coverage_fraction:.8,status:'limited',verified:false,reasons:[]}});
    return route.fulfill({json:{status:'ready',query:{ra:150,dec:20},grid:{width_px:20,height_px:20},
      frames:[frame(0),frame(1)],warnings:['Registration unverified; different passbands do not prove variability.'],processing:null}});
  });
  await page.route('**/cache/processing-qa.svg', route => route.fulfill({contentType:'image/svg+xml',
    body:'<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"><rect width="20" height="20" fill="gray"/></svg>'}));
  await page.goto(origin+'/?ra=150&dec=20');
  await page.waitForFunction(()=>document.getElementById('info').textContent.startsWith('Frame'));
  assert((await page.locator('#info').textContent()).includes('candidate 16.9 px from prediction; target identity unverified'),'association label');
  assert((await page.locator('#frame-quality').textContent()).includes('Violet pixels have no usable data'),'missing-data legend');
  const warning = await page.locator('#warnings').textContent();
  await page.getByRole('checkbox',{name:'Show centroid / candidates'}).check();
  assert(await page.locator('.candidate').count()===1,'candidate hidden by far association');
  await page.getByRole('button',{name:'Next date →',exact:true}).click();
  assert((await page.locator('#frame-date').textContent()).includes('2026-01-02'),'date switch');
  assert(await page.locator('#warnings').textContent()===warning,'warning persistence');
  association = 'prediction_outside_field';
  await page.reload();
  await page.waitForFunction(()=>document.getElementById('info').textContent.includes('target prediction outside field'));
  assert(!(await page.locator('#info').textContent()).includes('centroid candidate: fair'),'outside field mislabeled as target centroid');
  await page.getByRole('checkbox',{name:'Show centroid / candidates'}).check();
  await page.getByRole('checkbox',{name:'Show expected position'}).check();
  assert(await page.locator('.candidate').count()===1,'outside-field candidates lost');
  assert(await page.locator('#marker').isHidden() && await page.locator('#measured').isHidden(),'outside-field target markers visible');
  const widths = [];
  for (const width of [1440,768,375]) {
    await page.setViewportSize({width,height:900});
    const actual = await page.locator('#frame').evaluate(el=>el.getBoundingClientRect().width);
    assert(actual <= width,'image overflow at '+width);
    assert(await page.evaluate(()=>document.documentElement.scrollWidth<=window.innerWidth),'page overflow at '+width);
    widths.push({viewport:width,image_width:actual});
  }
  await page.unroute('**/api/**');
  await page.unroute('**/cache/processing-qa.svg');
  return {provenance:'controlled synthetic browser fixtures; integrated frontend with outside-field association',
    passed:['far candidate labelled without identity claim','violet no-data legend','candidate retained','date switch','persistent warnings',
      'outside-field label','outside-field candidates retained','outside-field target markers hidden','page fits desktop and mobile'],
    widths};
}
