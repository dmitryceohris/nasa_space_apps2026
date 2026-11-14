async page => {
  await page.unroute('**/api/**'); await page.unroute('**/cache/fail.png');
  const base = new URL(page.url()).origin;
  const results=[];
  for (const target of ['barnard','epsindib','wise0855','rwaur','rdra','Proxima Cen','61 Cyg A','Vega']) {
    const q = new URLSearchParams({target,mode:'single',max_per_epoch:1,fov:6,det:target==='wise0855'?6:3});
    await page.goto(base+'/?'+q);
    await page.waitForFunction(()=>document.getElementById('info').textContent.startsWith('Frame'),null,{timeout:30000});
    await page.waitForFunction(()=>document.getElementById('frame').complete && document.getElementById('frame').naturalWidth>0);
    const before=await page.locator('#info').textContent();
    await page.getByRole('button',{name:'→',exact:true}).click();
    if (await page.locator('#info').textContent() === before) throw new Error('Date did not change for '+target);
    results.push({target,before,after:await page.locator('#info').textContent(),loaded:await page.locator('#frame').evaluate(el=>el.naturalWidth),warnings:await page.locator('#warnings').textContent()});
  }
  return {provenance:'real API and real PNGs; server cache allowed',results};
}
