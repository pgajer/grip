// Run after render-S4.R. PLAYWRIGHT_MODULE and CHROMIUM_EXECUTABLE are optional.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const path=require('path');
const fs=require('fs');
(async()=>{
  const root=path.resolve('papers/grip-software-paper');
  const browser=await chromium.launch({headless:true,
    ...(process.env.CHROMIUM_EXECUTABLE ? {executablePath:process.env.CHROMIUM_EXECUTABLE}:{}),
    args:['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
  const page=await browser.newPage({viewport:{width:720,height:640},deviceScaleFactor:2});
  const errors=[]; page.on('pageerror',e=>errors.push(String(e)));
  await page.route(/^https?:/,r=>r.abort()); // All assets must work offline.
  const output=path.join(root,'reproducibility/figures/saddle');
  fs.mkdirSync(output,{recursive:true});
  for(const name of ['mesh1','mesh2','mesh3','overlay','displacement']) {
    await page.goto('file://'+path.join(root,'build/saddle-widgets',name+'.html'));
    await page.waitForFunction(()=>document.querySelector('canvas')?.width>0);
    await page.waitForTimeout(2500);
    await page.screenshot({path:path.join(output,name+'.png'),
      clip:{x:100,y:125,width:520,height:390}});
  }
  await page.goto('file://'+path.join(root,'supplement/S4-interactive-saddle.html'));
  await page.waitForFunction(()=>document.querySelectorAll('canvas').length===8);
  await page.waitForTimeout(2000);
  // Exercise every scene offline, including below-the-fold views.
  for(let i=0;i<8;i++) {
    const canvas=page.locator('canvas').nth(i);
    await canvas.scrollIntoViewIfNeeded();
    await page.waitForTimeout(150);
    const before=await canvas.screenshot();
    const box=await canvas.boundingBox();
    await page.mouse.move(box.x+box.width*.55,box.y+box.height*.55);
    await page.mouse.down();
    await page.mouse.move(box.x+box.width*.75,box.y+box.height*.35,{steps:8});
    await page.mouse.up(); await page.waitForTimeout(150);
    const after=await canvas.screenshot();
    if(before.equals(after)) throw new Error('Scene '+i+' did not respond to rotation');
  }
  if(errors.length) throw new Error(errors.join('\n'));
  console.log('Captured five views; all eight offline scenes respond to rotation, no JavaScript errors.');
  await browser.close();
})().catch(e=>{console.error(e);process.exit(1)});
