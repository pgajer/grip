// Usage: node tools/pkg/capture-metric-mds-comparison.cjs <render-dir> <image-dir>
// PLAYWRIGHT_MODULE and CHROMIUM_EXECUTABLE may locate local installations.
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const fs = require('fs');
const path = require('path');
const {pathToFileURL} = require('url');
(async () => {
  const [source, destination] = process.argv.slice(2).map(p => path.resolve(p));
  if (!source || !destination) throw new Error('Supply render and image directories');
  fs.mkdirSync(destination, {recursive: true});
  const browser = await chromium.launch({headless: true,
    ...(process.env.CHROMIUM_EXECUTABLE ? {executablePath: process.env.CHROMIUM_EXECUTABLE} : {}),
    args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader']});
  const results = [];
  try {
    for (const file of fs.readdirSync(source).filter(f => f.endsWith('.html')).sort()) {
      const panels = file.endsWith('-panels.html');
      const page = await browser.newPage({viewport: {width: panels ? 1220 : 920, height: panels ? 450 : 540}, deviceScaleFactor: 1});
      const errors = [];
      page.on('pageerror', e => errors.push(String(e)));
      await page.route(/^https?:/, route => route.abort());
      await page.goto(pathToFileURL(path.join(source, file)).href);
      await page.waitForFunction(() => {
        const el = document.querySelector('.rglWebGL');
        return el?.rglinstance && el.querySelector('canvas')?.width > 0;
      });
      const expectedViews = await page.locator('.rglWebGL').count();
      if (expectedViews < 1 || expectedViews > (panels ? 3 : 1)) throw new Error(`Unexpected view count: ${file}`);
      await page.waitForFunction(expected => [...document.querySelectorAll('.rglWebGL')].filter(el => el.rglinstance && el.querySelector('canvas')?.width > 0).length === expected, expectedViews);
      await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
      if (panels) {
        const cameras = await page.evaluate(() => [...document.querySelectorAll('.rglWebGL')].map(el => {
          const scene=el.rglinstance, par=scene.getObj(scene.scene.rootSubscene).par3d;
          return JSON.stringify({camera:Array.from(par.userMatrix.getAsArray()),zoom:par.zoom,FOV:par.FOV,bbox:par.bbox});
        }));
        if (new Set(cameras).size !== 1) throw new Error(`Panel cameras or bounds differ: ${file}`);
      }
      const output = path.join(destination, file.replace(/\.html$/, '.png'));
      await page.evaluate(() => window.scrollTo(0, 0));
      await page.evaluate(() => document.fonts.ready);
      await page.screenshot({path: output, fullPage: true});
      // Check rotation changes the camera, without changing saved coordinates.
      const readCamera = () => page.evaluate(() => {
        const s = document.querySelector('.rglWebGL').rglinstance;
        return Array.from(s.getObj(s.scene.rootSubscene).par3d.userMatrix.getAsArray());
      });
      const before = await readCamera();
      const canvas = await page.locator('canvas').first().boundingBox();
      await page.mouse.move(canvas.x + canvas.width * .45, canvas.y + canvas.height * .45);
      await page.mouse.down();
      await page.mouse.move(canvas.x + canvas.width * .65, canvas.y + canvas.height * .55, {steps: 8});
      await page.mouse.up();
      const after = await readCamera();
      if (JSON.stringify(before) === JSON.stringify(after)) throw new Error(`Rotation failed: ${file}`);
      if (errors.length) throw new Error(errors.join('\n'));
      results.push({file, views: expectedViews, offline: true, rotation: true, javascriptErrors: errors});
      await page.close();
      console.log(`Captured and checked ${file}`);
    }
  } finally { await browser.close(); }
  fs.writeFileSync(path.join(source, 'browser-validation.json'), JSON.stringify(results, null, 2) + '\n');
})().catch(e => {console.error(e); process.exit(1);});
