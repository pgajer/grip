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
      if (await page.locator('.ivue-legend details').count()) throw new Error(`Legend table remains: ${file}`);
      await page.evaluate(() => document.fonts.ready);
      const wrappedLabels = await page.locator('.ivue-legend > div > span:last-child').evaluateAll(nodes =>
        nodes.filter(node => {
          const range = document.createRange(); range.selectNodeContents(node);
          return range.getClientRects().length > 1;
        }).map(node => node.textContent));
      if (wrappedLabels.length) throw new Error(`Wrapped legend entries: ${file}: ${wrappedLabels.join(', ')}`);
      const selector = page.locator('.grip-layout-selector');
      const hasSelector = await selector.count();
      if (hasSelector) {
        if (await page.locator('.ivue-tools button').count()) throw new Error(`Old view buttons remain: ${file}`);
        const choices = await selector.locator('option').allTextContents();
        if (JSON.stringify(choices) !== JSON.stringify(['metric-MDS', 'metric-MDS + edge-KK', 'both']))
          throw new Error(`Incorrect layout choices: ${file}`);
        const inspect = () => page.evaluate(() => {
          const el = document.querySelector('.rglWebGL'), scene = el.rglinstance;
          const root = scene.getObj(scene.scene.rootSubscene);
          const payload = JSON.parse(document.querySelector('script[data-for="' + el.id + '"]').textContent);
          const data = payload.jsHooks.render.find(hook => hook.data?.ids).data;
          return {ids: data.ids, shown: root.objects,
            geometry: JSON.stringify(Object.values(data.ids).flat().map(id => scene.getObj(id).vertices)),
            camera: JSON.stringify({matrix: root.par3d.userMatrix.getAsArray(), zoom: root.par3d.zoom, bbox: root.par3d.bbox}),
            legends: Array.from(el.querySelectorAll('.ivue-legend > div')).filter(row => row.style.display !== 'none').map(row => row.textContent)};
        });
        const initial = await inspect();
        for (const mode of ['mds', 'refined', 'both', 'mds', 'both']) {
          await selector.selectOption(mode);
          const state = await inspect();
          for (const [name, ids] of Object.entries(state.ids)) {
            const visible = name === 'reference' || mode === 'both' || mode === name;
            if (ids.some(id => state.shown.includes(id) !== visible)) throw new Error(`Wrong object visibility: ${file} ${mode} ${name}`);
          }
          const expectedLegends = (state.ids.reference ? 1 : 0) + (mode === 'both' ? 2 : 1);
          if (state.legends.length !== expectedLegends) throw new Error(`Wrong legend visibility: ${file} ${mode}`);
          if (state.camera !== initial.camera || state.geometry !== initial.geometry)
            throw new Error(`Switching layouts changed camera or coordinates: ${file}`);
        }
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
      results.push({file, views: expectedViews, offline: true, rotation: true, layoutSelector: Boolean(hasSelector), legendTables: false, javascriptErrors: errors});
      await page.close();
      console.log(`Captured and checked ${file}`);
    }
  } finally { await browser.close(); }
  fs.writeFileSync(path.join(source, 'browser-validation.json'), JSON.stringify(results, null, 2) + '\n');
})().catch(e => {console.error(e); process.exit(1);});
