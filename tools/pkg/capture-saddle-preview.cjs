// Run after render-saddle-preview.R. Set PLAYWRIGHT_MODULE or
// CHROMIUM_EXECUTABLE when these dependencies are outside the default paths.
// --animation also captures the 240-frame GIF for local visual inspection.
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const path = require('path');
const fs = require('fs');
const {pathToFileURL} = require('url');
(async () => {
  const out = path.resolve('output/readme-saddle');
  const browser = await chromium.launch({headless: true,
    ...(process.env.CHROMIUM_EXECUTABLE ? {executablePath: process.env.CHROMIUM_EXECUTABLE} : {}),
    args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader']});
  try {
    const page = await browser.newPage({viewport: {width: 1848, height: 820}});
    const errors = [];
    page.on('pageerror', e => errors.push(String(e)));
    await page.route(/^https?:/, r => r.abort());
    await page.goto(pathToFileURL(path.join(out, 'saddle-rotation.html')).href + '?paused');
    await page.waitForFunction(() => document.body.dataset.ready === 'true');
    await page.waitForFunction(() => [...document.querySelectorAll('#figure canvas')]
      .length === 3 && [...document.querySelectorAll('#figure canvas')].every(c => c.width > 0));
    const figure = page.locator('#figure');
    async function frame(i) {
      await page.evaluate(i => window.saddlePreview.setFrame(i), i);
      await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
      await page.evaluate(() => {
        const matrices = window.saddlePreview.scenes.map(s =>
          Array.from(s.getObj(s.scene.rootSubscene).par3d.userMatrix.getAsArray()));
        if (matrices.some(m => JSON.stringify(m) !== JSON.stringify(matrices[0])))
          throw new Error('The three camera matrices differ');
      });
    }
    await frame(0);
    fs.mkdirSync('man/figures', {recursive: true});
    await figure.screenshot({path: 'man/figures/readme-saddle-reference.png'});
    for (const i of [60, 120, 180]) {
      await frame(i);
      await figure.screenshot({path: path.join(out, `quarter-${i}.png`)});
    }
    // Exercise the ivue player itself, rather than checking only direct calls.
    await page.getByRole('button', {name: 'Play', exact: true}).click();
    await page.waitForFunction(() => Number(document.body.dataset.frame) !== 180);
    await page.getByRole('button', {name: 'Pause', exact: true}).click();
    await page.getByRole('button', {name: 'Reset', exact: true}).click();
    await page.waitForFunction(() => document.body.dataset.frame === '0');
    const slider = page.locator('input[type=range]');
    await slider.fill('80');
    await page.waitForFunction(() => document.body.dataset.frame === '80');
    if (process.argv.includes('--animation')) {
      const dir = path.join(out, 'frames');
      fs.mkdirSync(dir, {recursive: true});
      for (let i = 0; i < 240; i++) {
        await frame(i);
        await figure.screenshot({path: path.join(dir, `frame-${String(i).padStart(3, '0')}.png`)});
        if (i % 40 === 0) console.log(`Captured ${i + 1}/240 frames`);
      }
    }
    if (errors.length) throw new Error(errors.join('\n'));
    console.log('Verified three synchronized offline scenes; playback, pause, reset and scrubbing; no JavaScript errors.');
  } finally { await browser.close(); }
})().catch(e => {console.error(e); process.exit(1);});
