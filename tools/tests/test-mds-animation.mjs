// Test the rendered vignette in isolated headless Chrome.
// Usage: node tools/tests/test-mds-animation.mjs path/to/metric-mds-animation.html
// Set CHROME_BIN when Chrome is installed elsewhere.
import {spawn} from 'node:child_process';
import {readFile,writeFile,mkdtemp,rm} from 'node:fs/promises';
import {join,dirname,resolve} from 'node:path';
import {tmpdir} from 'node:os';
import {pathToFileURL} from 'node:url';
import assert from 'node:assert/strict';
const file=resolve(process.argv[2] || 'output/vignette-previews/doc/metric-mds-animation.html');
const out=dirname(file),url=pathToFileURL(file).href;

const profile=await mkdtemp(join(tmpdir(),'presentation-assets-qa-'));
const chrome=spawn(process.env.CHROME_BIN || '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',['--headless=new','--enable-unsafe-swiftshader','--no-first-run','--no-default-browser-check','--disable-background-networking','--disable-extensions','--remote-debugging-port=0',`--user-data-dir=${profile}`,'about:blank'],{stdio:'ignore'});
let ws;const checks=[],errors=[];const pause=ms=>new Promise(r=>setTimeout(r,ms));
async function until(fn){for(let i=0;i<200;i++){try{const v=await fn();if(v)return v;}catch{}await pause(100);}throw Error('Browser timed out');}
try{
 const port=await until(async()=>Number((await readFile(join(profile,'DevToolsActivePort'),'utf8')).split('\n')[0]));
 const info=await(await fetch(`http://127.0.0.1:${port}/json/version`)).json();
 ws=new WebSocket(info.webSocketDebuggerUrl);await new Promise((r,j)=>{ws.onopen=r;ws.onerror=j;});
 let seq=0;const pending=new Map();
 ws.onmessage=e=>{const m=JSON.parse(e.data);if(m.id){const p=pending.get(m.id);pending.delete(m.id);m.error?p.reject(Error(JSON.stringify(m.error))):p.resolve(m.result);}else if(m.method==='Runtime.exceptionThrown')errors.push(m.params.exceptionDetails);};
 const call=(method,params={},sessionId)=>new Promise((resolve,reject)=>{const id=++seq;pending.set(id,{resolve,reject});ws.send(JSON.stringify({id,method,params,...(sessionId?{sessionId}:{})}));});
 const {targetId}=await call('Target.createTarget',{url:'about:blank'});const {sessionId}=await call('Target.attachToTarget',{targetId,flatten:true});
 const page=(m,p)=>call(m,p,sessionId);
 const evaluate=async expression=>{const r=await page('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});if(r.exceptionDetails)throw Error(JSON.stringify(r.exceptionDetails));return r.result.value;};
 await page('Runtime.enable');await page('Page.enable');await page('Network.enable');await page('Network.emulateNetworkConditions',{offline:true,latency:0,downloadThroughput:0,uploadThroughput:0});await page('Emulation.setDeviceMetricsOverride',{width:1440,height:1050,deviceScaleFactor:1,mobile:false});await page('Page.navigate',{url});

 await until(async()=>await evaluate('!!document.querySelector(".rglWebGL")?.surfaceAnimation'));
 await evaluate('window.demo=document.querySelector(".rglWebGL").surfaceAnimation;window.control=role=>document.querySelector(`[data-role="${role}"]`)');
 await evaluate('document.querySelector(".grip-surface-animation").scrollIntoView()');
 for(const name of ['paraboloid','saddle']) {
   await evaluate(`control('surface').value='${name}';control('surface').dispatchEvent(new Event('change'))`);
   assert.equal(await evaluate('demo.state.surface'),name);
   assert.equal(await evaluate('demo.state.playing'),false);
   assert.equal(await evaluate('demo.datasets[demo.state.surface].truth.length'),250);
   assert.equal(await evaluate('demo.datasets[demo.state.surface].chords.length'),1000);
   assert.equal(await evaluate('demo.scene.getObj(demo.datasets[demo.state.surface].ids.points).type'),'spheres');
   assert.equal(await evaluate('demo.datasets[demo.state.surface].epochs'),name==='paraboloid'?6:5);
   for(let stage=0;stage<=5;stage++) {
     await evaluate(`demo.setStage(${stage})`);
     assert.equal(await evaluate('demo.state.stage'),stage);
   }
   for(let epoch=0;epoch<(name==='paraboloid'?6:5);epoch++) {
     await evaluate(`(()=>{const D=demo.datasets[demo.state.surface];const i=D.updates.findIndex(u=>u>${epoch}*D['pair.count']);control('frame').value=i;control('frame').dispatchEvent(new Event('input'));})()`);
     assert(await evaluate('JSON.stringify(demo.state.shown)===JSON.stringify(demo.datasets[demo.state.surface].frames[demo.state.frame])'));
     assert(await evaluate(`(()=>{const X=demo.state.shown,Y=demo.datasets[demo.state.surface].truth;let sum=0,n=0;for(let i=0;i<X.length;i++)for(let j=0;j<i;j++){const d=Math.hypot(...Y[i].map((v,k)=>v-Y[j][k]));const r=Math.hypot(...X[i].map((v,k)=>v-X[j][k]))/d-1;sum+=r*r;n++;}return Math.abs(sum/n-demo.state.error)<1e-12;})()`));
     assert(await evaluate(`(()=>{const D=demo.datasets[demo.state.surface],o=demo.scene.getObj(D.ids.points),v=o.vOffsets;return demo.state.shown.every((p,i)=>p.every((x,j)=>Math.abs(x-o.values[i*v.stride+v.vofs+j])<1e-5));})()`));
     assert(await evaluate('document.querySelector("path[data-series=stress]").getAttribute("d").split("L").length===demo.state.frame+1'));
   }
   await evaluate('demo.setStage(3)');
   await pause(200);
   await writeFile(join(out,`review-trace-${name}-start.png`),Buffer.from((await page('Page.captureScreenshot',{format:'png'})).data,'base64'));
   await evaluate('demo.setStage(4);control("speed").value="16";control("play").click()');
   await until(async()=>await evaluate('demo.state.frame>5'));
   await evaluate('control("play").click()');
   const frozen=await evaluate('JSON.stringify(demo.state)');
   await pause(250);assert.equal(await evaluate('JSON.stringify(demo.state)'),frozen);
   await evaluate('control("play").click()');
   await until(async()=>await evaluate('demo.state.stage===5&&!demo.state.playing'));
   assert(await evaluate('100*demo.state.error<=1e-12'));
   await pause(200);
   await writeFile(join(out,`review-trace-${name}-final.png`),Buffer.from((await page('Page.captureScreenshot',{format:'png'})).data,'base64'));
   checks.push(`${name}: all stages, actual states and GPU positions across every epoch, recomputed stress, full curve, pause/resume and completion`);
 }
 assert(await evaluate('demo.datasets.paraboloid.truth.every((p,i)=>p[0]===demo.datasets.saddle.truth[i][0]&&p[1]===demo.datasets.saddle.truth[i][1])'));
 await evaluate('control("reset").click()');assert.equal(await evaluate('demo.state.frame'),0);
 assert(await evaluate('document.documentElement.scrollWidth<=innerWidth'));
 checks.push('Matched sample positions, reset, desktop width, offline playback');
 await page('Emulation.setDeviceMetricsOverride',{width:390,height:844,deviceScaleFactor:1,mobile:true});await pause(200);
 assert(await evaluate('document.documentElement.scrollWidth<=innerWidth'));
 await writeFile(join(out,'review-trace-mobile.png'),Buffer.from((await page('Page.captureScreenshot',{format:'png'})).data,'base64'));
 assert.deepEqual(errors,[]);checks.push('Responsive layout and no browser exceptions');
 await writeFile(join(out,'mds-animation-browser-validation.json'),JSON.stringify({status:'passed',checks},null,2));
 console.log(JSON.stringify({status:'passed',checks}));
}finally{if(ws)ws.close();chrome.kill();await pause(300);await rm(profile,{recursive:true,force:true});}
