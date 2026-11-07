function(el, x, data) {
  'use strict';
  const box = document.createElement('section');
  box.className = 'grip-surface-animation';
  box.innerHTML = `
    <style>
      .grip-surface-animation {border:1px solid #d5e1e5;border-radius:8px;padding:12px;margin:20px 0;background:white;color:#21343e}
      .grip-surface-animation .controls {display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin:8px 0}
      .grip-surface-animation button,.grip-surface-animation select {font:inherit;padding:5px 8px}
      .grip-surface-animation .views {display:grid;grid-template-columns:minmax(0,1.7fr) minmax(230px,1fr);align-items:center}
      .grip-surface-animation .scene {height:490px;min-width:0}
      .grip-surface-animation svg {width:100%;font:12px sans-serif}
      .grip-surface-animation .error {font-size:22px;color:#007e87}
      .grip-surface-animation input[type=range] {flex:1;min-width:150px}
      @media(max-width:750px){.grip-surface-animation .views{grid-template-columns:1fr}}
    </style>
    <div class="controls">
      <label>Surface <select data-role="surface"><option value="paraboloid">Paraboloid</option><option value="saddle">Saddle</option></select></label>
      <button data-role="play">Start</button><button data-role="stop">Stop</button>
      <button data-role="reset">Reset</button>
      <label>Playback <select data-role="speed"><option value="1">Normal</option><option value="4">4×</option><option value="16">16×</option></select></label>
    </div>
    <div class="controls">
      <label>Stage <select data-role="stage"><option value="0">Surface</option><option value="1">Sample points</option><option value="2">Ambient chords</option><option value="3">Aligned random start</option><option value="4">SGD updates</option><option value="5">Final reconstruction</option></select></label>
      <label><input type="checkbox" data-role="mesh" checked> Surface</label>
      <label><input type="checkbox" data-role="truth"> Reference points</label>
      <label><input type="checkbox" data-role="chords"> 1,000 chords</label>
    </div>
    <div class="views"><div class="scene"></div><div>
      <p>Mean squared relative error × 100</p><p class="error" data-role="error">—</p>
      <svg viewBox="0 0 350 240" role="img" aria-label="Relative stress by epoch"></svg>
      <p data-role="readout"></p>
    </div></div>
    <div class="controls"><label>Saved update</label><input aria-label="Saved update" data-role="frame" type="range" min="0" value="0" step="1"></div>
    <p data-role="description" aria-live="polite"></p>`;
  el.parentNode.insertBefore(box, el);
  const $ = role => box.querySelector(`[data-role="${role}"]`);
  const host = box.querySelector('.scene');
  host.appendChild(el);
  const scene = el.rglinstance, root = scene.scene.rootSubscene;
  const sub = scene.getObj(root), originalBox = sub.par3d.bbox.slice();
  const radius = b => Math.sqrt(((b[1]-b[0])/2)**2 + ((b[3]-b[2])/2)**2 + ((b[5]-b[4])/2)**2)*1.1;
  const eyeRatio = sub.par3d.observer[2]/radius(originalBox);
  const datasets = data.surfaces;
  let D = datasets.paraboloid, stage = 0, frame = 0, playing = false, timer;
  let shown = D.truth;
  function visible(id, on) {
    if (on) scene.addToSubscene(id, root); else scene.delFromSubscene(id, root);
  }
  function coordinates(points) {
    const n = points.length;
    const vals = [0,1,2].flatMap(k => points.map(p => p[k]));
    scene.vertexSetter(null, {
      objid:D.ids.points, vertices:Array.from({length:n*3},(_,i)=>i%n),
      attributes:Array.from({length:n*3},(_,i)=>['x','y','z'][Math.floor(i/n)]),
      param:[-1,0,1],values:[vals,vals,vals],value:0,interp:false
    });
    scene.getObj(D.ids.points).vertices = points.map(p=>p.slice());
    shown = points;
  }
  function chart() {
    const svg = box.querySelector('svg');
    const L=52, R=335, T=15, B=205, lo=-15, hi=3;
    const px = update => L+(R-L)*update/(D.epochs*D['pair.count']);
    const py = error => T+(B-T)*(hi-Math.log10(Math.max(1e-15,100*error)))/(hi-lo);
    let s='';
    for(let y=lo;y<=hi;y+=3) {
      const pos=T+(B-T)*(hi-y)/(hi-lo);
      s+=`<line x1="${L}" x2="${R}" y1="${pos}" y2="${pos}" stroke="#e4e9eb"/><text x="${L-5}" y="${pos+4}" text-anchor="end">10^${y}</text>`;
    }
    for(let i=0;i<=D.epochs;i++)
      s+=`<text x="${px(i*D['pair.count'])}" y="223" text-anchor="middle">${i}</text>`;
    const threshold=py(1e-14);
    s+=`<line x1="${L}" x2="${R}" y1="${threshold}" y2="${threshold}" stroke="#8b969e" stroke-dasharray="4 3"/>`;
    if(stage>=3) {
      const path=D.errors.slice(0,frame+1).map((v,i)=>(i?'L':'M')+px(D.updates[i])+','+py(v)).join(' ');
      s+=`<path data-series="stress" d="${path}" fill="none" stroke="#007e87" stroke-width="2"/><circle cx="${px(D.updates[frame])}" cy="${py(D.errors[frame])}" r="3" fill="#007e87"/>`;
    }
    svg.innerHTML=s+'<text x="194" y="239" text-anchor="middle">Epoch</text>';
  }
  function render() {
    const coords=stage<3?D.truth:D.frames[frame];
    coordinates(coords);
    const bounds=D.bounds.flat();
    sub.par3d.bbox=bounds; sub.par3d.observer[2]=eyeRatio*radius(bounds);
    visible(D.ids.points,stage>0);
    visible(D.ids.surface,stage===0||$('mesh').checked);
    visible(D.ids.chords,stage===2||$('chords').checked);
    visible(D.ids.truth,stage>=3&&$('truth').checked);
    scene.drawScene();
    $('stage').value=stage; $('frame').max=D.frames.length-1; $('frame').value=frame;
    $('play').textContent=playing?'Pause':stage===5?'Replay':'Start';
    $('stop').disabled=!playing;
    $('error').textContent=stage>=3?(100*D.errors[frame]).toExponential(4):'—';
    const updates=D.updates[frame], count=D['pair.count'];
    $('readout').textContent=stage>=3?
      `${updates.toLocaleString()} pair updates · ${(updates/count).toFixed(3)} epochs. RMS relative error: ${(100*Math.sqrt(D.errors[frame])).toExponential(3)}%.`:
      `${D.truth.length} points; ${count.toLocaleString()} target pairs.`;
    const descriptions=[
      'Both surfaces have two intrinsic dimensions. Their ambient distances are measured in three dimensions.',
      'Points are sampled uniformly over surface area. Colors retain original height; drag to rotate and scroll to zoom.',
      'Only 1,000 chords are drawn. All 31,125 pairs are used for fitting and evaluation.',
      'The random start is rigidly aligned to the truth, without rescaling.',
      'Actual saved SGD states, independently aligned to the truth. No coordinates or errors are interpolated.',
      `Stopped after ${D.epochs} epochs (${D['stop.reason']}). The dashed line is the 100 S = 10⁻¹² stopping threshold.`
    ];
    $('description').textContent=descriptions[stage]; chart();
  }
  function stop(){playing=false;clearTimeout(timer);}
  function advance() {
    if(!playing)return;
    let delay=1600;
    if(stage===4) {
      const epoch=Math.floor(D.updates[frame]/D['pair.count']);
      const steps=D.updates.filter(u=>u>epoch*D['pair.count']&&u<=(epoch+1)*D['pair.count']).length;
      delay=(epoch===0?16000:4000)/Math.max(1,steps);
    }
    timer=setTimeout(()=>{
      if(stage<4)stage++;
      else if(frame<D.frames.length-1)frame++;
      else{stage=5;stop();}
      render();advance();
    },Math.max(4,delay/Number($('speed').value)));
  }
  function setStage(value) {
    stop();stage=value;
    if(stage<=3)frame=0;
    if(stage===5)frame=D.frames.length-1;
    render();
  }
  $('play').onclick=()=>{
    if(playing)stop();
    else{if(stage===5){stage=0;frame=0;}playing=true;advance();}
    render();
  };
  $('stop').onclick=()=>{stop();render();};
  $('reset').onclick=()=>setStage(0);
  $('stage').onchange=()=>setStage(Number($('stage').value));
  $('frame').oninput=()=>{stop();frame=Number($('frame').value);stage=frame===0?3:4;render();};
  $('surface').onchange=()=>{
    stop();Object.values(D.ids).forEach(id=>visible(id,false));
    D=datasets[$('surface').value];stage=0;frame=0;render();
  };
  ['mesh','truth','chords'].forEach(role=>$(role).onchange=render);
  $('speed').onchange=()=>{if(playing){clearTimeout(timer);advance();}};
  Object.values(datasets.saddle.ids).forEach(id=>visible(id,false));
  const resize=()=>{
    const bounds=host.getBoundingClientRect();
    el.style.width=bounds.width+'px';el.style.height=bounds.height+'px';
    scene.canvas.width=bounds.width;scene.canvas.height=bounds.height;scene.drawScene();
  };
  new ResizeObserver(resize).observe(host);
  // Inspection is attached to this widget, so multiple instances do not clash.
  el.surfaceAnimation={datasets,scene,setStage,get state(){return {surface:D.key,stage,frame,playing,shown,error:D.errors[frame],updates:D.updates[frame]};}};
  render();resize();
}
