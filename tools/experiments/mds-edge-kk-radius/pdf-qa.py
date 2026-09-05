#!/usr/bin/env python3
"""Render final PDF pages and contact sheets for visual review (not deliverables)."""
from pathlib import Path
import subprocess
from PIL import Image,ImageDraw
HERE=Path(__file__).resolve().parent;OUT=HERE.parents[2]/'output/mds-edge-kk-radius/pdf-qa';OUT.mkdir(parents=True,exist_ok=True)
for old in list(OUT.glob('final-*.png'))+list(OUT.glob('contact-*.jpg')):old.unlink()
subprocess.run(['pdftoppm','-scale-to','1000','-png',str(HERE/'report.pdf'),str(OUT/'final')],check=True)
files=sorted(OUT.glob('final-*.png'))
for start in range(0,len(files),8):
 canvas=Image.new('RGB',(2000,4*800),'#dddddd');draw=ImageDraw.Draw(canvas)
 for j,file in enumerate(files[start:start+8]):
  im=Image.open(file).convert('RGB');im.thumbnail((990,765));x=(j%2)*1000;y=(j//2)*800
  canvas.paste(im,(x,y+25));draw.text((x+10,y+5),file.stem,fill='black')
 canvas.save(OUT/f'contact-{start//8+1}.jpg',quality=90)
print(len(files),'pages rendered;',len(range(0,len(files),8)),'contact sheets')
