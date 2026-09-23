#!/usr/bin/env python3
"""Check the source collection, archived snapshots, and local preview links."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
import csv, hashlib, re, sys
root=Path(__file__).resolve().parents[2]
v=root/'vignettes'
assert len(list(v.glob('*.Rmd')))==14, 'Expected ten guides and four appendices'
assert not list((v/'articles').glob('*.Rmd')), 'Website-only tutorial sources remain'
for row in csv.DictReader((v/'archives/manifest.tsv').open(),delimiter='\t'):
    path=v/'archives'/row['archive_path']
    if path.exists():
        assert hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'],path
    else:
        assert path.suffix in ('.html','.gif'), 'Missing archived source: '+str(path)
for p in v.glob('*.Rmd'):
    text=p.read_text()
    assert '%\\VignetteEngine{knitr::rmarkdown}' in text,p
    assert not re.search(r'website-only article|depth_decay\s*=|branch_spread\s*=',text),p
class Page(HTMLParser):
    def __init__(self):super().__init__();self.links=[];self.ids=set()
    def handle_starttag(self,t,a):
        a=dict(a)
        if 'id' in a:self.ids.add(a['id'])
        if t=='a' and 'href' in a:self.links.append(a['href'])
out=Path(sys.argv[1]) if len(sys.argv)>1 else root/'output/vignette-previews'
cache={}
def page(p):
    p=p.resolve()
    if p not in cache:
        x=Page();x.feed(p.read_text());cache[p]=x
    return cache[p]
errors=[]
files=[out/'index.html',*sorted((out/'doc').glob('*.html')),*sorted((out/'articles').glob('*.html'))]
if not (out/'doc').exists():files += sorted((out/'reference').glob('*.html'))
for f in files:
    for href in page(f).links:
        u=urlsplit(href)
        if u.scheme or u.netloc or u.path.startswith('/'):continue
        target=f.parent/unquote(u.path) if u.path else f
        if not target.exists():errors.append((f.name,href,'missing file'))
        elif u.fragment and target.suffix=='.html' and unquote(u.fragment) not in page(target).ids:
            errors.append((f.name,href,'missing anchor'))
if errors:
    for e in errors:print(*e,sep=': ')
    raise SystemExit(1)
print('Verified 14 installed sources, archived checksums, and documentation links/anchors.')
