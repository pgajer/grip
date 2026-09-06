#!/usr/bin/env python3
"""Check citation keys and the associated source-support verification records."""
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
import re

P=Path(__file__).resolve().parents[1]
tex=(P/'geodesic_mds.tex').read_text()
keys=set()
for block in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\])*\{([^}]+)\}',tex):
    keys.update(x.strip() for x in block.split(','))
bib=set(re.findall(r'@\w+\s*\{\s*([^,]+),',(P/'geodesic_mds.bib').read_text()))
class Evidence(HTMLParser):
    def __init__(self):
        super().__init__();self.rows=[];self.current=None
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if 'data-citation-key' in a:
            self.current=dict(key=a['data-citation-key'],status=a.get('data-status'),links=[])
            self.rows.append(self.current)
        if self.current is not None and tag=='a' and 'data-source-link' in a and a.get('href'):
            self.current['links'].append(a['href'])
    def handle_endtag(self,tag):
        if tag=='tr':self.current=None
e=Evidence();e.feed((P/'citation_verification.html').read_text())
assert keys and keys<=bib,(keys,bib)
counts=Counter(r['key'] for r in e.rows)
assert set(counts)==keys and all(n==1 for n in counts.values()),counts
assert all(r['status']=='verified' and r['links'] for r in e.rows)
log=P/'build/geodesic_mds.log'
if log.exists():
    text=log.read_text(errors='replace')
    assert not re.search(r'(Citation .*undefined|There were undefined citations)',text)
print(f'Citation records passed for {len(keys)} cited sources; substantive support is recorded in HTML.')
