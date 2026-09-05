#!/usr/bin/env python3
"""Narrow citation contract check, with no private-tree dependency."""
from pathlib import Path
from html.parser import HTMLParser
import re
p=Path(__file__).resolve().parent
keys=set(k.strip() for group in re.findall(r'\\cite\w*\{([^}]+)\}',(p/'report.tex').read_text()) for k in group.split(','))
bib=set(re.findall(r'@\w+\s*\{\s*([^,]+)',(p/'references.bib').read_text()))
class Evidence(HTMLParser):
    def __init__(self):super().__init__();self.rows={};self.current=None
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if 'data-citation-key' in a:
            key=a['data-citation-key'];assert key not in self.rows
            self.rows[key]={'status':a.get('data-status'),'source':False};self.current=key
        if tag=='a' and self.current and 'data-source-link' in a and a.get('href'):
            self.rows[self.current]['source']=True
    def handle_endtag(self,tag):
        if tag=='tr':self.current=None
h=Evidence();h.feed((p/'citation_verification.html').read_text())
assert keys==set(h.rows) and keys<=bib
assert all(v['status']=='verified' and v['source'] for v in h.rows.values())
log=p/'build'/'report.log'
if log.exists():assert not re.search(r'Citation .* undefined|There were undefined references',log.read_text())
print('Citation gate passed:',', '.join(sorted(keys)))
