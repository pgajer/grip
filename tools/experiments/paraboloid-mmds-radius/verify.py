#!/usr/bin/env python3
"""Check the Markdown citation contract and the generated experiment manifest."""
import csv
import hashlib
import json
from html.parser import HTMLParser
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent


class Citations(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = {}
        self.active = None

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        if tag == "tr" and "data-citation-key" in attrs:
            self.active = attrs["data-citation-key"]
            assert self.active not in self.rows, "Duplicate citation row"
            self.rows[self.active] = dict(status=attrs["data-status"], linked=False)
        if self.active and tag == "a" and "data-source-link" in attrs:
            assert attrs.get("href", "").startswith("https://")
            self.rows[self.active]["linked"] = True

    def handle_endtag(self, tag):
        if tag == "tr":
            self.active = None


def main():
    used = set(re.findall(r"<!-- cite:([\w-]+) -->", (ROOT/"README.md").read_text()))
    bib = set(re.findall(r"@\w+\{([\w-]+),", (ROOT/"references.bib").read_text()))
    evidence = Citations()
    evidence.feed((ROOT/"citation_verification.html").read_text())
    assert used == bib == set(evidence.rows), "Citation coverage mismatch"
    assert all(a["status"] == "verified" and a["linked"] for a in evidence.rows.values())
    out = ROOT.parents[2]/"output/paraboloid-mmds-radius"
    manifest = json.loads((out/"manifest.json").read_text())
    assert hashlib.sha256((ROOT/"experiment.py").read_bytes()).hexdigest() == manifest["script_sha256"]
    for filename, checksum in manifest["files"].items():
        assert hashlib.sha256((out/filename).read_bytes()).hexdigest() == checksum, filename
    metrics = list(csv.DictReader((out/"metrics.csv").open()))
    runs = list(csv.DictReader((out/"optimizer_runs.csv").open()))
    assert len(metrics) == 100 and len(runs) == 120
    assert sum(a["selected"] == "True" for a in runs) == 20
    assert manifest["checks"]["selected_runs_success"]
    assert manifest["checks"]["max_3d_relative_recovery"] < 1e-9
    print("Verified: two supported citations, artifact hashes, 100 metric rows, 120 starts, 20 selected fits.")


if __name__ == "__main__":
    main()
