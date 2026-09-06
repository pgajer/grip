#!/usr/bin/env python3
"""Read-only verification of portable evidence bytes and design coverage."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
EVIDENCE = PAPER / "evidence"


def table(study, name):
    with gzip.open(EVIDENCE / study / "summary" / (name + ".csv.gz"), "rt") as stream:
        return list(csv.DictReader(stream))


def main():
    manifest = json.loads((EVIDENCE / "source-manifest.json").read_text())
    for record in manifest["files"]:
        payload = (EVIDENCE / record["export_path"]).read_bytes()
        assert len(payload) == record["export_bytes"]
        assert hashlib.sha256(payload).hexdigest() == record["export_sha256"]
        raw = gzip.decompress(payload) if record["export_path"].endswith(".gz") else payload
        assert len(raw) == record["source_bytes"]
        assert hashlib.sha256(raw).hexdigest() == record["source_sha256"]
    expected = {"study-initializers": {"scores": 125, "starts": 150},
                "study-radius": {"scores": 16416, "graphs": 1216, "starts": 3648,
                                 "optimizer-sensitivity": 512, "objective-validation": 176,
                                 "coordinate-validation": 1}}
    for study, tables in expected.items():
        for name, rows in tables.items():
            actual = len(table(study, name))
            assert actual == rows, (study, name, actual, rows)
    coordinate_check = table("study-radius", "coordinate-validation")[0]
    assert int(coordinate_check["candidates"]) == 16928
    assert float(coordinate_check["max_scaled_discrepancy"]) < 1e-10
    routes = EVIDENCE/"retained-route-validation"
    route_manifest = json.loads((routes/"source-manifest.json").read_text())
    assert route_manifest["configurations"] == 24 and route_manifest["routes_per_configuration"] == 28680
    assert hashlib.sha256((routes/"independent-check.csv").read_bytes()).hexdigest() == route_manifest["output_sha256"]
    print(f"Verified {len(manifest['files'])} evidence files and all declared table counts.")


if __name__ == "__main__":
    main()
