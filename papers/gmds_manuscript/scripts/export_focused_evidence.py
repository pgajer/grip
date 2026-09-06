#!/usr/bin/env python3
"""Freeze verified grip evidence for a portable focused-manuscript build.

Import is explicit; normal paper rendering does not require a grip checkout.
CSV exports are byte-preserving deterministic gzip containers. Source hashes
refer to decompressed bytes; package provenance stays within each study.
"""
import argparse
import gzip
import hashlib
import json
import subprocess
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
STUDIES = {
    "surface-theory": {
        "root": "tools/reports/geodesic-mds",
        "files": ["report.tex", "references.bib", "citation_verification.html"],
    },
    "study-initializers": {
        "root": "papers/grip-software-paper/reproducibility/experiments/mds-initializer-sensitivity",
        "files": ["README.md", "summary/results.md", "summary/validation.json",
                  "summary/scores.csv", "summary/selected-scores-percent.csv",
                  "summary/graph-selection.csv", "summary/graph-statistics.csv",
                  "summary/starts.csv", "summary/surface-scores.csv",
                  "summary/additional-budget-scores.csv", "summary/manifest.json"],
    },
    "study-radius": {
        "root": "tools/experiments/mds-edge-kk-radius",
        "files": ["README.md", "PROTOCOL.md", "PILOT.md", "summary/RESULTS.md",
                  "summary/validation.json", "summary/scores.csv", "summary/graphs.csv",
                  "summary/starts.csv", "summary/initializer-comparison.csv",
                  "summary/size-comparison.csv", "summary/scale-by-radius.csv",
                  "summary/locality.csv", "summary/snapshot-coordinates.csv",
                  "summary/optimizer-sensitivity.csv", "summary/objective-validation.csv",
                  "summary/coordinate-validation.csv"],
    },
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grip-source", required=True, type=Path)
    args = parser.parse_args()
    repo = args.grip_source.resolve()
    commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    records = []
    for study, spec in STUDIES.items():
        for rel in spec["files"]:
            source_rel = str(Path(spec["root"]) / rel)
            source = repo / source_rel
            raw = source.read_bytes()
            # Never silently freeze uncommitted or untracked source modifications.
            committed = subprocess.check_output(["git", "-C", str(repo), "show", f"{commit}:{source_rel}"])
            if raw != committed:
                raise RuntimeError(f"Source differs from committed evidence: {source_rel}")
            target_rel = str(Path(study) / rel) + (".gz" if source.suffix == ".csv" else "")
            target = PAPER / "evidence" / target_rel
            target.parent.mkdir(parents=True, exist_ok=True)
            payload = gzip.compress(raw, mtime=0) if source.suffix == ".csv" else raw
            target.write_bytes(payload)
            records.append(dict(study=study, source_path=source_rel, export_path=target_rel,
                                source_bytes=len(raw), source_sha256=hashlib.sha256(raw).hexdigest(),
                                export_bytes=len(payload), export_sha256=hashlib.sha256(payload).hexdigest()))
    manifest = dict(repository="https://github.com/pgajer/grip", commit=commit,
                    note="Independent studies; do not pool settings or reused samples as independent replications.",
                    files=records)
    (PAPER / "evidence/source-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Exported {len(records)} committed evidence files from {commit}.")


if __name__ == "__main__":
    main()
