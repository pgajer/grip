#!/usr/bin/env python3
"""Export highlighted coordinates from pinned grip blobs; normal builds are portable."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[1]
COMMIT = "a43edfef1bbee70e825e790677abb35f12dfff53"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grip-source", type=Path, required=True)
    args = parser.parse_args()
    records, frames = [], []
    with tempfile.TemporaryDirectory() as temporary:
        tmp = Path(temporary)
        for surface in ["paraboloid", "saddle"]:
            for rep in range(1, 4):
                case = f"{surface}-disk-rep{rep}-n240-r64"
                blobs = {}
                for subdir, suffix in [("coordinates", ".rds"), ("inputs", ".npz")]:
                    path = f"tools/experiments/mds-edge-kk-radius/summary/{subdir}/{case}{suffix}"
                    raw = subprocess.check_output(["git", "-C", str(args.grip_source), "show", f"{COMMIT}:{path}"])
                    records.append(dict(source_path=path, sha256=hashlib.sha256(raw).hexdigest()))
                    blobs[suffix] = raw
                (tmp/"coordinates.rds").write_bytes(blobs[".rds"])
                script = '''args <- commandArgs(TRUE)
v <- readRDS(args[1]); rows <- list()
for (fit in v) {
  id <- fit$scores[1,]
  if (id$k != 32) next
  for (method in c("stress", "stress_fixed_primary")) {
    z <- fit$candidates[[method]]
    if (is.null(z)) next
    rows[[length(rows)+1L]] <- data.frame(case=id$case,regime=id$regime,method=method,
      vertex=seq_len(nrow(z)),x=z[,1],y=z[,2],z=z[,3])
  }
}
write.csv(do.call(rbind,rows),args[2],row.names=FALSE)
'''
                (tmp/"export.R").write_text(script)
                subprocess.run(["Rscript", str(tmp/"export.R"), str(tmp/"coordinates.rds"), str(tmp/"coordinates.csv")], check=True)
                frames.append(pd.read_csv(tmp/"coordinates.csv"))
                truth = np.load(io.BytesIO(blobs[".npz"]))["truth"]
                frame = pd.DataFrame(truth, columns=["x", "y", "z"])
                frame["case"], frame["regime"], frame["method"] = case, "reference", "original"
                frame["vertex"] = np.arange(1, 241)
                frames.append(frame)
    data = pd.concat(frames, ignore_index=True).sort_values(["case", "regime", "method", "vertex"])
    assert len(data) == 7200
    target = P/"evidence/coordinate-components"
    target.mkdir(exist_ok=True)
    raw = data.to_csv(index=False, float_format="%.17g").encode()
    payload = gzip.compress(raw, mtime=0)
    (target/"coordinates.csv.gz").write_bytes(payload)
    (target/"source-manifest.json").write_text(json.dumps(dict(
        repository="https://github.com/pgajer/grip", commit=COMMIT,
        selection="Both surfaces; three base-disk replicates; n240 r64 k32; both graph regimes; stress before and after fixed-scale density continuation, plus truth.",
        sources=records, export_sha256=hashlib.sha256(payload).hexdigest(),
        rows=len(data), export_script="scripts/export_coordinate_components.py"), indent=2)+"\n")
    print(f"Exported {len(data)} coordinate rows from {len(records)} pinned blobs.")


if __name__ == "__main__":
    main()
