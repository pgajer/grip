#!/usr/bin/env python3
"""Keep private manuscript workspaces out of the public package tree."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
paths = subprocess.check_output(
    ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
    cwd=root,
).decode().split("\0")
private_roots = ("papers/", "dev/papers/", "dev/grip_paper/", "tools/experiments/",
                 "tools/benchmarks/", "tools/reports/", "tools/figures/")
# Missing files are unstaged removals during migration; a clean checkout has none.
bad = sorted({p for p in paths if p.startswith(private_roots) and (root / p).exists()})
if bad:
    print("Manuscript files must live outside the public package:", file=sys.stderr)
    print("\n".join(bad[:20]), file=sys.stderr)
    sys.exit(1)
print("No private manuscript workspaces in the public package tree.")
