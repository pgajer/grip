#!/usr/bin/env python3
"""Restore fitting CSV inputs from tracked compact data; no geodesic recomputation."""
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;OUT=HERE.parents[2]/'output/mds-edge-kk-radius/inputs';OUT.mkdir(parents=True,exist_ok=True)
for f in sorted((HERE/'summary/inputs').glob('*.npz')):
 c=np.load(f);np.savetxt(OUT/(f.stem+'-distance.csv'),c['d'],delimiter=',',fmt='%.17g');np.savetxt(OUT/(f.stem+'-truth.csv'),c['truth'],delimiter=',',fmt='%.17g')
print('Fitting inputs restored')
