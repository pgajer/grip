#!/usr/bin/env python3
"""One-cloud adapter around the Figure 7 numerical-reference implementation."""
import argparse
import csv
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from importlib.metadata import version

spec = importlib.util.spec_from_file_location(
    "pilot_reference", Path(__file__).with_name("surface-reference.py"))
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)


def run_reference(out, cloud, x, resolutions, sources):
    """The pilot's mesh/solver, with explicit sources instead of its n=1000 sentinel."""
    targets = np.arange(len(x), dtype=np.int32)
    for resolution in resolutions:
        dest = out / f"reference-r{cloud:02d}-m{resolution}-s{len(sources)}.npz"
        if dest.exists() and dest.with_suffix(".csv").exists():
            print(f"Existing {dest.name}", flush=True)
            continue
        started = time.perf_counter()
        vertices, faces = pilot.mesh(x, resolution)
        solver = pilot.geodesic.PyGeodesicAlgorithmExact(vertices, faces)
        values = np.empty((len(sources), len(x)))
        for row, source in enumerate(sources):
            values[row], _ = solver.geodesicDistances(
                np.array([source], dtype=np.int32), targets)
            if row == 0 or (row + 1) % 16 == 0:
                print(f"cloud={cloud} mesh={resolution} sources={row+1}/{len(sources)} "
                      f"elapsed={time.perf_counter()-started:.1f}s", flush=True)
        assert np.isfinite(values).all() and (values >= -1e-12).all()
        assert np.max(np.abs(values[np.arange(len(sources)), sources])) < 1e-12
        elapsed = time.perf_counter() - started
        np.savez_compressed(dest, distances=values, sources=sources,
                            mesh_resolution=resolution, mesh_vertices=len(vertices),
                            mesh_faces=len(faces), elapsed=elapsed)
        np.savetxt(dest.with_suffix(".csv"), np.column_stack((sources+1, values)),
                   delimiter=",", fmt="%.17g")
        print(f"SAVED {dest.name}: {elapsed:.2f}s", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--cloud", type=int, required=True)
    parser.add_argument("--grids", required=True)
    parser.add_argument("--sources", type=int, required=True)
    parser.add_argument("--fine-grid", type=int, required=True)
    parser.add_argument("--fine-sources", type=int, required=True)
    args = parser.parse_args()
    x = np.loadtxt(args.out / f"cloud-{args.cloud:02d}.csv", delimiter=",", skiprows=1)
    sources = np.atleast_1d(np.loadtxt(
        args.out / f"sources-{args.cloud:02d}.txt", dtype=np.int32))
    assert 2 <= args.fine_sources <= args.sources <= len(x)
    run_reference(args.out, args.cloud, x, list(map(int, args.grids.split(","))), sources[:args.sources])
    run_reference(args.out, args.cloud, x, [args.fine_grid], sources[:args.fine_sources])

    # Reuse the pilot's independent plane and smooth-surface controls, with n
    # taken from this cloud rather than the original five-cloud constants.
    count = min(40, len(x))
    vertices, faces = pilot.mesh(x[:count], 25, amplitude=0)
    solver = pilot.geodesic.PyGeodesicAlgorithmExact(vertices, faces)
    distance, _ = solver.geodesicDistances(
        np.array([0], dtype=np.int32), np.arange(count, dtype=np.int32))
    plane_error = float(np.max(np.abs(
        distance - np.linalg.norm(vertices[:count] - vertices[0], axis=1))))
    assert plane_error < 1e-9, plane_error
    rng = np.random.default_rng(4211000 + args.cloud)
    rows = []
    for source in sources[:args.sources]:
        target = int(rng.choice(np.delete(np.arange(len(x)), source)))
        distance, status, extent = pilot.smooth_distance(x[source, :2], x[target, :2])
        if status != 0 or extent > 1 + 1e-8:
            raise RuntimeError(f"Smooth control failed: {source}, {target}, {status}, {extent}")
        rows.append([source + 1, target + 1, distance, status, extent])
    with (args.out / "smooth-checks.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["i", "j", "distance", "status", "max_abs_parameter"])
        writer.writerows(rows)
    metadata = {name: version(name) for name in ("numpy", "scipy", "pygeodesic")}
    metadata.update(plane_max_absolute_error=plane_error, smooth_checks=len(rows))
    (args.out / "reference-environment.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Plane and {len(rows)} smooth-surface controls passed", flush=True)


if __name__ == "__main__":
    main()
