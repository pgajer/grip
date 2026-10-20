"""Dependency-aware ten-worker stage comparison; 50 GiB process-tree RSS watchdog."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time

STAGES = ("E11", "E12", "E14", "E16", "E18", "P0", "P2", "P5", "P8", "P14")
DIMS = (50, 25, 20, 10, 5, 4, 3)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for data in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    for name in ("root", "input-dir", "baseline", "paused-run", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    for key in ("R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE"):
        env.pop(key, None)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
        env[key] = "1"
        os.environ[key] = "1"
    import numpy as np
    script = Path(__file__).with_suffix(".R")
    limit = 50 * 1024 ** 3
    began = time.time()
    peak = 0
    all_jobs = []
    active = []
    done = []
    resources_handle = open(args.out / "resources.csv", "w", buffering=1)
    resources = csv.writer(resources_handle)
    resources.writerow(["unix_time", "rss_bytes", "cpu_percent", "active_workers"])
    pause_before = {f: sha256(args.paused_run / f) for f in ("status.json", "latest.json")}

    def publish(status, **kw):
        data = dict(status=status, updated=time.time(), started=began, pid=os.getpid(),
                    peak_rss_bytes=peak, rss_limit_bytes=limit, planned_fits=70,
                    completed_fits=sum(t[0] == "fit" for t in done),
                    completed_preparations=sum(t[0] == "prepare" for t in done),
                    active=[j[0] for j in active], **kw)
        temp = args.out / "status.json.tmp"
        temp.write_text(json.dumps(data, indent=2) + "\n")
        temp.replace(args.out / "status.json")

    def sample():
        nonlocal peak
        lines = subprocess.check_output(["ps", "-axo", "pid=,ppid=,rss=,%cpu="], text=True)
        table = {}
        for line in lines.splitlines():
            p, parent, rss, cpu = line.split()
            table[int(p)] = (int(parent), int(rss) * 1024, float(cpu))
        family = {os.getpid()}
        while True:
            extra = {p for p, row in table.items() if row[0] in family}
            if extra <= family:
                break
            family |= extra
        rss = sum(table[p][1] for p in family if p in table)
        cpu = sum(table[p][2] for p in family if p in table)
        peak = max(peak, rss)
        resources.writerow([time.time(), rss, cpu, len(active)])
        if rss > limit:
            raise MemoryError("Aggregate RSS exceeded 50 GiB; completed checkpoints retained")

    def launch(task):
        mode, stage, dim = task
        handle = open(args.out / stage / (mode + ("-" + str(dim) if dim else "") + ".log"), "w")
        command = ["Rscript", str(script), mode, str(args.root), str(args.out / stage)]
        if dim:
            command.append(str(dim))
        process = subprocess.Popen(command, env=env, stdout=handle,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        job = (task, process, handle)
        active.append(job); all_jobs.append(job)

    try:
        assert json.loads((args.paused_run / "status.json").read_text())["status"] == "paused_by_user"
        publish("validating_input")
        manifest = list(csv.DictReader(open(args.baseline / "source-manifest.csv")))
        assert all(sha256(row["path"]) == row["sha256"] for row in manifest)
        source_rows = [row for row in manifest if "/vignettes/" not in row["path"]]
        source_rows += [dict(path=str(p), sha256=sha256(p)) for p in (script, Path(__file__))]
        with open(args.out / "source-manifest.csv", "w") as handle:
            writer = csv.DictWriter(handle, fieldnames=("path", "sha256"))
            writer.writeheader(); writer.writerows(source_rows)
        projection = json.loads((args.input_dir / "projection.json").read_text())
        source = args.input_dir / "pc1814.npy"
        assert sha256(source) == projection["projection_sha256"]
        annotations = list(csv.DictReader(open(args.input_dir / "annotations.csv")))
        assert len(annotations) == 120804 and len({r["cell"] for r in annotations}) == 120804
        assert {r["stage"] for r in annotations} == set(STAGES)
        array = np.load(source, mmap_mode="r")
        assert array.shape == (120804, 1814)
        for j, stage in enumerate(STAGES):
            stage_dir = args.out / stage
            stage_dir.mkdir()
            seed = 2026102001 + j
            pool = np.array([i for i, row in enumerate(annotations) if row["stage"] == stage])
            assert len(pool) >= 1500
            indices = np.sort(np.random.default_rng(seed).choice(pool, 1500, replace=False))
            X = np.asarray(array[indices], dtype="<f8", order="C")
            assert np.isfinite(X).all()
            X.tofile(stage_dir / "input.bin")
            with open(stage_dir / "cells.csv", "w") as handle:
                writer = csv.DictWriter(handle, fieldnames=("index", *annotations[0].keys()))
                writer.writeheader()
                writer.writerows(dict(index=int(i) + 1, **annotations[i]) for i in indices)
            config = dict(stage=stage, available=len(pool), n=1500, input_dimension=1814,
                          sampling_seed=seed, sampling_rng="numpy.default_rng PCG64",
                          preparation_seed=20261101 + j, fit_seed=20190619,
                          input_sha256=sha256(stage_dir / "input.bin"),
                          cells_sha256=sha256(stage_dir / "cells.csv"),
                          projection_sha256=projection["projection_sha256"])
            (stage_dir / "input.json").write_text(json.dumps(config, indent=2) + "\n")
            sample()
        del array, X, annotations
        (args.out / "provenance.json").write_text(json.dumps(dict(
            projection=projection, annotations_sha256=sha256(args.input_dir / "annotations.csv"),
            baseline=str(args.baseline), paused_run=str(args.paused_run), pause_before=pause_before,
            dimensions=DIMS, stages=STAGES, maximum_workers=10, threads_per_worker=1,
            native_workspace_bytes=4 * 1024 ** 3, rss_limit_bytes=limit), indent=2) + "\n")
        pending = [("prepare", stage, None) for stage in STAGES]
        while pending or active:
            for job in list(active):
                task, process, handle = job
                code = process.poll()
                if code is None:
                    continue
                handle.close(); active.remove(job)
                if code:
                    raise RuntimeError("Task failed: " + str(task) + "; see its log")
                done.append(task)
                if task[0] == "prepare":
                    pending.extend(("fit", task[1], d) for d in DIMS)
            # Preparations unlock work; among ready fits, larger dimensions start first.
            pending.sort(key=lambda t: (0 if t[0] == "prepare" else 1, -(t[2] or 0), STAGES.index(t[1])))
            while pending and len(active) < 10:
                launch(pending.pop(0))
            sample(); publish("running")
            if active:
                time.sleep(1)
        assert len(done) == 80
        assert all(sha256(row["path"]) == row["sha256"] for row in source_rows)
        assert pause_before == {f: sha256(args.paused_run / f) for f in pause_before}
        publish("completed", full_run_unchanged=True, elapsed_seconds=time.time() - began)
    except BaseException as error:
        for _, process, _ in active:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
        for _, process, _ in active:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL); process.wait()
        publish("failed", message=str(error))
        raise
    finally:
        for _, _, handle in all_jobs:
            handle.close()
        resources_handle.close()


if __name__ == "__main__":
    main()
