"""Ten single-threaded fits with a sampled aggregate 50 GiB RSS watchdog."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--prior-protocol", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
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
    peak = 0
    jobs = []
    resource_file = open(args.out / "resources.csv", "w", buffering=1)
    resources = csv.writer(resource_file)
    resources.writerow(["unix_time", "stage", "total_rss_bytes", "total_cpu_percent", "processes"])

    def publish(state, **kwargs):
        data = dict(status=state, updated=time.time(), pid=os.getpid(),
                    peak_rss_bytes=peak, rss_limit_bytes=limit, **kwargs)
        temp = args.out / "status.json.tmp"
        temp.write_text(json.dumps(data, indent=2) + "\n")
        temp.replace(args.out / "status.json")

    def sample(stage):
        nonlocal peak
        lines = subprocess.check_output(["ps", "-axo", "pid=,ppid=,rss=,%cpu="], text=True)
        table = {}
        for line in lines.splitlines():
            fields = line.split()
            if len(fields) == 4:
                table[int(fields[0])] = (int(fields[1]), int(fields[2]) * 1024, float(fields[3]))
        descendants = {os.getpid()}
        while True:
            extra = {pid for pid, row in table.items() if row[0] in descendants}
            if extra <= descendants:
                break
            descendants |= extra
        rss = sum(table[pid][1] for pid in descendants if pid in table)
        cpu = sum(table[pid][2] for pid in descendants if pid in table)
        peak = max(peak, rss)
        resources.writerow([time.time(), stage, rss, cpu, len(descendants)])
        if rss > limit:
            raise MemoryError("Combined process-tree RSS exceeded 50 GiB; completed checkpoints retained")

    def launch(mode, extra=()):
        name = "-".join([mode] + list(map(str, extra))) if mode != "export" else mode
        handle = open(args.out / (name + ".log"), "w")
        process = subprocess.Popen(["Rscript", str(script), mode, str(args.root),
                                    str(args.out), *map(str, extra)], env=env,
                                   stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
        jobs.append((name, process, handle))
        return process

    def await_jobs(stage):
        while any(p.poll() is None for _, p, _ in jobs):
            sample(stage)
            failed = [name for name, p, _ in jobs if p.poll() not in (None, 0)]
            if failed:
                raise RuntimeError("Failed subprocess: " + ", ".join(failed))
            publish(stage, active=[name for name, p, _ in jobs if p.poll() is None])
            time.sleep(2)
        sample(stage)
        failed = [name for name, p, _ in jobs if p.returncode != 0]
        if failed:
            raise RuntimeError("Failed subprocess: " + ", ".join(failed))

    try:
        publish("exporting_subset")
        launch("export", [args.prior_protocol]); await_jobs("exporting_subset")
        subset = list(csv.DictReader(open(args.out / "subset.csv")))
        annotations = list(csv.DictReader(open(args.input_dir / "annotations.csv")))
        indices = np.array([int(row["index"]) - 1 for row in subset])
        assert all(row["cell"] == annotations[i]["cell"] for row, i in zip(subset, indices))
        manifest = json.loads((args.input_dir / "projection.json").read_text())
        source = args.input_dir / "pc1814.npy"
        assert sha256(source) == manifest["projection_sha256"]
        array = np.load(source, mmap_mode="r")
        assert array.shape == (120804, 1814) and len(indices) == 1500
        X = np.asarray(array[indices], dtype="<f8", order="C")
        assert np.isfinite(X).all()
        X.tofile(args.out / "subset-pc1814.bin")
        (args.out / "input-provenance.json").write_text(json.dumps(dict(
            projection=manifest, source=str(source), prior_protocol=str(args.prior_protocol),
            prior_protocol_sha256=sha256(args.prior_protocol),
            subset_sha256=sha256(args.out / "subset-pc1814.bin"),
            cell_order_verified=True), indent=2) + "\n")
        del array, X, annotations
        launch("prepare"); await_jobs("preparing_shared_constraints")
        for dim in (300, 200, 100, 50, 25, 20, 10, 5, 4, 3):
            launch("fit", [dim])
        await_jobs("fitting")
        publish("completed", dimensions=[300, 200, 100, 50, 25, 20, 10, 5, 4, 3], epochs=100)
    except BaseException as error:
        # Each R worker owns a process group; never signal unrelated retinal work.
        for _, process, _ in jobs:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
        for _, process, _ in jobs:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        publish("failed", message=str(error))
        raise
    finally:
        for _, _, handle in jobs:
            handle.close()
        resource_file.close()


if __name__ == "__main__":
    main()
