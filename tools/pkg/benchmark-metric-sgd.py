#!/usr/bin/env python3
"""Bounded paired comparison; standalone Python standard library supervisor.

Run after independent audit against an installed candidate, for example:
GRIP_BENCHMARK_LIBRARY=<candidate-library> python3 tools/pkg/benchmark-metric-sgd.py <output> <candidate.tar.gz>
Output directories must be fresh. This never edits reference fixtures.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
WORKER = Path(__file__).with_suffix(".R")
OUT = Path(sys.argv[1]).resolve()
ARTIFACT = Path(sys.argv[2]).resolve()
if not os.environ.get("GRIP_BENCHMARK_LIBRARY"):
    raise SystemExit("Set GRIP_BENCHMARK_LIBRARY to the installed candidate library")
if not ARTIFACT.is_file():
    raise SystemExit("Supply the source tarball installed in the candidate library")
if OUT.exists():
    raise SystemExit("Use a fresh output directory; original runs must be retained")
OUT.mkdir(parents=True)
protocol = dict(version=1, sizes=[64,256], application="bundled karate-club graph (34 vertices)",
    dimensions=[2,3], seeds=[1,2,3,4,5], budgets=[10,30,100],
    backends=["smacof","sgd"], per_process_seconds=30, per_process_rss_bytes=1024**3,
    campaign_seconds=1200, rss_poll_seconds=.05, execution="sequential alternating backend order",
    defaults="provisional calibration choices; no parameter tuning in this pilot",
    elapsed_thresholds_seconds=[.01,.03,.1,.3,1,3],
    timing="fit_seconds includes wrapper, normalization and checkpoint scoring; excludes process startup, common preloading of SMACOF in both workers, and frozen distance preparation",
    memory="sampled whole-worker RSS, including R and loaded dependencies; not an exact peak",
    comparison="retrospective best completed-fit envelope over three separate budgets, not a single anytime run; mark missing candidates unavailable; also report summed effort over all three fits",
    r_worker_sha256=hashlib.sha256(WORKER.read_bytes()).hexdigest(),
    supervisor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    candidate_tarball_sha256=hashlib.sha256(ARTIFACT.read_bytes()).hexdigest(),
    candidate_tarball_name=ARTIFACT.name,
    source_git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
    candidate_source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
        for p in ["R/metric_mds.R","R/metric_mds_sgd.R","src/metric_mds_sgd.cpp","DESCRIPTION"]})
(OUT/"protocol.json").write_text(json.dumps(protocol,indent=2)+"\n")
with (OUT/"prepare.log").open("w") as log:
    try:
        subprocess.run(["Rscript","--vanilla",str(WORKER),"prepare",str(OUT)],cwd=ROOT,
                       stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
    except (subprocess.CalledProcessError,subprocess.TimeoutExpired) as error:
        (OUT/"preparation-status.json").write_text(json.dumps(dict(status="preparation_failed",error=str(error)),indent=2)+"\n")
        raise SystemExit("Fixture preparation failed; no fits started. See prepare.log")
cases = list(csv.DictReader((OUT/"cases.csv").open()))
hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob("*.rds")}
(OUT/"fixture-sha256.json").write_text(json.dumps(hashes,indent=2)+"\n")
jobs=[]
for index, case in enumerate(cases):
    for budget in protocol["budgets"]:
        backends=protocol["backends"] if index%2==0 else protocol["backends"][::-1]
        for backend in backends:
            jobs.append(dict(**case,backend=backend,budget=budget))
(OUT/"jobs.json").write_text(json.dumps(jobs,indent=2)+"\n")
fields=list(jobs[0])+["status","fit_seconds","rmse","raw_stress","initial_raw_stress",
    "iterations","pair_updates","termination","warning","error","process_seconds","sampled_peak_rss_bytes"]
campaign_start=time.monotonic()
with (OUT/"results.csv").open("w",newline="") as stream:
    writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); stream.flush()
    for index, job in enumerate(jobs):
        row=dict(job)
        if time.monotonic()-campaign_start > protocol["campaign_seconds"]:
            row.update(status="not_run_campaign_budget",error="Campaign time allowance exhausted")
            writer.writerow(row); stream.flush(); continue
        stem=f'{job["case"]}-{job["backend"]}-{job["budget"]}'
        destination=OUT/(stem+".csv")
        with (OUT/(stem+".log")).open("w") as log:
            started=time.monotonic(); peak=0; status=None
            process=subprocess.Popen(["Rscript","--vanilla",str(WORKER),"fit",str(OUT),job["case"],
                job["backend"],str(job["budget"]),str(destination)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            while process.poll() is None:
                observed=subprocess.run(["ps","-o","rss=","-p",str(process.pid)],capture_output=True,text=True)
                if observed.stdout.strip():
                    peak=max(peak,int(observed.stdout.strip())*1024)
                if peak > protocol["per_process_rss_bytes"]:
                    status="memory_limit"
                elif time.monotonic()-started > protocol["per_process_seconds"]:
                    status="time_limit"
                elif time.monotonic()-campaign_start > protocol["campaign_seconds"]:
                    status="campaign_time_limit"
                if status:
                    process.terminate()
                    try: process.wait(timeout=2)
                    except subprocess.TimeoutExpired: process.kill(); process.wait()
                    break
                time.sleep(protocol["rss_poll_seconds"])
            if destination.exists():
                try:
                    with destination.open() as result_file:
                        payload=next(csv.DictReader(result_file))
                    if not payload.get("status"): raise ValueError("Result lacks status")
                    row.update(payload)
                except (OSError,ValueError,StopIteration,csv.Error) as error:
                    row.update(status="invalid_result",error=str(error))
            if status or process.returncode != 0:
                row.update(status=status or "process_error",error=f"worker exit {process.returncode}; see {stem}.log")
            elif not destination.exists(): row.update(status="missing_result",error="Worker produced no result")
            row.update(process_seconds=time.monotonic()-started,sampled_peak_rss_bytes=peak)
        writer.writerow(row); stream.flush()
        if (index+1)%20==0: print(f"Recorded {index+1}/{len(jobs)} fits",flush=True)
print(f"Completed manifest: {OUT/'results.csv'}",flush=True)
