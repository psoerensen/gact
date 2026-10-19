"""Preserve this timing study without duplicating the genotype/large LD files."""
import csv
import hashlib
import sys
import zipfile
from pathlib import Path

root = Path.cwd().resolve()
evidence, report = (Path(arg).resolve() for arg in sys.argv[1:])
assert evidence.is_relative_to(root) and report.is_relative_to(root)
archive = report / "source-and-evidence.zip"
assert not archive.exists(), "Refusing to replace existing evidence"
files = {}

def include(path, prefix):
    if path.is_dir():
        for file in sorted(path.rglob("*")):
            if file.is_file():
                files[f"{prefix}/{file.relative_to(path).as_posix()}"] = file
    else:
        assert path.is_file(), path
        files[prefix] = path

include(evidence, "provenance")
for file in sorted(report.iterdir()):
    if file.is_file():
        include(file, f"report/{file.name}")
with (evidence / "runs.csv").open(encoding="utf-8-sig", newline="") as stream:
    for row in csv.DictReader(stream):
        run = Path(row["run"]).resolve()
        assert run.is_relative_to(root / "build" / "benchmarks")
        include(run, f"runs/{run.name}")
for package in ("gbase", "gmap", "gbayes"):
    include(root / "build" / "task-packages" / "library" / package, f"installed/{package}")
for script in ("run_manuscript.ps1", "simulated_finemapping_comparison_stage.R",
               "time_ld_swap.ps1", "summarize_ld_swap.R", "audit_ld_swap_example.R",
               "archive_ld_swap.py"):
    include(root / "benchmarks" / "genomic_20k_50k" / script, f"scripts/{script}")
include(root / "build/benchmarks/manuscript-20261005-142319/website_sim_stats-1/output.rds",
        "fixtures/frozen-statistics.rds")
include(root / "build/benchmarks/manuscript-20261005-231349/cmp_suite_prepare-1/human-validated-map.bim",
        "fixtures/human-validated-map.bim")
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
    for name, path in files.items():
        bundle.write(path, name)
with zipfile.ZipFile(archive) as bundle:
    assert bundle.testzip() is None, "ZIP integrity check failed"
digest = hashlib.file_digest(archive.open("rb"), "sha256").hexdigest()
(report / "source-and-evidence.sha256").write_text(f"{digest}  {archive.name}\n", encoding="utf-8")
print(f"Archive: {archive}\nBytes: {archive.stat().st_size}\nEntries: {len(files)}\nSHA256: {digest}")
