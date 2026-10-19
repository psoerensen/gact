"""Summarize completed candidates and preserve capped-run status separately."""
import csv
import json
from pathlib import Path
import re
import sys

root = Path(sys.argv[1])
records, candidates = [], []
for resource in sorted(root.glob("*/resources.json")):
    r = json.loads(resource.read_text())
    directory = resource.parent
    text = (directory / "stderr.log").read_text()
    phases = re.findall(r"phase=(\S+) elapsed=([\d.eE+-]+)", text)
    values = {}
    metrics = directory / "fit.metrics.tsv"
    if metrics.exists():
        for line in metrics.read_text().splitlines()[1:]:
            key, value = line.split("\t", 1)
            values[key] = value
    row = {k: r[k] for k in ("label", "backend", "timeout", "exit_code",
        "wall_seconds", "cpu_seconds", "peak_resident_mib", "sampled_private_mib")}
    row.update(chunk_size=int(r["command"][9]),
               gram_cache_mib=int(r["command"][8]),
               threads=int(r["command"][4]))
    row.update(last_completed_phase=phases[-1][0] if phases else "none",
        fit_status=values.get("fit_status", "not returned"),
        global_fit_seconds=values.get("global_fit_seconds", ""),
        complete_fit_seconds=values.get("fit_seconds", ""),
        calibration_seconds=values.get("calibration_seconds", ""),
        scan_seconds=values.get("scan_seconds", ""))
    records.append(row)
    for line in text.splitlines():
        if "candidate_log_delta=" not in line:
            continue
        fields = dict(re.findall(r"(\w+)=(\S+)", line))
        fields["label"] = r["label"]
        candidates.append(fields)
for name, rows in (("probe-summary.csv", records), ("completed-candidates.csv", candidates)):
    if rows:
        with (root / name).open("w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
print(json.dumps(records, indent=2))
findings = {
    "scope": "Bounded native diagnosis of the standalone glma infinitesimal LOCO path",
    "frozen_input": {"individuals": 5000, "markers": 37991,
                     "chromosomes": 22, "genotype_blocks_copied_verbatim": True},
    "small_input": {"individuals": 1000, "markers": 10000,
                    "chromosomes": 3, "fixture": "build/glma/scaling/diagnostic"},
    "controls": {"simulation_pairs": 15, "ridge_responses_per_candidate": 16,
                 "relative_tolerance": 1e-7, "chunk_size": "256, with one 128-width diagnostic",
                 "gram_cache_mib": 64, "loco_batch_size": 3,
                 "threads": 1},
    "limitations": [
        "One sequential probe per configuration on a shared laptop; no precise speed ratio or large-data throughput claim.",
        "Timeouts are incomplete fits, not successful timings.",
        "CLI sweep cap is 1000; R fitting cap is 100. Common completed candidate traces are comparable only below 100 sweeps.",
        "The earlier 88-minute R run had no phase observers, so its exact stopping phase cannot be reconstructed.",
        "Peak resident/private memory is measured; workflow budget is an admission estimate.",
        "No numerical implementation or statistical default was changed."
    ],
    "source_findings": [
        {"finding": "Portable standalone profile disables MKL. Chunk ridge rejects threads>1 without MKL and uses serial Eigen products at one thread.",
         "owners": ["gsuite", "gsolve"],
         "sources": ["packages/glma/src/Makevars.win:31", "../gsolve/src/chunk_ridge.cpp:62"]},
        {"finding": "R algorithm selects the projection route; global and LOCO variance estimation always call the chunk-ridge variance-ratio fitter.",
         "owners": ["gsuite", "glma", "greml"],
         "sources": ["packages/glma/src/glma.cpp:354", "packages/glma/src/glma.cpp:414", "../glma/src/mixed_scan.cpp:1757"]},
        {"finding": "Each candidate solves phenotype plus 15 simulations. Every chromosome has its own scalar refit; three-model batching shares traversals but retains each model's arithmetic and convergence.",
         "owners": ["glma", "greml", "gsolve"],
         "sources": ["../greml/src/association_variance_ratio.cpp:66", "../glma/src/mixed_scan.cpp:1849"]},
        {"finding": "LOCO sessions inherit the global variance ratio, but start with empty coefficients. Calibration and scan separately construct chromosome projection states and solve the phenotype without warm coefficients.",
         "owners": ["glma", "greml"],
         "sources": ["../greml/src/association_variance_ratio.cpp:49", "../glma/src/mixed_scan.cpp:1917", "../glma/src/mixed_scan.cpp:2004", "../glma/src/mixed_scan.cpp:2025"]},
        {"finding": "The R block width 256 requires 74.08 MiB of within-chunk Grams for 37991 markers, exceeding the fixed 64 MiB cache. The last 5223 markers rebuild Grams and use coordinate fallback on each sweep. Width 128 requires about 37.08 MiB. The capped comparison is a diagnostic, not a completed workflow timing.",
         "owners": ["gsuite", "gsolve"],
         "sources": ["packages/glma/R/glma.R:272", "../gsolve/src/chunk_ridge.cpp:217", "../gsolve/src/chunk_ridge.cpp:334"]}
    ],
    "recommended_order": [
        "Expose optional native phase/candidate diagnostics through R, including sweeps, backend and cache counts.",
        "Measure a solver chunk width/cache budget that retains all within-chunk Grams, separately from scan width, without changing tolerances or the model.",
        "Reuse accepted LOCO phenotype projections and plans across calibration and scanning under a bounded memory budget; compare all association statistics against the current route.",
        "Add a supported threaded BLAS route for portable standalone packages in gsolve, or qualify an explicitly optional MKL profile.",
        "Measure global/LOCO variance fitting after these changes before selecting a different fitting controller; faster approximate prediction fitting has different acceptance semantics."
    ]
}
def gram_bytes(markers, chunk):
    blocks, tail = divmod(markers, chunk)
    return 8 * (blocks * chunk * chunk + tail * tail)
findings["gram_calculations"] = {
    "width_256_bytes": gram_bytes(37991, 256),
    "width_128_bytes": gram_bytes(37991, 128),
    "uncached_markers_at_width_256": 37991 - (64 * 2**20 // (8 * 256**2)) * 256
}
(root / "findings.json").write_text(json.dumps(findings, indent=2))
reference = [r for r in candidates if r["label"] == "frozen-mkl"]
chunk128 = [r for r in candidates if r["label"] == "frozen-mkl-chunk128"]
if reference and chunk128:
    for key in ("candidate_log_delta", "equation", "sweeps"):
        assert float(reference[0][key]) == float(chunk128[0][key])
    findings["first_candidate_cache_comparison"] = {
        "width256_solve_seconds": float(reference[0]["solve_s"]),
        "width128_solve_seconds": float(chunk128[0]["solve_s"]),
        "width256_product_seconds": float(reference[0]["product_s"]),
        "width128_product_seconds": float(chunk128[0]["product_s"]),
        "candidate_log_delta": float(reference[0]["candidate_log_delta"]),
        "equation": float(reference[0]["equation"]),
        "sweeps": int(reference[0]["sweeps"]),
        "agreement": "Recorded first candidate agrees at printed precision; incomplete full fits."
    }
    (root / "findings.json").write_text(json.dumps(findings, indent=2))
# Independent route check: matching log-ratios/equations on common completed
# candidates establishes numerical agreement, not whole-workflow equivalence.
labels = sorted(set(r["label"] for r in candidates))
if "small-portable" in labels and "small-mkl" in labels:
    a = [r for r in candidates if r["label"] == "small-portable"]
    b = [r for r in candidates if r["label"] == "small-mkl"]
    for x, y in zip(a, b):
        assert x["phase"] == y["phase"]
        for key in ("candidate_log_delta", "equation", "sweeps"):
            assert float(x[key]) == float(y[key]), (key, x, y)
    print(f"Common completed candidate traces agree at printed precision: {min(len(a),len(b))}")
