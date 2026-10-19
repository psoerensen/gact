# Standalone glma LOCO investigation, 6 October 2026

Internal bounded measurements; these are not manuscript performance claims.
Both standalone joint workflows now complete on the full frozen gact example:
Gaussian/ridge-only infinitesimal analysis in 57.00 seconds and prior-selected
analysis in 124.12 seconds. The slow standalone infinitesimal variance-search
route has been replaced; its measurements below remain historical evidence.

## Earlier LDAK comparison

The preserved native joint full-panel prediction comparison used 5,000 samples,
50,000 markers, one thread and fixed marker power -0.5. Whole-process CPU
seconds were 25.98/28.23/28.75 for glma and 22.71/22.40/32.23 for LDAK-KVIK.
Their timing ranges overlapped. This compared Windows-native glma with
WSL/Linux KVIK and did not establish exact estimator or runtime equivalence.
Glma's workflow returned predictions, including global/LOCO refits, rather
than association calibration and scanning. See
[the native comparison](../../../glma/docs/prediction-background.md#joint-full-panel-fitting-and-external-comparison-17-september-2026).

The standalone `prior_selected_loco` wrapper selects that joint full-panel
fitting mode. The earlier infinitesimal route made separate scalar
variance-ratio searches for the global model and chromosome exclusions.
Its replacement uses the joint engine with one Gaussian/ridge prior and a
shared global variance estimate. These are different estimators and amounts
of computational work; no old/new estimator equivalence is asserted.

## Earlier variance-search changes and checks

The infinitesimal fit uses ridge chunks no wider than 128, keeping all
within-chunk Grams within the existing cache for this fixture. Accepted LOCO
phenotype projections can be reused for calibration/scanning within the
memory allowance; stricter projection tolerances retain the recomputation
route. Portable gsolve now supports bounded OpenMP dense products for
multi-thread chunk ridge, and prior initialization forwards the thread count.
Statistical defaults and acceptance criteria were retained.

Focused native association/mixed-model checks, the portable dense/ridge
oracles, standalone glma R tests and clean-load/DLL closure checks passed.
This is partial runtime verification; a complete source-package check was
not run. The isolated standalone library was refreshed; the normal user R
library was not installed.

## Completed and capped measurements

The small fixture retains its covariate and phenotype: 1,000 samples,
10,000 markers and three chromosomes. The full frozen example has 5,000
samples, 37,991 retained markers and 22 chromosomes. Each measurement used
a fresh R process; actual native DLL paths were recorded.

| Input and method | Threads | Analysis seconds | Process seconds | Peak resident MiB | Outcome |
| --- | ---: | ---: | ---: | ---: | --- |
| Small infinitesimal, preserved baseline | 1 | 51.19 | 52.78 | 189.26 | Completed |
| Small infinitesimal, updated | 1 | 30.94 | 32.49 | 177.70 | Completed |
| Small infinitesimal, updated | 4 | 32.43 | 34.32 | 178.66 | Completed |
| Full prior-selected, before strong-marker fix | 1 | 24.26 to error | 25.84 | 318.66 | Failed |
| Full prior-selected, before strong-marker fix | 4 | 23.28 to error | 25.04 | 318.98 | Failed |
| Full infinitesimal, variance-search implementation | 4 | Incomplete | 600.08 | 622.07 | Timed out |
| Full prior-selected, strong-marker fix | 4 | 124.12 | 129.18 | 391.79 | Completed |
| Full infinitesimal, Gaussian joint replacement | 4 | 57.00 | 59.74 | 386.11 | Completed |

The small before/after maximum absolute association-statistic difference was
1.9549393e-7; one/four-thread updated results differed by 3.5527137e-14.
Four threads did not improve small-fixture elapsed time. The completed small
comparison cannot establish a full-data speedup.

The two earlier full prior-selected attempts failed with
`blocked RHE high-effect chunk requires phenotype adjustment`. The native
method then rejected strong-effect chunks; the adjustment was missing at that
point. These failures are not completed fit timings. The full infinitesimal
run used updated DLLs, consumed 837.98 CPU seconds and stopped at its explicit
600-second elapsed limit. All timing processes ended. The earlier user-stopped
full baseline was not repeated and supplies no completed speedup denominator.

Evidence is retained under `build/glma/scaling/r-lmm-investigation-20261006/`
and in `build/benchmarks/final-summaries/`, with study prefixes
`manuscript-20261006-181205--`, `181309--`, `181343--`, `181420--`,
`181448--` and `181515--` (each uses the full date-prefixed study identifier).
The run index is `build/gact-revision-cleanup/timing-fix-runs.json`; the timing
and agreement helpers are beside this report. Measurement logs and provenance
remain in the retained summaries; obsolete run directories have been removed.
No commit or push was made at that measurement stage.

## Completed joint RHE/prior-selected fit

The strong-marker fix adjusts the starting RHE calculation only. It identifies
strong marginal marker directions, removes their contribution for a second
blocked RHE pass, and combines that contribution with the remaining background
heritability estimate. The final joint candidate/Monte Carlo fit, global/LOCO
refits and association scan retain the original phenotype and all requested
markers. Original-design Gram cache seeds remain reusable. Inputs without
detected strong markers retain the single-pass initializer.

One fresh four-thread run used the frozen 5,000-individual, 37,991-marker,
22-chromosome fixture, the existing prior grid and automatic scaling. It
completed with all 37,991 association rows having status `ok` and finite
p-values, with no warnings. The selected prior was
`elastic_p0.010000_f0.100000`. Automatic scaling classified this panel as
low-structure and used scaling factor 1; the timing does not include the
additional ridge-reference calibration required for strongly structured data.

| Phase | Seconds |
| --- | ---: |
| Initial variance preparation and RHE | 28.93 |
| Joint prior candidates and Monte Carlo fitting | 17.15 |
| Global and chromosome-excluded refits | 58.46 |
| Total native prediction fit | 104.56 |
| Association scan | 9.79 |
| Complete R analysis | 124.12 |
| Complete process | 129.18 |

The fit total includes the first three phases; these rows must not be summed
with that total. R analysis also includes genotype preparation, structure
assessment and wrapper work. Process CPU time was 163.64 seconds. Peak
resident memory was 391.79 MiB; sampled peak private committed memory was
1,189.91 MiB. These are different Windows memory measures.

The largest measured phase is the global/LOCO refit batch. This successful run
does not establish an exact speedup over the unfinished infinitesimal run:
the methods perform different fitting work. No further timing or test campaign
was added after this completed run.

Evidence is in
`build/benchmarks/final-summaries/`, with the study prefix
`manuscript-20261006-190318--`; this retains the process metrics, timing and
diagnostics. The run index and summary are
`build/gact-revision-cleanup/rhe-timing-runs.json` and
`build/gact-revision-cleanup/rhe-timing-summary.csv`. Loaded native paths were
recorded. The measured `libglma.dll` came from the refreshed isolated standalone
library and matched the portable native prefix, SHA-256
`01DD536ED2984763C518474D82ED90E655B1CB2F8E602B76BAD0C90BC68D743A`.
All owned build and timing helpers have finished. The normal user R library
has not been replaced; changes were uncommitted at that measurement stage.

## Gaussian-only replacement for standalone infinitesimal LOCO

`glma(..., method = "infinitesimal_loco")` and `glma_infinitesimal_loco()` now
dispatch exclusively to native `joint_full_panel_ridge`, introduced in native
glma 0.1.1 (ABI 0). The standalone C++ bridge no longer calls
`fit_association()` or the repeated variance-ratio controller. Retired
projection algorithms and variance-search controls produce explicit errors.
Low-level native interfaces remain for existing native consumers and historical
comparisons; neither standalone R LOCO method invokes that slow controller.

The replacement has one Gaussian/ridge candidate and fixed marker power -1,
giving equal prior variances on empirically standardized markers. Blocked RHE
initializes the shared three-point Monte Carlo heritability grid. Global
variance estimates are held fixed across global and chromosome-excluded refits;
retained marker prior variances are normalized by the active marker count.
Warm starts use the closest full-response grid fit and reuse the original
projected chunk Grams. The association calculation is the existing scaled
residual-offset test, with automatic structure handling. The
[public workflow guide](../../docs/glma-workflow.qmd#account-for-a-genomic-background-with-loco)
documents the grid, stopping rule, no-bracket fallback and variance-scale
calculation, plus the meaning of returned effects and standard errors.

One four-thread measurement used the same frozen 5,000-individual,
37,991-marker, 22-chromosome gact input. All 37,991 associations returned `ok`
with finite p-values and no warnings. The selected prior was `ridge`, candidate
count 1, marker power -1 and heritability .1251001. The grid bracketed a zero;
`heritability_grid_used_initial` was false. The estimated global genetic and
residual variances were .1438813220 and 1.0062480675; all 22 LOCO records report
that same pair. Automatic scaling classified this fixture as low-structure
and used scaling factor 1. Strongly structured panels can require additional
ridge-reference calibration; that work is absent from this timing.

| Gaussian joint phase | Seconds |
| --- | ---: |
| Initial variance preparation and RHE | 18.00 |
| Joint ridge/grid fitting | 7.67 |
| Global and chromosome-excluded refits | 21.26 |
| Total native prediction fit | 46.93 |
| Association scan | 3.12 |
| Complete R analysis | 57.00 |
| Complete process | 59.74 |

The fit total includes the first three phases. R analysis also includes
preparation, structure assessment and wrapper work. Process CPU time was
86.77 seconds. Peak resident memory was 386.11 MiB; sampled peak private
committed memory was 1,185.20 MiB. These Windows measures are distinct.
This is one bounded full-example measurement, not a repeated benchmark study
or evidence of large-data throughput. It supplies no exact speedup denominator
for the unfinished old variance-search workflow.

The existing focused native prediction check passed in 9.53 seconds, including
an independent dense Gaussian/ridge solution comparison for global and
chromosome-excluded effects. Standalone package checks passed, including
Gaussian-only dispatch, retired-control rejection, multi-trait preparation,
streamed background subsets, strong-marker association and thread agreement.
Clean-load and portable DLL closure checks passed. No complete source-package
check or broad null-tail qualification was run.

Fresh evidence is in
`build/benchmarks/final-summaries/`, with the study prefix
`manuscript-20261006-192735--`. The run index and
summary are `build/gact-revision-cleanup/gaussian-timing-runs.json` and
`build/gact-revision-cleanup/gaussian-timing-summary.csv`; build output is in
`gaussian-build2.stdout.log`/`gaussian-build2.stderr.log` beside them. Native
runtime version 0.1.1 and actual loaded paths were recorded. The measured
`libglma.dll` matched the installed portable prefix, SHA-256
`002448C2DD4A24C374C991E01DEBA4852DF085D0C44F8B3AF39F432370571235`.

The isolated standalone library and binary were refreshed. The normal user R
library remains untouched. The public documentation candidate rendered and
validated 7,442 internal references across 92 HTML pages with no failures.
At measurement completion, the required gtools updater refused to overwrite uncommitted managed
files, so the fresh candidate remains in
`build/gtools-public/candidate-1dso7_8t/repository/` and gtools was preserved.
All owned build, documentation and timing processes had finished. Changes
were uncommitted at that measurement stage.
