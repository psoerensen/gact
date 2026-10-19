# CAD/T2D gcorr workflow

This bounded workload was requested for a gact tutorial. It uses the existing
CAD (GWAS2, PMID26343387) and T2D without UK Biobank (GWAS1, PMID30297969)
summary statistics, stored ordinary EUR LD scores and the existing 503-person
EUR BED reference. No GWAS/genotype download, database ingestion or persistent
LD/eigenvector files are part of the workload.

The worker adapts previously computed qgg genotype metadata into the public
gbase Glist descriptor, then validates it against current BED/BIM/FAM files.
It checks database summary-statistic allele orientation explicitly against the
reference. Follow the published gact reference preparation: BED is filtered to
`GAlist$rsids`, then the LD workflow uses all markers in that filtered Glist.
Select common reference/summary/LD-score rows. Apply no additional MAF,
missingness, HWE, palindrome, indel or duplicate-position filter. Retain
duplicate positions when marker IDs differ. The input audit counts retained
low-MAF, palindromic and position-duplicate markers without excluding them.

The workflow preserves extracted `stat$n`, including marker-specific values
when present. The audited cache is constant: GWAS1 49071.273294961 and GWAS2
40743.152404981. Both raw database summary files lack an `n` column, so extraction
uses study metadata `neff`, which for these binary traits is
`n_case*n_control/(n_case+n_control)`. Ingestion can also infer missing summary
sample sizes with `round(median(1/(2*af*(1-af)*seb^2)))`; extraction retains any
existing field. No factor-four, total-participant or liability conversion is
inferred. Results use the database information-size convention.

The 6 October runs incorrectly replaced this field with total participants
(456236 and 184305) and added a separate QC mask. Their numerical results are
superseded and preserved as historical evidence. The first 7 October diagnostic
changed only n, retaining that extra mask and full-reference normalization;
it is also superseded by the website-filtered run.

The current ordinary LDSC fit applies qgg's default squared-Z threshold
`max(.001*max(stat$n),80)`, which is 80 here. Native filtering is inclusive at
the boundary; qgg's is strict. No such Z cutoff is applied to chromosome-set
LDSC, as in qgg's set branch. Chromosome marker counts are matched set lengths.
The ordinary joint gcorr fit uses one common pre-filter marker-universe count,
whereas qgg uses separate post-filter counts per marginal/pair regression.
gcorr retains weighted regression and coordinated delete-block uncertainty;
qgg's regression is unweighted. Thus reference-filter alignment is not a claim
that the estimators are identical. Full BED counts are separate provenance.
The common-row mask is inherited from the existing three-study GWAS1/GWAS2/
GWAS6 LDSC cache. Only GWAS1 and GWAS2 are fitted here; GWAS6 contributes to
marker availability. The downloadable helper reproduces this mask explicitly.

The regional partition is the published Berisa/Pickrell EUR table at
`ac125e47bf7ff3e90be31f278a7b6a61daaba0dc`, retrieved from
<https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/EUR/fourier_ls-all.bed>.
The accompanying BED is a maintained small reference input, not generated LD.
The reference BIM rs429358 coordinate 19:45411941 agrees with GRCh37, as
documented by [NCBI](https://www.ncbi.nlm.nih.gov/variation/view/help/).
Regions use start-inclusive/stop-exclusive coordinates. Uncovered markers are
counted explicitly and excluded only from the joint regional modeled universe.

LDSC estimates marginal and cross intercepts. Ordinary GNOVA conditions on
the same-run fitted LDSC cross intercept; its uncertainty does not include
estimation of that nuisance input. SumHer illustrates an equal-marker baseline
using legacy ordinary scores for architecture and tagging; it is not a
MAF/LD-weighted LDAK architecture or an upstream LDAK software comparison.
Legacy score window/correction provenance is not reconstructed by recomputation.
Treat these genome-wide estimates as workflow illustrations conditional on the
stored score definition, rather than a freshly qualified manuscript result.

The corrected chr22 HESS/rho-HESS pilot uses K <= 50; the corrected genome
illustration explicitly uses K <= 20. Both use eigenvalues > 1 and relative
tolerance 1e-10. The
503-by-marker standardized regional genotype matrix defines complete signed
LD; the smaller Gram is decomposed. All traits share that decomposition.
Compact quadratic/cross moments are passed to the existing joint correction
and uncertainty formulas. Chr22 is a numerical/resource pilot. Full-genome
regional results assume negligible inter-region LD and the supplied regions
as the modeled genetic universe. Negative estimates and missing uncertainty
are retained, with their reasons; no PSD repair is applied.

The historical full-genome projection retained 83728 components in total.
That exceeds both corrected information sizes. The joint correction requires
`1-total_rank/n > 0`, so this full-genome fit is inadmissible with these inputs.
Do not bypass the guard or interpret the old total-N estimates as corrected
results. The corrected genome illustration caps every region at 20 components,
selected before inspecting trait projections. With 1702 regions, its conservative
rank bound is 34040, giving a correction-denominator bound above 0.16 for CAD.
The worker checks this bound before preparation and saves actual retained ranks
and denominators afterward. This choice changes the spectral model; it is not
an automatic rank selection rule or a change to package defaults. The
[HESS FAQ](https://huwenboshi.github.io/hess/faq/) recommends sample sizes above
50000 and notes possible downward bias when smaller samples require stronger
truncation. These extracted information sizes are below that threshold. Neither
the lower cap nor a successful numerical solve establishes statistical calibration.

The completed corrected genome illustration retains 33971 components and passes
rank admission, but returns a negative CAD total (-0.187175), undefined total
rg and unavailable conditional uncertainty. The chr22 pilot likewise returns
a negative CAD contribution and undefined total rg. These raw results are
preserved with native reasons; they do not provide interpretable CAD/T2D
genetic-correlation inference under the current HESS specification. Completed
GNOVA and SumHer estimates and all timing/memory records accompany the
corrected tutorial downloads and the local remaining-analyses evidence record.

Overlap counts and phenotypic correlation are not recorded in the database.
Any rho-HESS examples must state their overlap scenario explicitly; an LDSC
cross intercept is not a phenotypic correlation. Empirical claims requiring
these values remain conditional. Preparation is reusable for sensitivity fits.
The user subsequently stated their understanding that CAD and T2D samples do
not overlap. The chromosome-set analysis therefore uses zero shared samples
as its working assumption and fixes the cross-trait LDSC intercept to zero;
this is not independent verification of cohort membership. A free cross
intercept can also reflect shared confounding, so this constraint is a model
choice beyond the sample-overlap assertion. Marginal intercepts remain free.

Chromosomes are represented as an ordinary named `sets` list of marker IDs.
`gcorr_set_scores.R` constructs disjoint set columns by masking the existing
ordinary LD scores, with matched set counts and diagonal binary
annotation overlap. The existing `gcorr_ldsc_partitioned()` jointly fits those
columns and reports set-specific h2, covariance and rg with coordinated block
jackknife uncertainty. This construction assumes negligible LD across set
boundaries; chromosome groups use the within-chromosome LD approximation.
Arbitrary interleaved gene/feature sets instead require their true annotation
LD scores. No annotation matrix is saved to disk. No new eigendecomposition,
genome-wide LD file, native solver or qgg comparison is part of this extension.

Run through the existing process-measured launcher:

```powershell
./benchmarks/genomic_20k_50k/run_manuscript.ps1 -GcorrComparison `
  -OnlyStages gcorr_inputs -AnalysisThreads 4 -TimeoutSeconds 3600
# Set FixtureRds to that run's gcorr_inputs-4/output.rds.
# For a historical cache lacking n, also pass -SummaryNRds with the audited
# aligned-database-n.rds. The worker rejects uncorrected legacy caches.
./benchmarks/genomic_20k_50k/run_manuscript.ps1 -GcorrComparison `
  -FixtureRds '<input RDS>' -OnlyStages gcorr_chr22 -AnalysisThreads 4
./benchmarks/genomic_20k_50k/run_manuscript.ps1 -GcorrComparison `
  -FixtureRds '<input RDS>' `
  -LdscEstimatesCsv '<corrected LDSC estimates.csv>' `
  -OnlyStages 'gcorr_gnova,gcorr_sumher' `
  -AnalysisThreads 4 -TimeoutSeconds 1200
./benchmarks/genomic_20k_50k/run_manuscript.ps1 -GcorrComparison `
  -FixtureRds '<input RDS>' -OnlyStages gcorr_genome `
  -HessMaximumComponents 20 -AnalysisThreads 4 -TimeoutSeconds 3600
./benchmarks/genomic_20k_50k/run_manuscript.ps1 -GcorrComparison `
  -FixtureRds '<input RDS>' -OnlyStages gcorr_sets `
  -AnalysisThreads 4 -TimeoutSeconds 900
```

The full-genome HESS stage requires the explicit smaller cap above. Historical
K <= 50 fails the corrected information sizes' total-rank condition. A same-input
corrected LDSC CSV supplies GNOVA's nuisance intercept without repeating LDSC;
the source path, file checksum and used intercept are recorded with the run.
The earlier overlap sensitivity counts used physical total N; the corrected
regional worker retains only the user's zero-overlap scenario.

Authoritative website sources:

- [Filter reference BED to database marker IDs](https://psoerensen.github.io/gact/Document/Process_1000G.html)
- [Compute LD for every marker in filtered Glist](https://psoerensen.github.io/gact/Document/Compute_sparseLD_1000G.html)
- [Use stat$n and named chromosome sets in LDSC](https://psoerensen.github.io/gact/Document/LD_score_regression.html)

The launcher records hardware and actual R-process peak working set/private
memory. Stage clocks exclude input RDS loading; `fit-phase.csv` includes
ordinary R adapter work. Regional preparation includes validation and reader
setup, with decode/product/eigen/projection clocks saved per region. Serialization
is separate. Record results only after successful completion with native gcorr
0.1.1 and the compatible isolated R installation. No broad calibration or
ecosystem qualification is implied.

Summarize the completed pilot and genome runs without rerunning analyses:

```sh
Rscript --vanilla benchmarks/genomic_20k_50k/summarize_gcorr_cad_t2d.R \
  '<pilot run>' '<score run>' '<regional run>' '<validated input-stage directory>' '<output directory>'
```

This writes compact CSV tables, two plots and tutorial text fragments. It does
not export GWAS marker-level statistics, genotypes, LD or eigenvectors.
