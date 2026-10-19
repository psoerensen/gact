param([switch]$GcorrComparison, [string]$SummaryNRds, [string]$LdscEstimatesCsv, [ValidateRange(1,2048)][int]$HessMaximumComponents=50, [switch]$GlmaComparison, [switch]$Database, [switch]$RealLD, [switch]$Qgg, [string]$ReferenceRds, [string]$CohortRds, [string]$FixtureRds, [string]$StartStage, [string]$QggReuseRun, [string]$OnlyStages, [ValidateRange(1,12)][int]$AnalysisThreads=1, [switch]$PackageComparison, [string]$QggCompareRun, [switch]$GenePreparation, [switch]$WebsiteExamples, [switch]$RemainingExamples, [switch]$SimComparison, [switch]$QggInterfaces, [ValidateRange(0,1)][double]$LdSwapProbability=0, [ValidateRange(0,3600)][int]$TimeoutSeconds=0)
$ErrorActionPreference = 'Stop'
$glmaBenchmark = [bool]$GlmaComparison
if ($GlmaComparison) { $env:GLMA_SIM_FIXTURE = if ($FixtureRds) { (Resolve-Path -LiteralPath $FixtureRds).Path } else { '' } }
if ($QggInterfaces -and -not $SimComparison) { throw 'QggInterfaces requires SimComparison' }
if ($LdSwapProbability -gt 0 -and -not ($SimComparison -and $QggInterfaces)) { throw 'LdSwapProbability requires SimComparison and QggInterfaces' }
$env:SIM_COMPARE_QGG_INTERFACES = if ($QggInterfaces) { '1' } else { '0' }
$env:SIM_COMPARE_LD_SWAP_PROBABILITY = $LdSwapProbability.ToString([Globalization.CultureInfo]::InvariantCulture)
if ($QggReuseRun) { $env:QGG_BENCH_REUSE_RUN = (Resolve-Path -LiteralPath $QggReuseRun).Path }
if ($QggCompareRun) { $env:QGG_COMPARE_RUN = (Resolve-Path -LiteralPath $QggCompareRun).Path }
if ($SimComparison -and $QggCompareRun) { $env:SIM_COMPARE_REUSE = (Resolve-Path -LiteralPath $QggCompareRun).Path }
if ($GcorrComparison) { $env:LC_ALL='C'; $env:LANG='C' }
if ($GcorrComparison -and $FixtureRds) { $env:GCORR_CAD_INPUT = (Resolve-Path -LiteralPath $FixtureRds).Path }
if ($SummaryNRds -and -not $GcorrComparison) { throw 'SummaryNRds requires GcorrComparison' }
$env:GCORR_SUMMARY_N_RDS = if ($SummaryNRds) { (Resolve-Path -LiteralPath $SummaryNRds).Path } else { '' }
if ($LdscEstimatesCsv -and -not $GcorrComparison) { throw 'LdscEstimatesCsv requires GcorrComparison' }
$env:GCORR_LDSC_ESTIMATES_CSV = if ($LdscEstimatesCsv) { (Resolve-Path -LiteralPath $LdscEstimatesCsv).Path } else { '' }
$env:GCORR_HESS_MAXIMUM_COMPONENTS = $HessMaximumComponents.ToString([Globalization.CultureInfo]::InvariantCulture)
$root = (Resolve-Path '.').Path
$runName = 'manuscript-' + (Get-Date -Format 'yyyyMMdd-HHmmss')
$runRoot = Join-Path $root ('build/benchmarks/work/' + $runName)
if (Test-Path -LiteralPath $runRoot) { throw 'Output destination exists' }
New-Item -ItemType Directory -Path $runRoot | Out-Null
$cpu = Get-CimInstance Win32_Processor
$system = Get-CimInstance Win32_ComputerSystem
$os = Get-CimInstance Win32_OperatingSystem
@{cpu=$cpu.Name; physicalCores=$cpu.NumberOfCores; logicalCores=$cpu.NumberOfLogicalProcessors; ramBytes=$system.TotalPhysicalMemory; os=$os.Caption; osVersion=$os.Version; storage=(Get-Volume -DriveLetter C | Select-Object DriveLetter,FileSystem,Size,SizeRemaining); memoryAPI='System.Diagnostics.Process.PeakWorkingSet64; R analysis process identified by Sys.getpid(), excluding launcher and helper processes'; processPrivateMetric='sampled maximum PrivateMemorySize64, 100ms interval'; timezone='Europe/Copenhagen'} | ConvertTo-Json -Depth 5 | Set-Content (Join-Path $runRoot 'hardware.json')
$stages = if ($Qgg) { @('qgg_download','qgg_glist','qgg_ld','qgg_inputs','qgg_vegas','qgg_joint','qgg_mapping','qgg_bayesc_fit','qgg_total_c','qgg_bayesr_fit','qgg_total_r','qgg_pathway') } elseif ($Database) { @('database') } elseif ($RealLD) { @('reference_import','founders','bed_conversion','preparation','vegas_real','joint_real','mapping_real','weights_real','total_real','pathway_real') } else { @('vegas','joint','mapping','total','pathway') }
if ($RealLD -and $ReferenceRds) { $env:GSIM_BENCH_REFERENCE_RDS = (Resolve-Path -LiteralPath $ReferenceRds).Path; $stages = @($stages | Where-Object { $_ -ne 'reference_import' }) }
if ($RealLD -and $CohortRds) { $env:GSIM_BENCH_COHORT_RDS = (Resolve-Path -LiteralPath $CohortRds).Path; $stages = @($stages | Where-Object { $_ -notin @('reference_import','founders','bed_conversion') }) }
if ($RealLD -and $FixtureRds) { $env:GSIM_BENCH_FIXTURE_RDS = (Resolve-Path -LiteralPath $FixtureRds).Path; $stages = @($stages | Where-Object { $_ -notin @('reference_import','founders','bed_conversion','preparation') }) }
if ($StartStage) { $index = [Array]::IndexOf($stages,$StartStage); if ($index -lt 0) { throw 'StartStage is absent from selected stages' }; $stages = $stages[$index..($stages.Length-1)] }
$threadCounts = if ($Qgg) { @($AnalysisThreads) } elseif ($Database -or $RealLD) { @(1) } else { @(1,4) }
if ($PackageComparison) { $stages=@('suite_joint','suite_mapping','suite_vegas'); $threadCounts=@($AnalysisThreads) }
if ($GenePreparation) { $stages=@('gene_blocks'); $threadCounts=@($AnalysisThreads) }
if ($WebsiteExamples) { $stages=@('website_inputs','website_pops_rr','website_pops_c','website_pops_r','website_magma_b_joint','website_magma_d_joint','website_magma_v_joint','website_magma_b_bayesr','website_magma_d_bayesr','website_magma_v_bayesr'); $threadCounts=@(1) }
if ($RemainingExamples) { $stages=@('website_ldsc_inputs','website_ldsc_h2','website_ldsc_rg','website_ldsc_chr','website_ldsc_reg','website_sim_glist','website_sim_filter','website_sim_ld','website_sim_stats','website_sim_regional','website_sim_genome','website_sim_credible'); $threadCounts=@(1) }
if ($SimComparison) { $stages=@('cmp_qgg_ld','cmp_qgg_regional','cmp_qgg_genome','cmp_qgg_credible','cmp_suite_prepare','cmp_suite_ld','cmp_suite_regional','cmp_suite_genome','cmp_suite_credible_oracle','cmp_compare'); $threadCounts=@(1) }
if ($glmaBenchmark) { $stages=if ($FixtureRds) { @('glma_linear','glma_loco','glma_prior') } else { @('glma_fixture','glma_linear','glma_loco','glma_prior') }; $threadCounts=@($AnalysisThreads) }
if ($GcorrComparison) { $stages=@('gcorr_inputs','gcorr_ldsc','gcorr_gnova','gcorr_sumher','gcorr_chr22','gcorr_genome','gcorr_sets'); $threadCounts=@($AnalysisThreads) }
if ($OnlyStages) { $selected = $OnlyStages.Split(','); if (@($selected | Where-Object { $_ -notin $stages }).Count) { throw 'OnlyStages contains an unknown stage' }; $stages = @($stages | Where-Object { $_ -in $selected }) }
$worker = if ($Qgg) { 'benchmarks/genomic_20k_50k/qgg_examples_stage.R' } elseif ($RealLD) { 'benchmarks/genomic_20k_50k/real_ld_stage.R' } else { 'benchmarks/genomic_20k_50k/manuscript_stage.R' }
if ($PackageComparison) { $worker='benchmarks/genomic_20k_50k/compare_qgg_packages_stage.R' }
if ($GenePreparation) { $worker='benchmarks/genomic_20k_50k/gene_blocks_stage.R' }
if ($WebsiteExamples) { $worker='benchmarks/genomic_20k_50k/website_analysis_stage.R' }
if ($RemainingExamples) { $worker='benchmarks/genomic_20k_50k/website_remaining_stage.R' }
if ($SimComparison) { $worker='benchmarks/genomic_20k_50k/simulated_finemapping_comparison_stage.R' }
if ($glmaBenchmark) { $worker='benchmarks/genomic_20k_50k/glma_simulated_stage.R' }
if ($GcorrComparison) { $worker='benchmarks/genomic_20k_50k/gcorr_cad_t2d_stage.R' }
$rows = @()
$keepAwake = $false
if ($GlmaComparison) {
 if (-not ('GlmaBenchmarkPower' -as [type])) {
 Add-Type -TypeDefinition @'
using System.Runtime.InteropServices;
public static class GlmaBenchmarkPower {
 [DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint flags);
}
'@
 }
 $previousExecutionState = [GlmaBenchmarkPower]::SetThreadExecutionState([uint32]2147483649)
 if ($previousExecutionState -eq 0) { throw 'Cannot keep the benchmark system awake' }
 $keepAwake = $true
}
try {
foreach ($threads in $threadCounts) {
 foreach ($stage in $stages) {
  $dest = Join-Path $runRoot ($stage + '-' + $threads)
  New-Item -ItemType Directory -Path $dest | Out-Null
  $workerSnapshot = Join-Path $dest 'worker.R'
  Copy-Item -LiteralPath (Join-Path $root $worker) -Destination $workerSnapshot
  $env:OMP_NUM_THREADS = "$threads"
  $env:OPENBLAS_NUM_THREADS = "$threads"
  $env:MKL_NUM_THREADS = "$threads"
  $env:OMP_DYNAMIC = 'FALSE'
  $sw = [Diagnostics.Stopwatch]::StartNew()
  $proc = Start-Process -FilePath 'C:/Program Files/R/R-4.4.1/bin/x64/Rscript.exe' -ArgumentList @('--vanilla',('"'+$workerSnapshot+'"'),$stage,"$threads",('"'+$dest+'"')) -WorkingDirectory $root -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $dest 'stdout.log') -RedirectStandardError (Join-Path $dest 'stderr.log')
  $peak = 0L; $private = 0L; $rCpuSeconds = 0.0
  $handle = $proc.Handle
  $rProc = $null; $timedOut = $false; $nativeModulesCaptured = $false
  while (-not $proc.HasExited) {
   $proc.Refresh()
   if ($null -eq $rProc -and (Test-Path (Join-Path $dest "r-process.pid"))) { try { $rProc = Get-Process -Id ([int](Get-Content (Join-Path $dest "r-process.pid"))); $rHandle = $rProc.Handle } catch {} }
   if ($null -ne $rProc) { try { $rProc.Refresh(); $peak = [Math]::Max($peak,$rProc.PeakWorkingSet64); $private = [Math]::Max($private,$rProc.PrivateMemorySize64); $rCpuSeconds = $rProc.TotalProcessorTime.TotalSeconds } catch {} }
   if (($glmaBenchmark -or $GcorrComparison) -and $null -ne $rProc -and -not $nativeModulesCaptured) {
    try {
     $nativeModules = @($rProc.Modules | Where-Object {$_.ModuleName -in @('gbase.dll','gcorr.dll','libgcorr.dll','glma.dll','libglma.dll','libgsolve.dll','libgreml.dll','libgbayes.dll','libgmat.dll','gbits.dll','gbits_gsolve_adapter.dll')} | Select-Object ModuleName,FileName)
     if ($nativeModules.ModuleName -contains $(if ($GcorrComparison) {'libgcorr.dll'}else{'libglma.dll'})) {
      $nativeModules | Export-Csv -LiteralPath (Join-Path $dest 'loaded-native-modules.csv') -NoTypeInformation
      $nativeModulesCaptured = $true
     }
    } catch {}
   }
   try { $peak = [Math]::Max($peak,$proc.PeakWorkingSet64); $private = [Math]::Max($private,$proc.PrivateMemorySize64) } catch {}
   if ($TimeoutSeconds -gt 0 -and $sw.Elapsed.TotalSeconds -ge $TimeoutSeconds) {
    $timedOut = $true
    if ($null -ne $rProc -and -not $rProc.HasExited) { $rProc.Kill() }
    if (-not $proc.HasExited) { $proc.Kill() }
    break
   }
   Start-Sleep -Milliseconds 100
  }
  $proc.WaitForExit(); $sw.Stop()
  if ($null -ne $rProc) { try { $rCpuSeconds = $rProc.TotalProcessorTime.TotalSeconds } catch {} }
  try { $peak = [Math]::Max($peak,$proc.PeakWorkingSet64) } catch {}
  $bytes = (Get-ChildItem -LiteralPath $dest -Recurse -File | Measure-Object Length -Sum).Sum
  $row = [pscustomobject]@{stage=$stage; threads=$threads; processWallSeconds=$sw.Elapsed.TotalSeconds; rProcessCpuSeconds=$rCpuSeconds; peakWorkingSetBytes=$peak; sampledPeakPrivateBytes=$private; exitCode=$proc.ExitCode; timedOut=$timedOut; timeoutSeconds=$TimeoutSeconds; stageDirectoryBytes=$bytes; path=$dest}
  $rows += $row
  $rows | Export-Csv (Join-Path $runRoot 'process-measurements.csv') -NoTypeInformation
  Write-Output ($row | ConvertTo-Json -Compress)
  if ($proc.ExitCode -ne 0) { throw ('Stage failed: ' + $stage + '; preserved at ' + $dest) }
 }
}
} finally {
 if ($keepAwake) { [void][GlmaBenchmarkPower]::SetThreadExecutionState($previousExecutionState) }
}
Write-Output ('RESULTS ' + $runRoot)
