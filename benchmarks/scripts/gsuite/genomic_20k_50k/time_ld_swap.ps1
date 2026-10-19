param([ValidateRange(1,10)][int]$Repeats=3)
$ErrorActionPreference = 'Stop'
$env:LC_ALL = 'C'
$env:LANG = 'C'
$root = (Resolve-Path '.').Path
$evidence = Join-Path $root ('build/benchmarks/work/ld-swap-' + (Get-Date -Format 'yyyyMMdd-HHmmss'))
if (Test-Path -LiteralPath $evidence) { throw 'Output destination exists' }
New-Item -ItemType Directory -Path $evidence | Out-Null
$sourceRows = @()
foreach ($repo in @('gsuite','gbayes','gmap','gmat','gsolve')) {
    $checkout = (Resolve-Path -LiteralPath (Join-Path $root ('../' + $repo))).Path
    $safe = 'safe.directory=' + $checkout.Replace('\','/')
    $head = & git -c $safe -C $checkout rev-parse HEAD
    if ($LASTEXITCODE -ne 0) { throw ('Cannot identify source: ' + $repo) }
    $diffPath = Join-Path $evidence ($repo + '-source.diff')
    & git -c $safe -C $checkout diff HEAD "--output=$diffPath"
    if ($LASTEXITCODE -ne 0) { throw ('Cannot preserve source diff: ' + $repo) }
    $sourceRows += [pscustomobject]@{repository=$repo; head=$head}
}
$sourceRows | Export-Csv (Join-Path $evidence 'source-commits.csv') -NoTypeInformation
$installed = foreach ($package in @('gbase','gmap','gbayes')) {
    Get-ChildItem -LiteralPath (Join-Path $root ('build/task-packages/library/' + $package)) -Recurse -File |
        Get-FileHash -Algorithm SHA256 | Select-Object Path,Algorithm,Hash
}
$installed | Export-Csv (Join-Path $evidence 'installed-package-sha256.csv') -NoTypeInformation
foreach ($script in @('run_manuscript.ps1','simulated_finemapping_comparison_stage.R','time_ld_swap.ps1')) {
    Copy-Item -LiteralPath (Join-Path $root ('benchmarks/genomic_20k_50k/' + $script)) -Destination $evidence
}
$runs = @()
foreach ($repeat in 1..$Repeats) {
    # Alternate serial order; do not benchmark competing R processes.
    $probabilities = if ($repeat % 2 -eq 1) { @(0.0,0.05) } else { @(0.05,0.0) }
    foreach ($probability in $probabilities) {
        $output = @(& (Join-Path $root 'benchmarks/genomic_20k_50k/run_manuscript.ps1') `
            -SimComparison -QggInterfaces -LdSwapProbability $probability `
            -QggCompareRun build/benchmarks/manuscript-20261005-231349 `
            -OnlyStages 'cmp_suite_regional,cmp_suite_genome')
        $output | Write-Output
        $result = @($output | Where-Object { $_ -like 'RESULTS *' })
        if ($result.Count -ne 1) { throw 'Cannot identify completed run' }
        $runs += [pscustomobject]@{repeat=$repeat; probability=$probability; run=$result[0].Substring(8)}
        $runs | Export-Csv (Join-Path $evidence 'runs.csv') -NoTypeInformation
    }
}
Write-Output ('EVIDENCE ' + $evidence)
