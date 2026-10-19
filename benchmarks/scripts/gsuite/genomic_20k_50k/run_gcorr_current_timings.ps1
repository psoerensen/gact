param([ValidateSet('ldsc','gnova','sumher','sets','hess_prepare','hess_raw','hess_partial')]
  [string[]]$Stages=@('ldsc','gnova','sumher','sets','hess_prepare','hess_raw','hess_partial'))
$ErrorActionPreference='Stop'
[System.Threading.Thread]::CurrentThread.CurrentCulture=[System.Globalization.CultureInfo]::InvariantCulture
$taskRoot=(Resolve-Path '.').Path
$taskRun=Join-Path $taskRoot ('build/benchmarks/work/gcorr-current-'+(Get-Date -Format 'yyyyMMdd-HHmmss'))
New-Item -ItemType Directory -Path $taskRun | Out-Null
$env:PATH=(Resolve-Path 'build/task-packages/prefix/bin').Path+';C:/rtools44/x86_64-w64-mingw32.static.posix/bin;'+$env:PATH
$env:LC_ALL='C'; $env:LANG='C'
$rscript='C:/Program Files/R/R-4.4.1/bin/x64/Rscript.exe'
Write-Host "GCORR_RUN $taskRun"
$rows=@()
foreach($stage in $Stages) {
  $out=Join-Path $taskRun $stage
  New-Item -ItemType Directory -Path $out | Out-Null
  $worker=Join-Path $out 'worker.R'
  Copy-Item -LiteralPath (Join-Path $taskRoot 'benchmarks/genomic_20k_50k/gcorr_memory_stage.R') -Destination $worker
  $child=Start-Process -FilePath $rscript -ArgumentList @('--vanilla',('"'+$worker+'"'),$stage,('"'+$out+'"')) -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $out 'stdout.log') -RedirectStandardError (Join-Path $out 'stderr.log')
  $null=$child.Handle
  $baseline=$null; $peakResident=0L; $peakPrivate=0L; $samples=0L
  $memoryRows=[System.Collections.Generic.List[object]]::new()
  $phaseClock=[System.Diagnostics.Stopwatch]::new()
  while(-not $child.HasExited) {
    $child.Refresh()
    if((Test-Path -LiteralPath (Join-Path $out 'ready')) -and $null -eq $baseline) {
      $baseline=@{resident=$child.WorkingSet64; private=$child.PrivateMemorySize64}
      $baseline | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $out 'baseline-memory.json')
      Set-Content -LiteralPath (Join-Path $out 'go') -Value 'go'
      $phaseClock.Start()
      Write-Host "$stage computation started"
    }
    if($null -ne $baseline -and -not (Test-Path -LiteralPath (Join-Path $out 'done'))) {
      $peakResident=[Math]::Max($peakResident,$child.WorkingSet64)
      $peakPrivate=[Math]::Max($peakPrivate,$child.PrivateMemorySize64)
      $samples++
      $memoryRows.Add([pscustomobject]@{seconds=$phaseClock.Elapsed.TotalSeconds;resident_bytes=$child.WorkingSet64;private_bytes=$child.PrivateMemorySize64})
    }
    Start-Sleep -Milliseconds 200
  }
  $child.WaitForExit()
  if(($null -ne $child.ExitCode -and $child.ExitCode -ne 0) -or -not (Test-Path -LiteralPath (Join-Path $out 'success'))) {
    Get-Content -LiteralPath (Join-Path $out 'stderr.log') -Tail 14
    throw "$stage failed with exit code $($child.ExitCode)"
  }
  $timing=Import-Csv -LiteralPath (Join-Path $out 'timing.csv')
  $memoryRows | Export-Csv -LiteralPath (Join-Path $out 'memory-samples.csv') -NoTypeInformation
  $rows+=[pscustomobject]@{stage=$stage;seconds=[double]$timing.seconds;baseline_resident_gib=$baseline.resident/1GB;peak_resident_gib=$peakResident/1GB;resident_growth_gib=[Math]::Max([double]0,[double](($peakResident-$baseline.resident)/1GB));baseline_private_gib=$baseline.private/1GB;peak_private_gib=$peakPrivate/1GB;samples=$samples}
  $rows | Export-Csv -LiteralPath (Join-Path $taskRun 'timings-memory.csv') -NoTypeInformation
  Write-Host "$stage completed: $($timing.seconds) seconds; peak $([Math]::Round($peakResident/1GB,3)) GiB"
}
Write-Host "GCORR_TIMINGS_COMPLETE $taskRun"
