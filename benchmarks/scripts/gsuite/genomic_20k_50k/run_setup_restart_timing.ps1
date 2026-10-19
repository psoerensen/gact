param()
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path '.').Path
$runName = 'manuscript-setup-restart-' + (Get-Date -Format 'yyyyMMdd-HHmmss')
$runRoot = Join-Path $root ('build/benchmarks/work/' + $runName)
if (Test-Path -LiteralPath $runRoot) { throw 'Output destination exists' }
New-Item -ItemType Directory -Path $runRoot | Out-Null
$cpu = Get-CimInstance Win32_Processor
$system = Get-CimInstance Win32_ComputerSystem
$os = Get-CimInstance Win32_OperatingSystem
@{cpu=$cpu.Name; physicalCores=$cpu.NumberOfCores; logicalCores=$cpu.NumberOfLogicalProcessors; ramBytes=$system.TotalPhysicalMemory; os=$os.Caption; osVersion=$os.Version; storage=(Get-Volume -DriveLetter C | Select-Object DriveLetter,FileSystem,Size,SizeRemaining); memoryAPI='System.Diagnostics.Process.PeakWorkingSet64; R analysis process identified by Sys.getpid(), excluding launcher and helper processes'; processPrivateMetric='sampled maximum PrivateMemorySize64, 100ms interval'; timezone='Europe/Copenhagen'} | ConvertTo-Json -Depth 5 | Set-Content (Join-Path $runRoot 'hardware.json')
$stages = @('database_fresh','database_reuse')
$threadCounts = @(1)
$worker = 'benchmarks/genomic_20k_50k/setup_restart_stage.R'
$gactSource = Join-Path (Split-Path $root) 'gact'
@{baselineCommit='2510035'; candidateFiles=@(Get-ChildItem (Join-Path $gactSource 'R') -File | ForEach-Object { @{name=$_.Name; sha256=(Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash} }); scope='Single fresh setup followed by checksum-verified reuse; network and operating-system caches are uncontrolled'; library='build/gact-revision-cleanup/setup-r-library'} | ConvertTo-Json -Depth 6 | Set-Content (Join-Path $runRoot 'source-provenance.json')
$rows = @()
foreach ($threads in $threadCounts) {
 foreach ($stage in $stages) {
  $dest = Join-Path $runRoot ($stage + '-' + $threads)
  New-Item -ItemType Directory -Path $dest | Out-Null
  $env:OMP_NUM_THREADS = "$threads"
  $env:OPENBLAS_NUM_THREADS = "$threads"
  $env:MKL_NUM_THREADS = "$threads"
  $env:OMP_DYNAMIC = 'FALSE'
  $sw = [Diagnostics.Stopwatch]::StartNew()
  $proc = Start-Process -FilePath 'C:/Program Files/R/R-4.4.1/bin/x64/Rscript.exe' -ArgumentList @('--vanilla',$worker,$stage,"$threads",('"'+$dest+'"')) -WorkingDirectory $root -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $dest 'stdout.log') -RedirectStandardError (Join-Path $dest 'stderr.log')
  $peak = 0L; $private = 0L
  $handle = $proc.Handle
  $rProc = $null
  while (-not $proc.HasExited) {
   $proc.Refresh()
   if ($null -eq $rProc -and (Test-Path (Join-Path $dest "r-process.pid"))) { try { $rProc = Get-Process -Id ([int](Get-Content (Join-Path $dest "r-process.pid"))); $rHandle = $rProc.Handle } catch {} }
   if ($null -ne $rProc) { try { $rProc.Refresh(); $peak = [Math]::Max($peak,$rProc.PeakWorkingSet64); $private = [Math]::Max($private,$rProc.PrivateMemorySize64) } catch {} }
   try { $peak = [Math]::Max($peak,$proc.PeakWorkingSet64); $private = [Math]::Max($private,$proc.PrivateMemorySize64) } catch {}
   Start-Sleep -Milliseconds 100
  }
  $proc.WaitForExit(); $sw.Stop()
  try { $peak = [Math]::Max($peak,$proc.PeakWorkingSet64) } catch {}
  $bytes = (Get-ChildItem -LiteralPath $dest -Recurse -File | Measure-Object Length -Sum).Sum
  $row = [pscustomobject]@{stage=$stage; threads=$threads; processWallSeconds=$sw.Elapsed.TotalSeconds; peakWorkingSetBytes=$peak; sampledPeakPrivateBytes=$private; exitCode=$proc.ExitCode; stageDirectoryBytes=$bytes; path=$dest}
  $rows += $row
  $rows | Export-Csv (Join-Path $runRoot 'process-measurements.csv') -NoTypeInformation
  Write-Output ($row | ConvertTo-Json -Compress)
  if ($proc.ExitCode -ne 0) { throw ('Stage failed: ' + $stage + '; preserved at ' + $dest) }
 }
}
Write-Output ('RESULTS ' + $runRoot)






