$ErrorActionPreference='Stop'
if(-not(Test-Path 'work/revision-benchmark/input-preparation/timing.json')){throw 'Run input preparation first'}
$taskRoot=Join-Path (Resolve-Path 'work/revision-benchmark').Path 'scoring-grid'
New-Item -ItemType Directory -Force -Path $taskRoot | Out-Null
$taskCases=@()
foreach($taskN in @(500,2000)){foreach($taskM in @(10000,50000)){foreach($taskT in @(1,4)){
 $taskCases += [pscustomobject]@{n=$taskN;m=$taskM;k=0;overlap=0;threads=$taskT}
 foreach($taskK in @(3,20)){foreach($taskO in @(0,0.5)){$taskCases += [pscustomobject]@{n=$taskN;m=$taskM;k=$taskK;overlap=$taskO;threads=$taskT}}}
}}}
foreach($taskCase in $taskCases){
 $taskName='n{0}-m{1}-k{2}-o{3}-t{4}' -f $taskCase.n,$taskCase.m,$taskCase.k,([int]($taskCase.overlap*100)),$taskCase.threads
 $taskOut=Join-Path $taskRoot $taskName
 if(Test-Path $taskOut){throw "Existing result: $taskOut"}
 New-Item -ItemType Directory -Path $taskOut | Out-Null
 $taskStart=[DateTime]::UtcNow
 $taskP=Start-Process -FilePath 'C:/Program Files/R/R-4.4.1/bin/Rscript.exe' -ArgumentList @('--vanilla','work/scoring_grid_worker.R',$taskCase.n,$taskCase.m,$taskCase.k,$taskCase.overlap.ToString([Globalization.CultureInfo]::InvariantCulture),$taskCase.threads,$taskOut) -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $taskOut 'stdout.log') -RedirectStandardError (Join-Path $taskOut 'stderr.log')
 $taskPeak=0L;$taskPrivate=0L
 while(-not $taskP.HasExited){
  if(Test-Path (Join-Path $taskOut 'pid.txt')){
   $taskRId=[int](Get-Content (Join-Path $taskOut 'pid.txt'));$taskR=Get-Process -Id $taskRId -ErrorAction SilentlyContinue
   if($null -ne $taskR){$taskPeak=[Math]::Max($taskPeak,$taskR.PeakWorkingSet64);$taskPrivate=[Math]::Max($taskPrivate,$taskR.PrivateMemorySize64)}
  }
  Start-Sleep -Milliseconds 100;$taskP.Refresh()
 }
 $taskP.WaitForExit()
 [pscustomobject]@{peak_working_set_bytes=$taskPeak;sampled_private_bytes=$taskPrivate;process_seconds=([DateTime]::UtcNow-$taskStart).TotalSeconds;exit_code=$taskP.ExitCode} | ConvertTo-Json | Set-Content (Join-Path $taskOut 'process.json')
 Write-Output "$taskName exit=$($taskP.ExitCode)"
 if($taskP.ExitCode -ne 0){throw "Failed scoring case $taskName"}
}

