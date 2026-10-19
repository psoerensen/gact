param([string[]]$Stages=@('single_blr_bayesr'))
$ErrorActionPreference='Stop'
$taskRun=Join-Path (Resolve-Path 'work/revision-benchmark').Path 'fits-single-blr'
New-Item -ItemType Directory -Force -Path $taskRun | Out-Null
foreach($taskStage in $Stages){
 $taskOut=Join-Path $taskRun $taskStage
 if(Test-Path -LiteralPath $taskOut){throw "Output already exists: $taskOut"}
 New-Item -ItemType Directory -Path $taskOut | Out-Null
 $taskStart=[DateTime]::UtcNow
 $taskProc=Start-Process -FilePath 'C:/Program Files/R/R-4.4.1/bin/Rscript.exe' -ArgumentList @('--vanilla','work/single_blr_worker.R',$taskStage,$taskOut) -RedirectStandardOutput (Join-Path $taskOut 'stdout.log') -RedirectStandardError (Join-Path $taskOut 'stderr.log') -WindowStyle Hidden -PassThru
 $taskPeak=0L;$taskPrivate=0L;$taskCpu=0
 while(-not $taskProc.HasExited){
  $taskPidFile=Join-Path $taskOut 'pid.txt'
  if(Test-Path $taskPidFile){
   $taskRId=[int](Get-Content $taskPidFile)
   $taskR=Get-Process -Id $taskRId -ErrorAction SilentlyContinue
   if($null -ne $taskR){$taskPeak=[Math]::Max($taskPeak,$taskR.PeakWorkingSet64);$taskPrivate=[Math]::Max($taskPrivate,$taskR.PrivateMemorySize64);$taskCpu=$taskR.CPU}
  }
  Start-Sleep -Milliseconds 100
  $taskProc.Refresh()
 }
 $taskProc.WaitForExit()
 $taskCode=$taskProc.ExitCode
 [pscustomobject]@{stage=$taskStage;process_seconds=([DateTime]::UtcNow-$taskStart).TotalSeconds;peak_working_set_bytes=$taskPeak;sampled_private_bytes=$taskPrivate;cpu_seconds=$taskCpu;exit_code=$taskCode;memory_definition='R PID lifetime peak working set observed at 100ms polling'} | ConvertTo-Json | Set-Content (Join-Path $taskOut 'process.json')
 Write-Output "$taskStage finished: exit=$taskCode; peak=$taskPeak"
 if($taskCode -ne 0){throw "Stage failed: $taskStage"}
}
