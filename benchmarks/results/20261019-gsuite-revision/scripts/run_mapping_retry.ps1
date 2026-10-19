$ErrorActionPreference='Stop'
$taskRoot=Join-Path (Resolve-Path 'work/revision-benchmark').Path 'remaining-modules'
New-Item -ItemType Directory -Force -Path $taskRoot | Out-Null
foreach($taskStage in @('regional_twochain','genome_twochain')){
 $taskOut=Join-Path $taskRoot $taskStage
 if(Test-Path $taskOut){throw "Existing result $taskOut"}
 New-Item -ItemType Directory -Path $taskOut | Out-Null
 $taskStart=[DateTime]::UtcNow
 $taskP=Start-Process -FilePath 'C:/Program Files/R/R-4.4.1/bin/Rscript.exe' -ArgumentList @('--vanilla','work/component_mapping_worker_retry.R',$taskStage,$taskOut) -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $taskOut 'stdout.log') -RedirectStandardError (Join-Path $taskOut 'stderr.log')
 $taskPeak=0L;$taskPrivate=0L
 while(-not $taskP.HasExited){
  if(Test-Path (Join-Path $taskOut 'pid.txt')){
   $taskRId=[int](Get-Content (Join-Path $taskOut 'pid.txt'));$taskR=Get-Process -Id $taskRId -ErrorAction SilentlyContinue
   if($null -ne $taskR){$taskPeak=[Math]::Max($taskPeak,$taskR.PeakWorkingSet64);$taskPrivate=[Math]::Max($taskPrivate,$taskR.PrivateMemorySize64)}
  }
  Start-Sleep -Milliseconds 100;$taskP.Refresh()
 }
 $taskP.WaitForExit()
 if(-not (Test-Path (Join-Path $taskOut 'timing.json'))){
  [pscustomobject]@{stage=$taskStage;status='failed';operation_seconds=$null;error=((Get-Content (Join-Path $taskOut 'stderr.log') -Tail 8) -join ' ');output_bytes=0} | ConvertTo-Json | Set-Content (Join-Path $taskOut 'timing.json')
 }
 [pscustomobject]@{peak_working_set_bytes=$taskPeak;sampled_private_bytes=$taskPrivate;process_seconds=([DateTime]::UtcNow-$taskStart).TotalSeconds;exit_code=$taskP.ExitCode} | ConvertTo-Json | Set-Content (Join-Path $taskOut 'process.json')
 Write-Output "$taskStage exit=$($taskP.ExitCode)"
 if($taskP.ExitCode -ne 0){Write-Output "Module $taskStage failed; preserving its diagnostics and continuing independent modules."}
}
