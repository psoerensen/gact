$ErrorActionPreference='Stop'
$taskOut=Join-Path (Resolve-Path 'work/revision-benchmark').Path 'input-preparation'
if(Test-Path $taskOut){throw 'Input preparation output already exists'}
New-Item -ItemType Directory -Path $taskOut | Out-Null
$taskStart=[DateTime]::UtcNow
$taskP=Start-Process -FilePath 'C:/Program Files/R/R-4.4.1/bin/Rscript.exe' -ArgumentList @('--vanilla','work/input_preparation_worker.R',$taskOut) -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $taskOut 'stdout.log') -RedirectStandardError (Join-Path $taskOut 'stderr.log')
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
if($taskP.ExitCode -ne 0){throw 'Input preparation failed'}
Write-Output 'Input preparation completed and matched the model fixture.'
