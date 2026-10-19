$ErrorActionPreference='Stop'
$taskOut=Join-Path (Resolve-Path 'work/summary-blr-benchmark').Path 'native-fit'
if(Test-Path $taskOut){throw 'Native result directory already exists'}
New-Item -ItemType Directory -Path $taskOut | Out-Null
$taskPrefix='C:/Users/au223366/Documents/GitHub/gsuite/build/task-packages/prefix'
$env:PATH="$taskPrefix/bin;C:/rtools44/x86_64-w64-mingw32.static.posix/bin;"+$env:PATH
$env:OMP_NUM_THREADS='1';$env:OPENBLAS_NUM_THREADS='1';$env:MKL_NUM_THREADS='1'
$taskStart=[DateTime]::UtcNow
$taskP=Start-Process -FilePath (Resolve-Path 'work/summary-blr-benchmark/fit_precomputed.exe').Path -ArgumentList @((Resolve-Path 'work/summary-blr-benchmark/input').Path,$taskOut) -WindowStyle Hidden -PassThru -RedirectStandardOutput (Join-Path $taskOut 'stdout.log') -RedirectStandardError (Join-Path $taskOut 'stderr.log')
$taskPeak=0L;$taskPrivate=0L;$taskCpu=0
while(-not $taskP.HasExited){
 $taskP.Refresh();$taskPeak=[Math]::Max($taskPeak,$taskP.PeakWorkingSet64);$taskPrivate=[Math]::Max($taskPrivate,$taskP.PrivateMemorySize64);$taskCpu=$taskP.CPU
 Start-Sleep -Milliseconds 100
}
$taskP.WaitForExit()
[pscustomobject]@{exit_code=$taskP.ExitCode;process_seconds=([DateTime]::UtcNow-$taskStart).TotalSeconds;cpu_seconds=$taskCpu;peak_working_set_bytes=$taskPeak;sampled_private_bytes=$taskPrivate;memory_scope='Native executable only; R preparation is a separate process'} | ConvertTo-Json | Set-Content (Join-Path $taskOut 'process.json')
Get-Content (Join-Path $taskOut 'stdout.log')
if($taskP.ExitCode -ne 0){Get-Content (Join-Path $taskOut 'stderr.log');throw 'Native sampler failed'}
