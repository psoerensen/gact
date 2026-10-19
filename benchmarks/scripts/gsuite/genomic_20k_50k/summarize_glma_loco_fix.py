"""Summarize completed and capped standalone LOCO runs without conflating them."""
import argparse
import csv
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('runs', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
rows = []
def number(value):
    # PowerShell Export-Csv follows the host's Danish decimal separator;
    # R's timing CSV uses the explicit C locale.
    return float(str(value).replace(',', '.'))
for run in json.loads(args.runs.read_text(encoding='utf-8-sig')):
    root = Path(run['path'])
    with (root / 'process-measurements.csv').open(encoding='utf-8-sig') as stream:
        process = next(csv.DictReader(stream))
    stage = Path(process['path'])
    timing = {}
    if (stage / 'timing.csv').exists():
        with (stage / 'timing.csv').open(encoding='utf-8-sig') as stream:
            timing = next(csv.DictReader(stream))
    timeout = process.get('timedOut', '').lower() == 'true'
    status = 'timeout' if timeout else ('completed' if process['exitCode'] == '0'
              and timing.get('status') == 'completed' else 'failed')
    rows.append(dict(label=run['label'], status=status, threads=run['threads'],
        analysis_seconds=timing.get('analysis_seconds', ''),
        process_seconds=number(process['processWallSeconds']),
        cpu_seconds=number(process['rProcessCpuSeconds']),
        peak_resident_mib=number(process['peakWorkingSetBytes']) / 2**20,
        sampled_private_mib=number(process['sampledPeakPrivateBytes']) / 2**20,
        path=str(stage)))
args.output.parent.mkdir(parents=True, exist_ok=True)
with args.output.open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
for row in rows:
    print(f"{row['label']}: {row['status']}, analysis={row['analysis_seconds']}s, "
          f"wall={float(row['process_seconds']):.2f}s, "
          f"peak={row['peak_resident_mib']:.2f} MiB")
