"""Bounded phase probes using the existing canonical glma_mixed build.

Run sequentially; record actual loaded DLLs, hashes, CPU, memory and timeout.
No package installation, numerical defaults, or machine settings are changed.
The native CLI allows 1000 sweeps versus the R wrapper's 100; compare only
accepted evaluations with <=100 sweeps. Other explicit controls match R.
"""
import argparse
import ctypes as C
from ctypes import wintypes as W
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_bed(root):
    """Copy retained physical marker blocks verbatim, retaining all 5000 rows."""
    bed, bim, fam = map(Path, (root / "frozen.paths").read_text().splitlines())
    markers = (root / "frozen.markers").read_text().splitlines()
    rows = [line.split() for line in bim.read_text().splitlines()]
    index = {r[1]: j for j, r in enumerate(rows)}
    assert len(markers) == 37991 and len(index) == 50000
    assert len(fam.read_text().splitlines()) == 5000
    target = root / "frozen"
    if target.with_suffix(".bed").exists():
        raise RuntimeError("Refusing to overwrite diagnostic input")
    with bed.open("rb") as src, target.with_suffix(".bed").open("wb") as dst:
        assert src.read(3) == bytes((0x6c, 0x1b, 1))
        dst.write(bytes((0x6c, 0x1b, 1)))
        for marker in markers:
            src.seek(3 + 1250 * index[marker])
            block = src.read(1250)
            assert len(block) == 1250
            dst.write(block)
    target.with_suffix(".bim").write_text(
        "".join("\t".join(rows[index[m]]) + "\n" for m in markers))
    shutil.copy2(fam, target.with_suffix(".fam"))
    (root / "input-hashes.json").write_text(json.dumps(
        {str(p.resolve()): digest(p) for p in [bed, bim, fam,
         target.with_suffix(".bed"), target.with_suffix(".bim"),
         target.with_suffix(".fam"), root / "frozen.phen"]}, indent=2))


def probe(root, label, prefix, backend, timeout, chunk, gram, threads):
    directory = root / label
    directory.mkdir()
    runtime = directory / "runtime"
    runtime.mkdir()
    shutil.copy2("build/glma/glma_mixed.exe", runtime / "glma_mixed.exe")
    native = Path("build/task-packages/prefix" if backend == "portable"
                  else "build/prefix")
    for dll in (native / "bin").glob("*.dll"):
        shutil.copy2(dll, runtime / dll.name)
    # Use the newly rebuilt glma for MKL; use installed production glma for
    # portable. These expose the same public ABI and current source.
    if backend == "mkl":
        shutil.copy2("build/glma/libglma.dll", runtime / "libglma.dll")
    env = os.environ.copy()
    env["PATH"] = ";".join([str(runtime.resolve()),
        "C:/rtools44/x86_64-w64-mingw32.static.posix/bin",
        "C:/rtools44/usr/bin",
        "C:/Program Files (x86)/Intel/oneAPI/mkl/2025.0/bin",
        "C:/Program Files (x86)/Intel/oneAPI/compiler/2026.1/bin",
        env["PATH"]])
    env.update(MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               OMP_NUM_THREADS=str(threads), LC_ALL="C", LANG="C")
    command = [str((runtime / "glma_mixed.exe").resolve()),
        str(prefix.resolve()), str(prefix.resolve()) + ".phen",
        str((directory / "fit").resolve()), str(threads), "association",
        "15", "1024", str(gram), str(chunk), "3", "1e-7"]
    psapi = C.WinDLL("psapi", use_last_error=True)
    kernel = C.WinDLL("kernel32", use_last_error=True)
    class Memory(C.Structure):
        _fields_ = [("cb", W.DWORD), ("faults", W.DWORD)] + [
            (name, C.c_size_t) for name in ("peak_ws", "ws", "peak_paged",
            "paged", "peak_nonpaged", "nonpaged", "pagefile", "peak_pagefile",
            "private")]
    psapi.GetProcessMemoryInfo.argtypes = [W.HANDLE, C.POINTER(Memory), W.DWORD]
    psapi.EnumProcessModulesEx.argtypes = [W.HANDLE, C.POINTER(W.HMODULE),
                                          W.DWORD, C.POINTER(W.DWORD), W.DWORD]
    psapi.GetModuleFileNameExW.argtypes = [W.HANDLE, W.HMODULE, W.LPWSTR, W.DWORD]
    kernel.GetProcessTimes.argtypes = [W.HANDLE] + [C.POINTER(W.FILETIME)] * 4
    modules, peak_ws, peak_private = set(), 0, 0
    with (directory / "stdout.log").open("w") as out, \
         (directory / "stderr.log").open("w") as err:
        start = time.perf_counter()
        proc = subprocess.Popen(command, env=env, stdout=out, stderr=err,
                                creationflags=subprocess.CREATE_NO_WINDOW)
        expired = False
        while True:
            handle = W.HANDLE(int(proc._handle))
            memory = Memory(); memory.cb = C.sizeof(memory)
            if psapi.GetProcessMemoryInfo(handle, C.byref(memory), memory.cb):
                peak_ws = max(peak_ws, memory.peak_ws)
                peak_private = max(peak_private, memory.private)
            array = (W.HMODULE * 512)(); needed = W.DWORD()
            if psapi.EnumProcessModulesEx(handle, array, C.sizeof(array),
                                         C.byref(needed), 3):
                for module in array[:min(512, needed.value // C.sizeof(W.HMODULE))]:
                    name = C.create_unicode_buffer(32768)
                    if psapi.GetModuleFileNameExW(handle, module, name, len(name)):
                        modules.add(name.value)
            if proc.poll() is not None:
                break
            if time.perf_counter() - start > timeout:
                expired = True
                proc.kill(); proc.wait()
                break
            time.sleep(.2)
        ft = [W.FILETIME() for _ in range(4)]
        cpu = None
        if kernel.GetProcessTimes(handle, *[C.byref(t) for t in ft]):
            cpu = sum((t.dwHighDateTime << 32) + t.dwLowDateTime
                      for t in ft[2:]) * 1e-7
        record = dict(label=label, backend=backend, command=command,
            timeout=expired, exit_code=proc.returncode,
            wall_seconds=time.perf_counter()-start, cpu_seconds=cpu,
            peak_resident_mib=peak_ws/2**20, sampled_private_mib=peak_private/2**20,
            loaded_modules=sorted(modules),
            dll_hashes={p.name: digest(p) for p in runtime.glob("*.dll")})
        (directory / "resources.json").write_text(json.dumps(record, indent=2))
        print(json.dumps({k: v for k, v in record.items()
                          if k not in ("loaded_modules", "dll_hashes", "command")}),
              flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("label")
    parser.add_argument("--fixture", choices=("diagnostic", "frozen"), default="diagnostic")
    parser.add_argument("--backend", choices=("portable", "mkl"), default="portable")
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--chunk", type=int, default=256)
    parser.add_argument("--gram", type=int, default=64)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    if args.fixture == "frozen":
        if not (args.root / "frozen.bed").exists():
            frozen_bed(args.root)
        prefix = args.root / "frozen"
    else:
        prefix = Path("build/glma/scaling/diagnostic")
    probe(args.root, args.label, prefix, args.backend,
          args.timeout, args.chunk, args.gram, args.threads)
