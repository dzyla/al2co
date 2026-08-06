#!/usr/bin/env python3
"""Benchmark al2co original vs new OpenMP version."""
import subprocess
import time
from pathlib import Path

AL2CO_DIR = Path("/home/dzyla/Code/al2co_show")
ALIGN_DIR = Path("/home/dzyla/Desktop/alignment_data/alignments")

def time_run(exe_path, infile, extra_args=""):
    """Time a single al2co run, return wall-clock seconds."""
    cmd = f"./{exe_path} -i {infile} -o /dev/null {extra_args}"
    start = time.perf_counter()
    subprocess.run(cmd, shell=True, cwd=str(AL2CO_DIR), capture_output=True)
    return time.perf_counter() - start

def benchmark(size_name, size_label):
    aln_dir = ALIGN_DIR / size_name
    if not aln_dir.exists():
        return
    print(f"\n{'='*60}")
    print(f"  {size_label}")
    print(f"{'='*60}")
    for fpath in sorted(aln_dir.glob("*.aln")):
        print(f"\n  {fpath.name}")
        orig_t = time_run("al2co_orig", str(fpath), "-p 1 -d 1 -n 2 -t 10")
        new_t  = time_run("al2co_new",  str(fpath), "-p 1 -d 1 -n 2 -t 10")
        speedup = orig_t / new_t if new_t > 0 else float('inf')
        print(f"    Original: {orig_t:.3f}s")
        print(f"    New:      {new_t:.3f}s")
        print(f"    Speedup:  {speedup:.2f}x")

if __name__ == "__main__":
    benchmark("300cols", "300-Column Alignments")
    benchmark("1000cols", "1000-Column Alignments")
    benchmark("3000cols", "3000-Column Alignments")
    print()
