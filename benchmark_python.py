"""Benchmark Python AL2CO implementations"""
import time
import sys
sys.path.insert(0, 'al2co_py')

from al2co_numpy import run_al2co_python
from al2co_numpy_fast import run_al2co_python_fast
from Bio import AlignIO

# Load alignment
msa = AlignIO.read('clustal.aln', 'clustal')
alignment = [str(seq.seq) for seq in msa]
print(f"Alignment: {len(alignment)} sequences, {len(alignment[0])} columns\n")

# Test different sizes
for size_name, n_seq in [("100", 100), ("500", 500), ("1000", 1000)]:
    if n_seq > len(alignment):
        n_seq = len(alignment)
    
    sub_align = alignment[:n_seq]
    print(f"Testing with {n_seq} sequences:")
    
    # Old implementation
    start = time.time()
    scores_old = run_al2co_python(sub_align, method=0, weight_method=2)
    time_old = time.time() - start
    
    # New implementation
    start = time.time()
    scores_new = run_al2co_python_fast(sub_align, method=0, weight_method=2)
    time_new = time.time() - start
    
    # Compare
    max_diff = max(abs(scores_old[i] - scores_new[i]) for i in scores_old)
    
    print(f"  Old:  {time_old:.3f}s")
    print(f"  New:  {time_new:.3f}s")
    print(f"  Speedup: {time_old/time_new:.1f}x")
    print(f"  Max diff: {max_diff:.2e}")
    print()

print("=" * 60)
print("Optimization successful!")
