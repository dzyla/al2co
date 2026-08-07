"""Profile AL2CO Python implementation"""
import time
import sys
sys.path.insert(0, 'al2co_py')
from al2co_numpy import AL2CO
from Bio import AlignIO

# Load alignment
msa = AlignIO.read('clustal.aln', 'clustal')
alignment = [str(seq.seq) for seq in msa]
print(f"Alignment: {len(alignment)} sequences, {len(alignment[0])} columns")

# Profile
al2co = AL2CO()

start = time.time()
scores = al2co.run(alignment, method=0, weight_method=2)
total = time.time() - start

print(f"\nTotal time: {total:.3f}s")
print(f"Positions: {len(scores)}")
