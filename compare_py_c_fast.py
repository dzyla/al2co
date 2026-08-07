"""Compare fast Python vs C AL2CO implementations"""
import sys
import os
import subprocess
import tempfile
from Bio import AlignIO

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'al2co_py'))
from al2co_numpy_fast import run_al2co_python_fast

def load_alignment(filename):
    msa = AlignIO.read(filename, "clustal")
    return [str(seq.seq) for seq in msa]

def run_c_implementation(alignment, method, weight_method):
    with tempfile.NamedTemporaryFile(mode='w', suffix='.aln', delete=False) as f:
        f.write("CLUSTAL W multiple sequence alignment\n\n")
        for i, seq in enumerate(alignment):
            if i > 0:
                f.write("\n")
            f.write(f"seq{i+1:04d}  {seq}\n")
        temp_file = f.name
    
    c_exe = os.path.join(os.path.dirname(__file__), 'al2co_v4')
    if not os.path.exists(c_exe):
        c_exe = os.path.join(os.path.dirname(__file__), 'al2co_avx512')
    
    if not os.path.exists(c_exe):
        return None
    
    temp_out = temp_file + '.out'
    cmd = [c_exe, '-i', temp_file, '-o', temp_out, '-f', str(weight_method), '-c', str(method)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        return None
    
    scores = {}
    with open(temp_out, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('*') or line.startswith('M:') or \
               line.startswith('al2co') or line.startswith('Input') or \
               line.startswith('Output') or line.startswith('Conservation') or \
               line.startswith('Gap') or line.startswith('Window') or \
               line.startswith('Weighting'):
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit():
                scores[int(parts[0])] = float(parts[2])
    
    os.remove(temp_file)
    if os.path.exists(temp_out):
        os.remove(temp_out)
    
    return scores

def main():
    print("=" * 70)
    print("Fast Python vs C AL2CO Comparison")
    print("=" * 70)
    
    alignment = load_alignment('clustal.aln')
    print(f"\nAlignment: {len(alignment)} sequences, {len(alignment[0])} columns\n")
    
    test_cases = [
        (0, 0, "Entropy + Unweighted"),
        (0, 2, "Entropy + Independent"),
        (1, 0, "Variance + Unweighted"),
        (1, 2, "Variance + Independent"),
    ]
    
    all_passed = True
    for method, weight_method, desc in test_cases:
        py_scores = run_al2co_python_fast(alignment, method=method, weight_method=weight_method)
        c_scores = run_c_implementation(alignment, method, weight_method)
        
        if c_scores is None:
            print(f"✗ {desc}: C failed")
            all_passed = False
            continue
        
        max_diff = max(abs(py_scores[pos] - c_scores[pos]) for pos in py_scores)
        
        if max_diff < 1e-4:
            print(f"✓ {desc}: Max diff = {max_diff:.2e}")
        else:
            print(f"✗ {desc}: Max diff = {max_diff:.2e} (TOO LARGE)")
            all_passed = False
    
    print("\n" + "=" * 70)
    if all_passed:
        print("✓ ALL TESTS PASSED - Fast Python implementation is accurate!")
    else:
        print("✗ SOME TESTS FAILED")
    print("=" * 70)

if __name__ == "__main__":
    main()
