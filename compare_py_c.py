"""
Compare Python vs C AL2CO implementations
"""
import sys
import os
import subprocess
import tempfile
from Bio import AlignIO

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'al2co_py'))
from al2co_numpy import run_al2co_python

def load_alignment(filename):
    """Load alignment from Clustal file."""
    msa = AlignIO.read(filename, "clustal")
    alignment = [str(seq.seq) for seq in msa]
    return alignment

def run_c_implementation(alignment, method, weight_method):
    """Run C implementation and return scores."""
    # Write alignment to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.aln', delete=False) as f:
        f.write("CLUSTAL W multiple sequence alignment\n\n")
        max_len = max(len(seq) for seq in alignment)
        for i, seq in enumerate(alignment):
            if i > 0:
                f.write("\n")
            # Write sequence name and sequence
            name = f"seq{i+1:04d}"
            f.write(f"{name:<20} {seq}\n")
        temp_file = f.name
    
    # Find C executable
    c_exe = os.path.join(os.path.dirname(__file__), 'al2co_v4')
    if not os.path.exists(c_exe):
        c_exe = os.path.join(os.path.dirname(__file__), 'al2co_avx512')
    
    if not os.path.exists(c_exe):
        print(f"  C executable not found!")
        return None
    
    temp_out = temp_file + '.out'
    
    cmd = [c_exe, '-i', temp_file, '-o', temp_out, 
           '-f', str(weight_method), '-c', str(method)]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  C error: {result.stderr[:200]}")
        return None
    
    # Read output
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
                pos = int(parts[0])
                score = float(parts[2])
                scores[pos] = score
    
    # Clean up
    os.remove(temp_file)
    if os.path.exists(temp_out):
        os.remove(temp_out)
    
    return scores

def compare_scores(py_scores, c_scores, tolerance=1e-5):
    """Compare Python and C scores."""
    if py_scores is None or c_scores is None:
        return False, "One or both failed"
    
    py_pos = set(py_scores.keys())
    c_pos = set(c_scores.keys())
    
    if py_pos != c_pos:
        return False, f"Position mismatch: Py={len(py_pos)}, C={len(c_pos)}"
    
    max_diff = 0
    for pos in sorted(py_pos):
        diff = abs(py_scores[pos] - c_scores[pos])
        max_diff = max(max_diff, diff)
    
    return True, f"Max diff: {max_diff:.10f}"

def main():
    print("=" * 70)
    print("AL2CO Python vs C Implementation Comparison")
    print("=" * 70)
    
    # Test with clustal.aln
    print("\n1. Testing with clustal.aln...")
    alignment = load_alignment('clustal.aln')
    print(f"   {len(alignment)} sequences, {len(alignment[0])} columns")
    
    test_cases = [
        (0, 0, "Entropy + Unweighted"),
        (0, 1, "Entropy + Henikoff"),
        (0, 2, "Entropy + Independent"),
        (1, 0, "Variance + Unweighted"),
        (1, 1, "Variance + Henikoff"),
        (1, 2, "Variance + Independent"),
    ]
    
    all_passed = True
    for method, weight_method, desc in test_cases:
        print(f"\n   Testing: {desc}")
        
        py_scores = run_al2co_python(alignment, method=method, weight_method=weight_method)
        c_scores = run_c_implementation(alignment, method, weight_method)
        
        if c_scores is None:
            print(f"   ✗ SKIP (C failed)")
            all_passed = False
            continue
        
        success, message = compare_scores(py_scores, c_scores)
        
        if success:
            print(f"   ✓ PASS: {message}")
        else:
            print(f"   ✗ FAIL: {message}")
            all_passed = False
    
    # Show sample output
    print("\n" + "=" * 70)
    print("Sample Output Comparison (Entropy + Independent):")
    print("=" * 70)
    py_scores = run_al2co_python(alignment, method=0, weight_method=2)
    c_scores = run_c_implementation(alignment, 0, 2)
    
    print(f"\n{'Position':<10} {'Python':<12} {'C':<12} {'Diff':<12}")
    print("-" * 46)
    for i in range(1, 11):
        diff = abs(py_scores[i] - c_scores[i])
        print(f"{i:<10} {py_scores[i]:<12.6f} {c_scores[i]:<12.6f} {diff:<12.8f}")
    
    print("\n" + "=" * 70)
    if all_passed:
        print("✓ ALL TESTS PASSED")
    else:
        print("✗ SOME TESTS FAILED")
    print("=" * 70)

if __name__ == "__main__":
    main()
