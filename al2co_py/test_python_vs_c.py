"""
Test Python implementation vs C implementation
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from Bio import AlignIO
from al2co_numpy import run_al2co_python
import numpy as np

def load_alignment(filename):
    """Load alignment from Clustal file."""
    msa = AlignIO.read(filename, "clustal")
    alignment = [str(seq.seq) for seq in msa]
    return alignment

def run_c_implementation(alignment, method, weight_method):
    """Run C implementation and return scores."""
    import subprocess
    import tempfile
    
    # Write alignment to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.aln', delete=False) as f:
        # Write Clustal format
        f.write("CLUSTAL W multiple sequence alignment\n\n")
        max_len = max(len(seq) for seq in alignment)
        for seq in alignment:
            f.write(f"{seq:<{max_len}}\n")
        temp_file = f.name
    
    # Run C implementation
    c_exe = os.path.join(os.path.dirname(__file__), '..', 'al2co_v4')
    if not os.path.exists(c_exe):
        c_exe = os.path.join(os.path.dirname(__file__), '..', 'al2co_avx512')
    
    temp_out = temp_file + '.out'
    
    cmd = [c_exe, '-i', temp_file, '-o', temp_out, 
           '-f', str(weight_method), '-c', str(method)]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"C implementation error: {result.stderr}")
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
        return False, "One or both implementations failed"
    
    # Check same positions
    py_pos = set(py_scores.keys())
    c_pos = set(c_scores.keys())
    
    if py_pos != c_pos:
        return False, f"Position mismatch: Python has {len(py_pos)}, C has {len(c_pos)}"
    
    # Compare scores
    max_diff = 0
    for pos in py_pos:
        diff = abs(py_scores[pos] - c_scores[pos])
        max_diff = max(max_diff, diff)
        if diff > tolerance:
            return False, f"Score mismatch at position {pos}: Python={py_scores[pos]:.6f}, C={c_scores[pos]:.6f}, diff={diff:.6f}"
    
    return True, f"Max difference: {max_diff:.8f}"

def main():
    print("Testing Python vs C AL2CO Implementation")
    print("=" * 60)
    
    # Test with small alignment
    print("\n1. Testing with clustal.aln...")
    alignment = load_alignment(os.path.join(os.path.dirname(__file__), '..', 'clustal.aln'))
    print(f"   Alignment: {len(alignment)} sequences, {len(alignment[0])} columns")
    
    # Test different methods
    for method in [0, 1, 2]:
        for weight_method in [0, 1, 2]:
            print(f"\n   Testing method={method}, weight_method={weight_method}...")
            
            # Python implementation
            py_scores = run_al2co_python(alignment, method=method, weight_method=weight_method)
            
            # C implementation
            c_scores = run_c_implementation(alignment, method, weight_method)
            
            # Compare
            success, message = compare_scores(py_scores, c_scores)
            
            if success:
                print(f"   ✓ PASS: {message}")
            else:
                print(f"   ✗ FAIL: {message}")
    
    # Test with larger alignment
    print("\n2. Testing with msa_out.aln (11K sequences)...")
    try:
        alignment_large = load_alignment(os.path.join(os.path.dirname(__file__), '..', 'msa_out.aln'))
        print(f"   Alignment: {len(alignment_large)} sequences")
        
        # Python implementation
        py_scores_large = run_al2co_python(alignment_large, method=0, weight_method=2)
        print(f"   Python completed: {len(py_scores_large)} positions")
        
        # C implementation (only if MAXSEQNUM >= 11000)
        c_scores_large = run_c_implementation(alignment_large, 0, 2)
        
        if c_scores_large is not None:
            success, message = compare_scores(py_scores_large, c_scores_large)
            if success:
                print(f"   ✓ PASS: {message}")
            else:
                print(f"   ✗ FAIL: {message}")
        else:
            print("   C implementation skipped (MAXSEQNUM limit)")
            
    except Exception as e:
        print(f"   Error: {e}")
    
    print("\n" + "=" * 60)
    print("Test completed!")

if __name__ == "__main__":
    main()
