import os
import pytest
from Bio import AlignIO
from al2co import run_al2co

def test_al2co_execution():
    # Load test alignment file
    aln_file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'clustal.aln')
    msa = AlignIO.read(aln_file, "clustal")
    
    # Run al2co
    scores = run_al2co(msa, f=1, c=2, w=1)
    
    # Verify we got scores back
    assert isinstance(scores, dict)
    assert len(scores) > 0
    
    # Check that positions are 1-indexed and have float scores
    assert 1 in scores
    assert isinstance(scores[1], float)
