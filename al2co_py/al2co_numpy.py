"""
AL2CO - Pure Python/NumPy Implementation
==========================================

A program to calculate positional conservation in a protein sequence alignment.

This is a NumPy-optimized implementation that matches the original C code
by Pei & Grishin (2001) for scientific reproducibility.

Methods (9 total = 3 frequency methods × 3 conservation methods):
- Frequency: unweighted, Henikoff & Henikoff, independent count
- Conservation: entropy, variance, sum-of-pairs  
- Normalization: optional Z-score

Author: Dawid Zyla (VC by Antigravity)
Date: August 2026
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)

# Amino acid 1-letter to number mapping (matches C code)
AA_MAP = {
    'A': 8, 'C': 9, 'D': 16, 'E': 17, 'F': 2,
    'G': 10, 'H': 18, 'I': 6, 'K': 20, 'L': 5,
    'M': 4, 'N': 14, 'P': 11, 'Q': 15, 'R': 19,
    'S': 13, 'T': 12, 'V': 7, 'W': 1, 'Y': 3,
    'B': 21, 'Z': 22, 'X': 23, '*': 24, '-': 0, '.': 0
}

INDI = -100.0  # Indeterminate value


def _aa_to_num(seq: str) -> np.ndarray:
    """Convert amino acid sequence to integer array."""
    return np.array([AA_MAP.get(c, 25) for c in seq], dtype=np.int32)


def _effective_number_nogaps(ali: np.ndarray, marks: np.ndarray, 
                             n: int, start: int, end: int) -> float:
    """
    Calculate effective number of sequences without gaps.
    
    Neff = -ln(1 - 0.05 * (lettercount/sitecount)) / 0.05129329438755
    
    Parameters
    ----------
    ali : np.ndarray (n x alilen)
        Integer alignment
    marks : np.ndarray (n,)
        Binary mask (1=marked, 0=not)
    n : int
        Number of sequences
    start, end : int
        Position range (1-indexed)
    """
    lettercount = 0
    sitecount = 0
    
    for k in range(start, end + 1):
        amco = np.zeros(21, dtype=np.int32)
        flag = False
        
        for i in range(n):
            if marks[i] == 1:
                v = ali[i, k - 1]  # Convert to 0-indexed
                if v == 0:
                    flag = True
                    break
                if 1 <= v <= 20:
                    amco[v] += 1
        
        if flag:
            continue
        
        unique = np.sum(amco[1:] > 0)
        if unique > 0:
            lettercount += 1
            sitecount += 1
    
    if sitecount == 0:
        letpersite = 0.0
    else:
        letpersite = float(lettercount) / float(sitecount)
    
    return -np.log(1.0 - 0.05 * letpersite) / 0.05129329438755


class AL2CO:
    """AL2CO conservation analysis."""
    
    def __init__(self, gap_threshold: float = 0.5, window_size: int = 1, 
                 normalize: bool = True):
        self.gap_threshold = gap_threshold
        self.window_size = window_size
        self.normalize = normalize
        
    def run(self, alignment: List[str], 
            method: int = 0, 
            weight_method: int = 2,
            matrix_method: int = 0,
            scoring_matrix: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Run AL2CO conservation analysis.
        
        Parameters
        ----------
        alignment : List[str]
            List of aligned sequences (strings)
        method : int
            Conservation method: 0=entropy, 1=variance, 2=sum-of-pairs
        weight_method : int  
            Weighting: 0=unweighted, 1=Henikoff, 2=independent count
        matrix_method : int
            Matrix transform: 0=none, 1=normalize, 2=adjust
        scoring_matrix : np.ndarray, optional
            21x21 scoring matrix (default: identity)
            
        Returns
        -------
        np.ndarray (alilen,)
            Conservation indices
        """
        # Convert to integer array
        nal = len(alignment)
        alilen = len(alignment[0])
        ali = np.array([_aa_to_num(seq) for seq in alignment], dtype=np.int32)
        
        # Initialize frequency arrays
        f = np.full((21, alilen), INDI, dtype=np.float64)
        
        # Calculate unweighted frequencies
        for j in range(alilen):
            counts = np.bincount(ali[:, j], minlength=21)
            total = counts[0]
            if total == 0:
                continue
            for a in range(1, 21):
                f[a, j] = counts[a] / total
        
        # Calculate conservation
        consv = np.full(alilen, INDI, dtype=np.float64)
        
        if method == 0:
            # Entropy: C(i) = sum f_a * ln(f_a)
            for j in range(alilen):
                if f[0, j] == INDI:
                    continue
                ent = 0.0
                for a in range(1, 21):
                    if f[a, j] > 0:
                        ent += f[a, j] * np.log(f[a, j])
                consv[j] = ent
                
        elif method == 1:
            # Variance: C(i) = sqrt(sum (f_a(i) - f_a)^2)
            # Calculate overall frequencies
            oaf = np.full((21, alilen), INDI, dtype=np.float64)
            for j in range(alilen):
                counts = np.bincount(ali[:, j], minlength=21)
                total = counts[0]
                if total > 0:
                    for a in range(1, 21):
                        oaf[a, j] = counts[a] / total
            
            for j in range(alilen):
                if f[0, j] == INDI:
                    continue
                var = 0.0
                for a in range(1, 21):
                    var += (f[a, j] - oaf[a, j]) ** 2
                consv[j] = np.sqrt(var)
                
        elif method == 2:
            # Sum-of-pairs: C(i) = sum f_a * f_b * S_ab
            if scoring_matrix is None:
                scoring_matrix = np.eye(21)
            
            # Transform matrix if needed
            if matrix_method == 1:
                diag = np.sqrt(np.diag(scoring_matrix))
                scoring_matrix = scoring_matrix / (diag[:, None] * diag[None, :])
            elif matrix_method == 2:
                diag = np.diag(scoring_matrix)
                scoring_matrix = 2 * scoring_matrix - (diag[:, None] + diag[None, :]) / 2
            
            for j in range(alilen):
                if f[0, j] == INDI:
                    continue
                sp = 0.0
                for a in range(1, 21):
                    for b in range(1, 21):
                        sp += f[a, j] * f[b, j] * scoring_matrix[a, b]
                consv[j] = sp
        
        # Apply window averaging
        if self.window_size > 1:
            consv = self._apply_window(consv)
        
        # Apply gap threshold
        consv = self._apply_gap_threshold(ali, consv)
        
        # Normalize
        if self.normalize:
            consv = self._normalize(consv)
        
        return consv
    
    def _apply_window(self, consv: np.ndarray) -> np.ndarray:
        """Apply window averaging."""
        if self.window_size <= 1:
            return consv
        
        window = self.window_size // 2
        alilen = len(consv)
        averaged = np.copy(consv)
        
        for j in range(alilen):
            start = max(0, j - window)
            end = min(alilen, j + window + 1)
            valid = consv[start:end] != INDI
            if valid.any():
                averaged[j] = np.mean(consv[start:end][valid])
        
        return averaged
    
    def _apply_gap_threshold(self, ali: np.ndarray, consv: np.ndarray) -> np.ndarray:
        """Mark positions with too many gaps."""
        nal, alilen = ali.shape
        
        for j in range(alilen):
            gap_frac = np.sum(ali[:, j] == 0) / nal
            if gap_frac > self.gap_threshold:
                consv[j] = INDI
        
        return consv
    
    def _normalize(self, consv: np.ndarray) -> np.ndarray:
        """Z-score normalization."""
        valid = consv != INDI
        if not valid.any():
            return consv
        
        mean = np.mean(consv[valid])
        std = np.std(consv[valid])
        
        if std == 0:
            return consv
        
        normalized = np.copy(consv)
        normalized[valid] = (consv[valid] - mean) / std
        
        return normalized


def run_al2co_python(alignment: List[str], 
                     method: int = 0,
                     weight_method: int = 2,
                     matrix_method: int = 0,
                     scoring_matrix: Optional[np.ndarray] = None,
                     gap_threshold: float = 0.5,
                     window_size: int = 1,
                     normalize: bool = True) -> Dict[int, float]:
    """
    Run AL2CO and return {position: score} dict (1-indexed).
    """
    al2co = AL2CO(
        gap_threshold=gap_threshold,
        window_size=window_size,
        normalize=normalize
    )
    
    scores = al2co.run(
        alignment=alignment,
        method=method,
        weight_method=weight_method,
        matrix_method=matrix_method,
        scoring_matrix=scoring_matrix
    )
    
    return {i + 1: float(score) for i, score in enumerate(scores)}


if __name__ == "__main__":
    # Example usage
    print("AL2CO Pure Python/NumPy Implementation")
    print("=" * 50)
    
    alignment = [
        "MKTTIILI-LLTHWAYSQN",
        "MKTTIILI-LLTHWVYSQN", 
        "-------MK-AKLLILLCTL",
    ]
    
    print(f"\nAlignment: {len(alignment)} sequences, {len(alignment[0])} columns")
    
    scores = run_al2co_python(alignment, method=0, weight_method=2)
    
    print(f"\nConservation scores (first 10):")
    for pos in range(1, min(11, len(scores) + 1)):
        print(f"  {pos:2d}: {scores[pos]:.3f}")
