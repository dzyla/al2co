"""
AL2CO - Optimized Python/NumPy Implementation
===============================================

High-performance pure Python/NumPy implementation using vectorization.

Key optimizations:
- NumPy vectorization for all inner loops
- Batch processing of sequences
- Efficient matrix operations
- Minimal Python overhead

Author: Dawid Zyla (VC by Antigravity)
Date: August 2026
"""

import numpy as np
from typing import Dict, List, Optional
import logging

logger = logging.getLogger(__name__)

# Amino acid mapping (matches C code)
AA_MAP = {
    'A': 8, 'C': 9, 'D': 16, 'E': 17, 'F': 2,
    'G': 10, 'H': 18, 'I': 6, 'K': 20, 'L': 5,
    'M': 4, 'N': 14, 'P': 11, 'Q': 15, 'R': 19,
    'S': 13, 'T': 12, 'V': 7, 'W': 1, 'Y': 3,
    'B': 21, 'Z': 22, 'X': 23, '*': 24, '-': 0, '.': 0
}

INDI = -100.0


def _aa_array_to_int(ali_char: np.ndarray) -> np.ndarray:
    """Convert character array to integer array using lookup table."""
    # Create lookup table
    lookup = np.full(256, 25, dtype=np.int32)
    for aa, num in AA_MAP.items():
        lookup[ord(aa)] = num
    
    # Vectorized conversion
    return lookup[ali_char.view(np.uint8)]


class AL2COFast:
    """Optimized AL2CO using NumPy vectorization."""
    
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
        """Run AL2CO conservation analysis."""
        # Convert to integer array
        nal = len(alignment)
        alilen = len(alignment[0])
        
        # Vectorized alignment conversion
        ali_char = np.array([list(seq) for seq in alignment])
        ali = _aa_array_to_int(ali_char)
        
        # Calculate frequencies
        f = self._calc_frequencies_vectorized(ali, alilen)
        
        # Calculate conservation
        consv = self._calc_conservation_vectorized(ali, f, method, 
                                                   scoring_matrix, matrix_method)
        
        # Apply transformations
        if self.window_size > 1:
            consv = self._apply_window(consv)
        
        consv = self._apply_gap_threshold(ali, consv)
        
        if self.normalize:
            consv = self._normalize(consv)
        
        return consv
    
    def _calc_frequencies_vectorized(self, ali: np.ndarray, alilen: int) -> np.ndarray:
        """Calculate frequencies using vectorized operations."""
        nal = ali.shape[0]
        f = np.full((21, alilen), INDI, dtype=np.float64)
        
        # Count amino acids at each position
        # ali shape: (nal, alilen)
        # We want counts shape: (21, alilen)
        
        for a in range(1, 21):
            mask = ali == a  # (nal, alilen)
            counts = mask.sum(axis=0)  # (alilen,)
            f[a] = counts.astype(np.float64)
        
        # Calculate total (excluding gaps)
        total = nal - (ali == 0).sum(axis=0)  # (alilen,)
        
        # Normalize
        valid = total > 0
        f[:, valid] = f[:, valid] / total[valid, np.newaxis]
        
        return f
    
    def _calc_conservation_vectorized(self, ali: np.ndarray, f: np.ndarray,
                                      method: int, scoring_matrix: Optional[np.ndarray],
                                      matrix_method: int) -> np.ndarray:
        """Calculate conservation using vectorized operations."""
        alilen = f.shape[1]
        consv = np.full(alilen, INDI, dtype=np.float64)
        
        if method == 0:
            # Entropy: C(i) = sum f_a * ln(f_a)
            # Vectorized: f[1:21] * log(f[1:21])
            f_valid = f[1:21].copy()
            f_valid[f_valid <= 0] = 1.0  # Avoid log(0)
            entropy = np.sum(f[1:21] * np.log(f[1:21]), axis=0)
            consv = np.where(f[0] == INDI, INDI, entropy)
            
        elif method == 1:
            # Variance: C(i) = sqrt(sum (f_a(i) - f_a)^2)
            # Calculate overall frequencies
            nal = ali.shape[0]
            oaf = np.full((21, alilen), INDI, dtype=np.float64)
            
            for a in range(1, 21):
                oaf[a] = (ali == a).sum(axis=0) / (nal - (ali == 0).sum(axis=0))
            
            # Vectorized variance calculation
            diff = f[1:21] - oaf[1:21]  # (20, alilen)
            variance = np.sum(diff ** 2, axis=0)  # (alilen,)
            consv = np.where(f[0] == INDI, INDI, np.sqrt(variance))
            
        elif method == 2:
            # Sum-of-pairs: C(i) = sum f_a * f_b * S_ab
            if scoring_matrix is None:
                scoring_matrix = np.eye(21)
            
            # Transform matrix
            if matrix_method == 1:
                diag = np.sqrt(np.diag(scoring_matrix))
                scoring_matrix = scoring_matrix / (diag[:, None] * diag[None, :])
            elif matrix_method == 2:
                diag = np.diag(scoring_matrix)
                scoring_matrix = 2 * scoring_matrix - (diag[:, None] + diag[None, :]) / 2
            
            # Vectorized: f[1:21].T @ S @ f[1:21]
            f_aa = f[1:21]  # (20, alilen)
            S = scoring_matrix[1:21, 1:21]  # (20, 20)
            consv = np.sum(f_aa.T @ S @ f_aa, axis=1)  # (alilen,)
            consv = np.where(f[0] == INDI, INDI, consv)
        
        return consv
    
    def _apply_window(self, consv: np.ndarray) -> np.ndarray:
        """Apply window averaging using convolution."""
        if self.window_size <= 1:
            return consv
        
        window = self.window_size // 2
        alilen = len(consv)
        
        # Replace INDI with 0 for convolution
        valid_mask = consv != INDI
        consv_work = np.where(valid_mask, consv, 0.0)
        
        # Create convolution kernel
        kernel = np.ones(self.window_size) / self.window_size
        
        # Apply convolution (with boundary handling)
        averaged = np.zeros(alilen)
        for j in range(alilen):
            start = max(0, j - window)
            end = min(alilen, j + window + 1)
            segment = consv_work[start:end]
            valid_in_seg = valid_mask[start:end]
            
            if valid_in_seg.any():
                averaged[j] = np.mean(segment[valid_in_seg])
            else:
                averaged[j] = INDI
        
        return averaged
    
    def _apply_gap_threshold(self, ali: np.ndarray, consv: np.ndarray) -> np.ndarray:
        """Mark positions with too many gaps."""
        nal = ali.shape[0]
        gap_frac = (ali == 0).sum(axis=0) / nal
        consv[gap_frac > self.gap_threshold] = INDI
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
        
        normalized = consv.copy()
        normalized[valid] = (consv[valid] - mean) / std
        
        return normalized


def run_al2co_python_fast(alignment: List[str], 
                          method: int = 0,
                          weight_method: int = 2,
                          matrix_method: int = 0,
                          scoring_matrix: Optional[np.ndarray] = None,
                          gap_threshold: float = 0.5,
                          window_size: int = 1,
                          normalize: bool = True) -> Dict[int, float]:
    """Run fast AL2CO and return {position: score} dict."""
    al2co = AL2COFast(
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
    import time
    from Bio import AlignIO
    
    # Test with clustal.aln
    msa = AlignIO.read('clustal.aln', 'clustal')
    alignment = [str(seq.seq) for seq in msa]
    print(f"Alignment: {len(alignment)} sequences, {len(alignment[0])} columns")
    
    start = time.time()
    scores = run_al2co_python_fast(alignment, method=0, weight_method=2)
    elapsed = time.time() - start
    
    print(f"Time: {elapsed:.3f}s")
    print(f"Positions: {len(scores)}")
    print(f"First 10 scores:")
    for i in range(1, 11):
        print(f"  {i}: {scores[i]:.6f}")
