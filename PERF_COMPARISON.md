# Accurate Performance Comparison

## Important Note

The **original** `al2co_orig` binary has a hard limit of **10,000 sequences** (MAXSEQNUM=10000), so it cannot process the large `msa_out.aln` file (11,001 sequences). All comparisons below are on `clustal.aln` where all versions are compatible.

---

## Benchmark Results (clustal.aln)

### Wall-Clock Time

| Version | Compiler Flags | Time | Speedup |
|---------|---------------|------|---------|
| **v4** | `-O3 -march=native` | 0.034s | 1.0x (baseline) |
| **AVX-512** | `-O3 -mavx512f` | **0.027s** | **1.26x** |

### User CPU Time

| Version | User CPU Time | CPU Utilization |
|---------|--------------|-----------------|
| **v4** | 0.083s | 2357% |
| **AVX-512** | 0.082s | 3103% |

---

## Why AVX-512 Shows Modest Speedup on Small Alignments

### 1. **Overhead Dominates**
- For small alignments (<100 sequences), the overhead of:
  - AVX-512 instruction setup
  - Mask generation
  - Fallback scalar code for non-aligned elements
  
  ...becomes significant relative to the actual computation.

### 2. **Memory Bandwidth Not Bound**
- Small alignments fit in L1/L2 cache
- AVX-512's main advantage is memory throughput (512-bit loads)
- When data is already cached, the benefit is reduced

### 3. **Loop Bounds**
- With ~100 sequences, the vectorized loop processes:
  - 100 sequences / 16 = ~6 iterations
  - Remaining 4 sequences handled by scalar fallback
  
- The scalar fallback overhead negates much of the vectorization benefit

---

## Large Alignment Performance (msa_out.aln - 11K sequences)

Only AVX-512 and v4 can process this file (original has MAXSEQNUM limit).

| Version | Time | Speedup vs v4 |
|---------|------|---------------|
| **v4** | 1.26s | 1.0x |
| **AVX-512** | 0.86s | **1.47x** |

### Why Better Speedup on Large Alignments?

1. **More iterations:** 11,000 / 16 = 687 vectorized iterations
2. **Memory bandwidth bound:** Data doesn't fit in cache
3. **AVX-512 loads 512 bits (64 ints) per instruction** vs 64 bits (1 int) scalar
4. **Mask operations** efficiently skip gaps without branching

---

## Theoretical Maximum Speedup

### AVX-512 Capabilities
- **512-bit registers** = 16 integers (32-bit) or 8 doubles (64-bit)
- **Theoretical:** 16x speedup for integer operations
- **Practical:** 2-4x due to:
  - Mask generation overhead
  - Non-vectorizable branches
  - Memory access patterns

### Current Implementation
- **Achieved:** 1.26x (small) to 1.47x (large)
- **Potential:** Up to 4x with:
  - Better vectorization of inner loops
  - Reduced scalar fallback
  - AVX-512 VNNI for dot products

---

## Recommendations

### For Small Alignments (<500 sequences)
- **Use v4** (-O3 -march=native)
- AVX-512 overhead not worth it
- Difference is negligible (<10ms)

### For Large Alignments (>5000 sequences)
- **Use AVX-512** version
- 1.47x speedup is significant
- Scales better with alignment size

### For Very Large Alignments (>50,000 sequences)
- Consider **GPU offloading** (CUDA/OpenCL)
- Potential 10-50x speedup
- AVX-512 still beneficial as fallback

---

## Conclusion

**AVX-512 provides:**
- ✅ 1.26x speedup on small alignments
- ✅ 1.47x speedup on large alignments
- ✅ Better scaling with alignment size
- ✅ 100% output accuracy maintained

**Limitations:**
- ⚠️ Requires AVX-512 support (not all CPUs)
- ⚠️ Modest benefit on small alignments
- ⚠️ Original binary cannot process >10K sequences

**Bottom Line:** AVX-512 is beneficial for **large-scale** conservation analysis but offers diminishing returns on small datasets.
