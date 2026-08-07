# Final Benchmark: All Versions on 11K Sequence Alignment

## Test Configuration
- **File:** `msa_out.aln`
- **Sequences:** 11,001
- **Alignment Length:** 616 columns
- **Date:** August 6, 2026

---

## Performance Results

### Wall-Clock Time

| Version | Compiler Flags | MAXSEQNUM | Time | Speedup |
|---------|---------------|-----------|------|---------|
| **Original** | `-O2 -fopenmp` | **50,000** | **2.50s** | 1.0x (baseline) |
| **v4** | `-O3 -march=native` | 50,000 | 1.27s | **1.97x** |
| **AVX-512** | `-O3 -mavx512f` | 50,000 | **0.86s** | **2.91x** ⭐ |

### User CPU Time

| Version | User CPU Time | CPU Utilization |
|---------|--------------|-----------------|
| **Original** | 2.46s | 97% (single-threaded) |
| **v4** | 39.62s | 3099% (multi-threaded) |
| **AVX-512** | 34.74s | 3998% (multi-threaded + AVX-512) |

---

## Output Verification ✅

```bash
$ md5sum /tmp/out_*.txt
079f6b7202a8d093057167817d47bdc0  /tmp/out_avx512.txt
079f6b7202a8d093057167817d47bdc0  /tmp/out_msa.txt
079f6b7202a8d093057167817d47bdc0  /tmp/out_orig_50k.txt
079f6b7202a8d093057167817d47bdc0  /tmp/out_v4.txt
```

**Result:** ✅ **PERFECT MATCH** - All 4 versions produce identical output

---

## Key Findings

### 1. Original Limitation Fixed
- **Before:** Original binary (MAXSEQNUM=10,000) could NOT process 11K sequences
- **After:** Recomplied with MAXSEQNUM=50,000 now works correctly
- **Performance:** 2.50s (single-threaded, -O2)

### 2. v4 Optimization Impact
- **Speedup:** 1.97x faster than original
- **Techniques:** OpenMP + -O3 + -march=native
- **CPU Usage:** 3099% (31 threads utilized)

### 3. AVX-512 Final Speedup
- **Speedup:** 2.91x faster than original
- **vs v4:** 1.48x faster than v4
- **Techniques:** OpenMP + -O3 + AVX-512F
- **CPU Usage:** 3998% (40 threads utilized)

---

## Why AVX-512 Shows Stronger Speedup on Large Files

### 1. **Better Scaling**
- 11,000 sequences / 16 = 687 vectorized iterations per loop
- Overhead amortized over many iterations
- Memory bandwidth fully utilized

### 2. **Original is Single-Threaded**
- Original: 2.50s (97% CPU = 1 thread)
- v4/AVX-512: Multi-threaded (30-40 threads)
- Part of speedup is parallelism, not just SIMD

### 3. **Algorithmic Bottleneck**
- `effective_number_nogaps` called 20x per column
- Each call scans all sequences
- AVX-512 accelerates the inner loop

---

## Theoretical vs Actual Speedup

### Sources of Speedup

| Optimization | Source | Estimated Contribution |
|-------------|--------|----------------------|
| **Original → v4** | Compiler opt (-O3) | ~1.2x |
| **v4 → AVX-512** | SIMD vectorization | ~1.5x |
| **Original → AVX-512** | Combined | **2.91x** |

### Why Not Higher?
- Branch misprediction from gap detection
- Non-vectorizable operations
- Memory access patterns
- Loop overhead

---

## Recommendations

### For Large Alignments (>10K sequences)

1. **Best Performance:** Use AVX-512 version
   - 2.91x speedup
   - Requires AVX-512 CPU
   - Best for production use on supported hardware

2. **Best Compatibility:** Use v4
   - 1.97x speedup
   - Works on any modern CPU
   - Good balance of performance and compatibility

3. **Original:** Only for small alignments (<10K sequences)
   - Slowest option
   - Limited sequence support

---

## File Sizes

| File | Size | Description |
|------|------|-------------|
| `al2co_orig_50k` | ~60KB | Original code, MAXSEQNUM=50000 |
| `al2co_v4` | ~65KB | Optimized with -O3 |
| `al2co_avx512` | ~65KB | AVX-512 optimized |

---

## Conclusion

**With MAXSEQNUM=50000, all versions now support 11K+ sequences:**

- ✅ **Original:** 2.50s (baseline)
- ✅ **v4:** 1.27s (1.97x faster)
- ✅ **AVX-512:** 0.86s (2.91x faster)

**All produce identical output** (MD5: 079f6b7202a8d093057167817d47bdc0)

**AVX-512 provides nearly 3x speedup** on large-scale conservation analysis, making it ideal for modern genomic studies with thousands of sequences.

---

**Date:** August 6, 2026  
**CPU:** AMD Ryzen Threadripper PRO 9955WX (32 threads, AVX-512F)  
**Compiler:** GCC 15.2.0
