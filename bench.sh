#!/bin/bash
# Benchmark script for al2co comparison
# Usage: bash bench.sh [300cols|1000cols|3000cols|all]

set -e
cd /home/dzyla/Code/al2co_show

SIZE="${1:-all}"

echo "=========================================="
echo "AL2CO Speed Benchmark"
echo "=========================================="
echo ""

if [ "$SIZE" = "300cols" ] || [ "$SIZE" = "all" ]; then
  echo "--- 300-column alignments ---"
  for f in /home/dzyla/Desktop/alignment_data/alignments/300cols/*.aln; do
    BASE=$(basename "$f")
    echo ""
    echo "File: $BASE"
    T1=$( { time ./al2co_orig -i "$f" -o /dev/null -p 1 -d 1 -n 2 -t 10 > /dev/null 2>&1; } 2>&1 | grep real | awk '{print $2}' )
    T2=$( { time ./al2co_new -i "$f" -o /dev/null -p 1 -d 1 -n 2 -t 10 > /dev/null 2>&1; } 2>&1 | grep real | awk '{print $2}' )
    echo "  Original: $T1"
    echo "  New:      $T2"
  done
fi

if [ "$SIZE" = "1000cols" ] || [ "$SIZE" = "all" ]; then
  echo ""
  echo "--- 1000-column alignments ---"
  for f in /home/dzyla/Desktop/alignment_data/alignments/1000cols/*.aln; do
    BASE=$(basename "$f")
    echo ""
    echo "File: $BASE"
    T1=$( { time ./al2co_orig -i "$f" -o /dev/null -p 1 -d 1 -n 2 -t 10 > /dev/null 2>&1; } 2>&1 | grep real | awk '{print $2}' )
    T2=$( { time ./al2co_new -i "$f" -o /dev/null -p 1 -d 1 -n 2 -t 10 > /dev/null 2>&1; } 2>&1 | grep real | awk '{print $2}' )
    echo "  Original: $T1"
    echo "  New:      $T2"
  done
fi

if [ "$SIZE" = "3000cols" ] || [ "$SIZE" = "all" ]; then
  echo ""
  echo "--- 3000-column alignments ---"
  for f in /home/dzyla/Desktop/alignment_data/alignments/3000cols/*.aln; do
    BASE=$(basename "$f")
    echo ""
    echo "File: $BASE"
    T1=$( { time ./al2co_orig -i "$f" -o /dev/null -p 1 -d 1 -n 2 -t 10 > /dev/null 2>&1; } 2>&1 | grep real | awk '{print $2}' )
    T2=$( { time ./al2co_new -i "$f" -o /dev/null -p 1 -d 1 -n 2 -t 10 > /dev/null 2>&1; } 2>&1 | grep real | awk '{print $2}' )
    echo "  Original: $T1"
    echo "  New:      $T2"
  done
fi
