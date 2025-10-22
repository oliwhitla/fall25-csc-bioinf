# Week 4 - Sequence Alignment

## Overview

Project: Implementing and benchmarking four sequence-alignment algorithms:


- Global (Needleman–Wunsch)

- Local (Smith–Waterman)

- Semi-Global (Fitting alignment)

- Affine-gap Global alignment



## Implementation Notes



- I used a single Align class containing all four algorithms

- Each method builds and fills a dynamic-programming (DP) matrix, returning the final alignment score only

- The Codon conversion was much easier this week – it compiled without type errors or Optional-handling issues

- I rounded runtimes to two decimal places instead of integers, since most runs complete within a few milliseconds; rounding to whole numbers produced mostly zeros
  


## Problems

### Why Affine Was Skipped for MT_human



The affine-gap version maintains three DP matrices

For large sequences (~16 500 bp × 16 500 bp), this triples both time and memory

The GitHub CLI runner has limited memory and running affine alignment on MT_human × MT_orang therefore caused out-of-memory or timeout errors (exit code 143)

I ran only the Global, Local, and Semi-Global alignments for the MT_human and MT_orang pair,
and skipped the Affine version for that case



## Time spent

~ day, 8hrs 




