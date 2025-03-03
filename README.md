# MAFFT Implemenration with MPI Parallelization

## Overview

This repository contains an implementation of Multiple Sequence Alignment (MSA) using both serial and MPI-based parallel approaches. The project aims to compare the performance and scalability between a traditional serial MSA algorithm and an optimized parallel version using MPI (Message Passing Interface).

The approach is based on the MAFFT algorithm and involves multiple stages including pairwise sequence alignment, guide tree construction, and progressive alignment.

The parallel version of the algorithm utilizes MPI to distribute computational tasks such as FFT (Fast Fourier Transform) application, pairwise alignments, and guide tree construction across multiple processes, significantly speeding up the alignment process for large datasets.

## Features

- **Serial MSA Approach**:
  - Pairwise sequence alignment using FFT.
  - Guide tree construction using UPGMA or Neighbor-Joining.
  - Progressive alignment and refinement.
  
- **MPI Parallel MSA Approach**:
  - Parallelized FFT and pairwise alignment computations.
  - Distributed guide tree construction.
  - Optimized inter-process communication and synchronization.
  
- **Benchmarking**:
  - Performance comparison between serial and parallel approaches.
  - Evaluation of speedup, scalability, and memory usage.
