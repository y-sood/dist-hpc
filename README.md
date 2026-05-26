# 🌐 HPC & Distributed Systems 

A collection of numerical and sorting algorithms engineered for scalability across multi-core and distributed-memory environments. This repository focuses on parallel architecture scaling, memory profiling (Cachegrind/Memcheck), and mitigating hardware bottlenecks like NUMA topology boundaries and atomic contention.

All performance evaluations, including Weak and Strong scaling tests, were executed on the SLURM-managed UPPMAX cluster.

---

## 🧮 Project Index

### 1. [N-Body Gravitational Simulation](./1_nbody_simulation)
**Focus:** $\mathcal{O}(N^2)$ algorithm optimization and cache locality.
- **Implementation:** Simulated the gravitational evolution of a galaxy using Newton's laws of motion.
- **Optimization Strategy:** Eliminated redundant computations using Newton's Third Law and precomputed inverse masses. Re-structured data into an Array of Structs (AoS) to drastically improve spatial locality.
- **Profiling:** Validated via Valgrind/Cachegrind, achieving a **0.00% Instruction Cache miss rate** and near-zero L1/L2 cache misses by aggressively managing memory access patterns and utilizing static inline functions.

### 2. [Parallel Quicksort (Shared Memory)](./2_parallel_quicksort)
**Focus:** Divide-and-conquer task parallelism and dynamic memory allocation.
- **Implementation:** Developed a parallel Quicksort using OpenMP, testing multiple global pivot selection strategies (e.g., median of local medians) to maintain load balance.
- **Memory Management:** Engineered an adaptive partitioning buffer that dynamically switches between stack allocation (for fast, zero-overhead memory access) and heap allocation based on an evaluated `MAX_STACK_SCRATCH` threshold (4MB).
- **Scaling Insights:** Analyzed strong and weak scaling fall-offs across 24 threads, isolating the performance degradation boundary crossing between NUMA nodes on the AMD Opteron CPU architecture. 

### 3. [Stochastic Epidemic Simulation (Distributed Memory)](./3_monte_carlo_mpi)
**Focus:** Monte Carlo methods, Gillespie SSA, and MPI distributed aggregation.
- **Implementation:** Simulated malaria transmission using the Gillespie Stochastic Simulation Algorithm (SSA) wrapped in a Monte Carlo framework. Scaled the simulation across 32 physical processes using MPI.
- **Scaling Results:** Achieved massive performance gains under strong scaling, actually hitting **super-linear speedups** (36.15x on 32 processes). 
- **Hardware Insights:** Investigated how the super-linear scaling was directly tied to the problem size shrinking enough to perfectly fit inside the L1/L2 caches of the localized processors, eliminating main memory pressure and cache contention entirely.

---

## 🔬 Core Competencies Demonstrated
- **Scalability Profiling:** Evaluating weak vs. strong scaling, Amdahl's Law limits, and NUMA node topology impacts.
- **Memory Diagnostics:** Deep utilization of Valgrind, Cachegrind, and Memcheck to trace memory leaks and optimize L1/L2 cache hit rates.
- **Distributed Computing:** Managing complex task partitioning and data gathering across MPI processes and OpenMP thread pools.