# Benchmark Results and Discussion

The solver was tested for grid sizes `N ∈ {16, 32, 64, 128, 256}` with forcing term `f(x,y) = 8π²sin(2πx)sin(2πy)` and tolerance `1e-8` (max 100,000 iterations).

---

## Raw Results

| N   | [1 MPI, 1 OMP] | [2 MPI, 1 OMP] | [4 MPI, 1 OMP] | [4 MPI, 2 OMP] | L2 Error (4 MPI)  | Iters (4 MPI) |
|-----|---------------|---------------|---------------|---------------|----------------|------------|
| 16  | 0.000433 s    | 0.000260 s    | 0.000778 s    | 0.1046 s      | 0.02856        | 180        |
| 32  | 0.00583 s     | 0.00449 s     | 0.00409 s     | 0.01551 s     | 0.009549       | 725        |
| 64  | 0.03388 s     | 0.02230 s     | 0.02090 s     | 0.10890 s     | 0.003287       | 2781       |
| 128 | 0.4280 s      | 0.2195 s      | 0.1704 s      | 0.2827 s      | 0.001134       | 10435      |
| 256 | 6.4796 s      | 3.4750 s      | 2.0516 s      | 2.3377 s      | 0.000339       | 38595      |

---

## MPI Scaling

MPI parallelism shows clear and consistent speedups, especially at larger grid sizes:

- **N = 128**: going from 1 to 4 MPI processes yields a **~2.5× speedup** (0.428 s → 0.170 s).
- **N = 256**: going from 1 to 4 MPI processes yields a **~3.2× speedup** (6.48 s → 2.05 s), approaching the theoretical maximum of 4× for 4 processes.

At small grids (N = 16), using 4 MPI processes is actually *slower* than the serial case. This is expected an due to the MPI communication overhead.

---

## OpenMP Behaviour

The hybrid configuration consistently *underperforms* the pure-MPI case across all grid sizes. This is limitation is probably due to the hardware configuration exploited to conduct the experiment.

---

The L2 error decreases as N grows (finer grid), while the number of iterations required for convergence grows rapidly with N. This is caused by the slow convergence rate of the Jacobi method.

For optimal performance, the code should be run on a cluster where each MPI rank is pinned to a separate physical node, allowing the OpenMP threads to exploit true shared-memory parallelism without resource contention.