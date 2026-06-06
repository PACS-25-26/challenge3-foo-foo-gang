# Parallel Jacobi Solver for the Laplace Equation

## Overview

This project implements a **parallel matrix-free solver** for the 2D Laplace equation on a unit square domain:

```
−∆u = f(x,y),   (x,y) ∈ (0,1)²
u = g(x,y),      on ∂Ω  (Dirichlet BCs)
```

The solver is based on the **Jacobi iterative method** and exploits **hybrid MPI + OpenMP parallelism**: the computational grid is partitioned into row-blocks distributed across MPI processes, while each process further parallelizes its local loops with OpenMP threads.

The problem is solved with forcing term `f(x,y) = 8π²sin(2πx)sin(2πy)`, whose exact solution is `u(x,y) = sin(2πx)sin(2πy)`. Both the forcing term and the Dirichlet boundary conditions can be specified as user-defined functions (via `muParser` expressions from the command line for `f`, and as lambdas in `main.cpp` for the BCs).

---

## Repository Structure

```
.
├── include/
│   ├── TypeTraits.hpp          # Type aliases
│   ├── LaplaceProblem.hpp      # LaplaceProblem struct and Domain class
│   ├── JacobiSolver.hpp        # JacobiSolver class
│   └── muParserInterface.hpp   # muParser wrapper for parsing f(x,y)
├── src/
│   ├── LaplaceProblem.cpp      # Domain methods: L2 error, VTK export
│   ├── JacobiSolver.cpp        # Core parallel Jacobi algorithm
│   └── main.cpp                # Main file
├── test/
│   ├── run_tests.sh            # Scalability benchmark script
|   ├── hw.info                 # Hardware info
│   └── data/                   # Collected benchmark results
└── Makefile
```

---

## Key Design Choices

**Row-block decomposition** — the `n-2` internal rows are split as evenly as possible among MPI ranks (remainder rows go to the first ranks). Each process holds a local submatrix with two extra *ghost rows* to store values received from neighboring ranks.

**Ghost-row exchange** — at each iteration, adjacent ranks exchange boundary rows via `MPI_Sendrecv`.

**Convergence check** — each rank computes its local squared-difference norm; `MPI_Allreduce` with `MPI_PROD` ensures global convergence only when *all* ranks satisfy the tolerance criterion.

**OpenMP parallelism** — the inner double loop over local rows is parallelized with `#pragma omp parallel for reduction(+:sq_diff)`.

**Solution gathering** — `MPI_Gatherv` collects all local solutions into rank 0's global matrix `d.U`, which is then exported to VTK format for visualization in ParaView.

---

The Makefile expects Eigen to be found either under `$(PACS_ROOT)/eigen3` or `/usr/include/eigen3`. Adjust `PACS_ROOT` at the top of the Makefile if needed.