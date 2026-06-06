# Reproducing the Benchmark Results

## Prerequisites

Make sure the following are installed and available on your system:

- Eigen3 headers (either via `pacs-examples` or system-wide at `/usr/include/eigen3`)
- `muParser` library

---

## Build

From the **project root**, adjust `PACS_ROOT` in the `Makefile` to point to your local `pacs-examples` directory (or leave it as-is if Eigen is installed system-wide at `/usr/include/eigen3`):

```bash
make
```

This compiles all sources under `src/` with `mpicxx -std=c++20 -O3 -fopenmp` and links against `muparser`, producing the `main` executable.

---

## Running the Scalability Benchmark

The script `test/run_tests.sh` automates the full benchmark. It recompiles the project from scratch, then runs the solver for grid sizes `N ∈ {16, 32, 64, 128, 256}` across four configurations:

| Configuration    | MPI processes | OMP threads |
|-----------------|:-------------:|:-----------:|
| [1 MPI, 1 OMP]  | 1             | 1           |
| [2 MPI, 1 OMP]  | 2             | 1           |
| [4 MPI, 1 OMP]  | 4             | 1           |
| [4 MPI, 2 OMP]  | 4             | 2           |

From the project root, run:

```bash
chmod +x test/run_tests.sh
cd test && ./run_tests.sh
```

The script sets `OMP_NUM_THREADS` before each `mpiexec` invocation, so no manual environment setup is needed.

Each run invokes the solver as:
```bash
mpiexec -n <P> ./main <N> "8*pi^2*sin(2*pi*x)*sin(2*pi*y)"
```

and prints a single line per configuration with wall-clock time, L2 error, and iteration count.

---

## Parameters

The solver is configured in `src/main.cpp` with:
- **Max iterations**: 100,000
- **Tolerance**: 1e-8
- **Boundary conditions**: `g(x,y) = 0` (modifiable in `main.cpp`)
- **Exact solution** (for L2 error): `u_ex(x,y) = sin(2πx)sin(2πy)`

To test with the inhomogeneous BCs, update the `dirichlet_bc` lambda and `u_ex` function in `main.cpp` accordingly, then recompile. For example: 
- `g(x,y) = x^2 + y^2`,
- `u_ex(x,y) = x^2 + y^2`, 
- `f(x,y) = -4`
 
should converge with L2 error near to zero.

---

## Hardware

Benchmarks were run on:

**Intel Core i7-1195G7 @ 2.90 GHz** — 4 physical cores, 8 logical threads (Hyper-Threading), 12 MB L2 cache.  
See `hw.info` in the project root for full details.