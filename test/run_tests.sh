#!/bin/bash

echo "Compilation..."
make distclean > /dev/null 2>&1
make > /dev/null 2>&1

echo "==================================================================================="
echo "                        BENCHMARK JACOBI SOLVER"
echo "==================================================================================="

FORMULA="8*pi^2*sin(2*pi*x)*sin(2*pi*y)"
N_VALUES=(16 32 64 128 256)

for N in "${N_VALUES[@]}"; do
    echo ">>> TEST N = $N"
    
    export OMP_NUM_THREADS=1
    echo -n "   [1 MPI, 1 OMP]  "
    mpiexec -n 1 ./main $N "$FORMULA" | grep "Time ="

    export OMP_NUM_THREADS=1
    echo -n "   [2 MPI, 1 OMP]  "
    mpiexec -n 2 ./main $N "$FORMULA" | grep "Time ="

    export OMP_NUM_THREADS=1
    echo -n "   [4 MPI, 1 OMP]  "
    mpiexec -n 4 ./main $N "$FORMULA" | grep "Time ="

    export OMP_NUM_THREADS=2
    echo -n "   [4 MPI, 2 OMP]  "
    mpiexec -n 4 ./main $N "$FORMULA" | grep "Time ="
    
done

echo "Benchmark completed"