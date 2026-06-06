#include "JacobiSolver.hpp"

#include <mpi.h> 
#include <omp.h>

size_t JacobiSolver::solve(Domain &d) const
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rows = d.n;
    int cols = d.n;
    // internal nodes are partitioned among MPI processes 
    int local_rows = (rows-2) / size;
    // if rows and processes are not multiples, compute the remaining rows
    int rest = (rows-2) % size;
    
    // global offset compuation for local indexes mapping into the whole grid
    int global_offset = (rank < rest) ? rank*(local_rows+1) : 
                                        (rank-rest)*local_rows + rest*(local_rows+1);
    if (rank < rest) ++local_rows;
    

    Matrix local_U = Matrix::Zero(local_rows+2, cols);

    Real h = d.h;
    Function f = d.problem.force;

    Matrix local_F = Matrix::Zero(local_rows+2, cols);
    for (size_t i = 1; i < local_rows+1; ++i){
        for (size_t j = 1; j < cols-1; ++j){
            
            Real x = d.getCoord(global_offset + i);
            Real y = d.getCoord(j);
            
            local_F(i,j) = f(x,y); 
        }
    }

    Function g_d = d.problem.dirichlet_bc;

    if (g_d){
        for (size_t i = 0; i < local_rows + 2; ++i) {
            Real x = d.getCoord(global_offset + i);
            local_U(i,0) = g_d(x, d.problem.lb); // y = 0
            local_U(i,cols-1) = g_d(x,d.problem.ub); // y = 1
        }
        if(rank == 0)
            for (size_t j = 0; j < cols; ++j)
                local_U(0,j) = g_d(d.problem.lb, d.getCoord(j)); // x = 0
        
        if(rank == size-1)
            for (size_t j = 0; j < cols; ++j)
                local_U(local_rows + 1, j) = g_d(d.problem.ub, d.getCoord(j)); // x = 1
    }

    Matrix local_U_prev = local_U;

    size_t it = 0;
    int global_converged = 0;
    Real err = tolerance + 100.;

    while (it <= max_iter && global_converged == 0){
        Real sq_diff = 0.;

        #pragma omp parallel for reduction(+:sq_diff)
        for (size_t i = 1; i < local_rows+1; ++i){
            for (size_t j = 1; j < cols-1; ++j){

                Real val = 0.25 * (local_U_prev(i-1,j) + local_U_prev(i+1,j) + 
                                local_U_prev(i,j-1) + local_U_prev(i,j+1) + local_F(i,j)*h*h);

                Real diff = val - local_U_prev(i,j);
                local_U(i,j) = val;
                sq_diff += diff * diff;
            }
        }

        if (rank != 0){
            int top_neigh = rank - 1;
            // Send actual rank's top row and Receive top_neigh's bottom row
            auto top_row = local_U.row(1);
            auto top_ghost = local_U.row(0);
            MPI_Sendrecv(top_row.data(), cols, MPI_DOUBLE, top_neigh, 0,
                         top_ghost.data(), cols, MPI_DOUBLE, top_neigh, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
        }

        if (rank != size-1){
            int bottom_neigh = rank + 1;
            // Send actual rank's bottom row and Receive bottom_neigh's top row
            auto bottom_row = local_U.row(local_rows);
            auto bottom_ghost = local_U.row(local_rows+1);
            MPI_Sendrecv(bottom_row.data(), cols, MPI_DOUBLE, bottom_neigh, 0,
                         bottom_ghost.data(), cols, MPI_DOUBLE, bottom_neigh, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        // computationally more efficient than assigning every loop local_U to local_U_prev
        local_U_prev.swap(local_U);
        
        err = sqrt(h * sq_diff);
        int local_converged = (err < tolerance) ? 1 : 0;
        MPI_Allreduce(&local_converged, &global_converged, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

        ++it;
    }

    // Gathering of local solutions
    int local_elements = local_rows * cols;

    std::vector<int> recv_counts(size);
    std::vector<int> displacements(size);
    
    if (rank == 0){
        int curr_displ = 0;

        for (int i = 0; i < size; ++i){
            int rank_rows = (rows-2) / size;
            if (i < rest) ++rank_rows;

            recv_counts[i] = rank_rows * cols;
            displacements[i] = curr_displ;

            curr_displ += recv_counts[i];
        }
    }

    MPI_Gatherv(local_U_prev.row(1).data(), local_elements, MPI_DOUBLE,
                d.U.row(1).data(), recv_counts.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0){
        for(size_t k = 0; k < d.n; ++k){
            Real coord = d.getCoord(k);
            d.U(0,k) = g_d(d.problem.lb, coord); // Top
            d.U(d.n-1, k) = g_d(d.problem.ub, coord); // Bottom
            d.U(k, 0) = g_d(coord, d.problem.lb); // Left
            d.U(k, d.n-1) = g_d(coord, d.problem.ub); // Right        
        }
    }

    return it;
}
