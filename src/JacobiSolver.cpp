#include "JacobiSolver.hpp"

#include <mpi.h> 
#include <omp.h>

void JacobiSolver::solve(Domain &d) const
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rows = d.n;
    int cols = d.n;
    int local_rows = (rows-2) / size;
    int rest = (rows-2) % size;
    
    int global_offset = (rank < rest) ? rank*(local_rows+1) : 
                                        (rank-rest)*local_rows + rest*(local_rows+1);
    if (rank < rest) ++local_rows;
    

    Matrix local_U = Matrix::Zero(local_rows+2, cols);
    Matrix local_U_prev = local_U;

    Real h = d.h;
    Function f = d.problem.force;

    size_t it = 0;
    int global_converged = 0;
    Real err = tolerance + 1.;

    while (it <= max_iter && global_converged == 0){
        Real sq_diff = 0.;
        for (size_t i = 1; i < local_rows+1; ++i){
            for (size_t j = 1; j < cols-1; ++j){
    
                Real x = d.getCoord(global_offset + i);
                Real y = d.getCoord(j);

                local_U(i,j) = 0.25 * (local_U_prev(i-1,j) + local_U_prev(i+1,j) + 
                                local_U_prev(i,j-1) + local_U_prev(i,j+1) + f(x,y)*h*h);

                sq_diff += (local_U(i,j) - local_U_prev(i,j)) * (local_U(i,j) - local_U_prev(i,j));
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

    MPI_Gatherv(local_U.row(1).data(), local_elements, MPI_DOUBLE,
                d.U.row(1).data(), recv_counts.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Se BCs sono = 0 va bene così, ma nisogna aggiungere riga 0 e n-1 se le considero Non Nulle

    std::cout << "Problem Solved" << std::endl;
}
