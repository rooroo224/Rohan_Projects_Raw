/** @file check_nodes.cpp
 *  @author alper
 *  @version 0.1
 *  @brief C/mpi code to ping the nodes.
 *  @date Jan, 5,2021
 */

#include <mpi.h>
#include <iostream>

int main(){

	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	std::cout << "Rank " << world_rank << " is online." << std::endl; 

	MPI_Finalize();

	return 0;
}
