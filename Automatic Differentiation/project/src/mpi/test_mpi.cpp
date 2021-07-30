/** @file main.cpp
 *  @author alper
 *  @version 0.1
 *  @brief MAster C code to run and divide the work between nodes.
 *  @date Jan, 5,2021
 */

#include <mpi.h>
#include <iostream>
#include <fstream>

#include "dco.hpp"

functions.hpp

#define INPUTS 4

int main(int argc, char** argv){

	assert(argc == INPUTS && "Wrong number of arguments");

	MPI_Init(&argc, &argv);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	std::cout << "Rank " << world_rank << " is online." << std::endl; 

	double glob_lower =  std::atof(argv[1]), glob_upper = std::atof(argv[2]);
	int glob_slices = std::atof(argv[3]) -1 ;
	double steps = (glob_upper - glob_lower) / glob_slices;

	int left_over_points = (glob_slices + 1) % world_size;
	int my_slices = int((glob_slices +1) / world_size); 
	 
	double my_lower = glob_lower + steps * my_slices * world_rank;
	double my_upper = my_lower + steps * (my_slices -1);

	if (world_rank == (world_size-1)){
	
		my_slices += (left_over_points); 
		my_upper += steps * (left_over_points);
		
	}

	
//work starts here
	MPI_Barrier (MPI_COMM_WORLD);
	double* data = new double[my_slices];

	fun(my_lower, my_upper, steps, my_slices, data);
	
	MPI_Barrier (MPI_COMM_WORLD);	
	std::cout << "Rank " << world_rank << ": my L, U, S is:" << my_lower<< "/"<< my_upper<< "/" << my_slices <<  std::endl; 

//write to file
	
	MPI_Barrier (MPI_COMM_WORLD);
	int writing_rank = 0;
	std::ofstream myfile;

	while (writing_rank < world_size) {
		if (world_rank== writing_rank) {
			std::cout << "Rank " << world_rank << " is writing." << std::endl; 	
			
			if (world_rank == 0) myfile.open ("data.txt");
			else myfile.open ("data.dat", std::fstream::app);
			
			for (int i = 0; i < my_slices; i++){
				std::cout << "Rank " << world_rank << ": " << data[i]<<" is data." << std::endl; 	
				myfile << data[i] << std::endl;
			}
			myfile.close();						
		}
		writing_rank ++;
		MPI_Barrier (MPI_COMM_WORLD);
	}

	std::cout << "Rank " << world_rank << " finished all." << std::endl; 
	MPI_Finalize();

	return 0;
}
