/** @file main.cpp
 *  @author alper
 *  @version 0.1
 *  @brief MAster C code to run and divide the work between nodes.
 *  @date Jan, 5,2021
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string>
#include <unistd.h>//sleep


#include "functions.hpp"

#define INPUTS 2

int main(int argc, char** argv){

// check for input params
	assert(argc == INPUTS && "Wrong number of arguments");


// init MPI
	MPI_Init(&argc, &argv);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	
	usleep(world_rank*10000 +10);
	std::cout << "Rank " << world_rank << " is online." << std::endl; 
	
// read form file
	
	std::vector<double> left_limits;
	std::vector<double> right_limits;
	std::vector<int> num_samples;
	
	int num_dimensions;
	
	bool is_grid;
	if( get_input(argv[1], num_dimensions, is_grid, left_limits,  right_limits,  num_samples)){
		std::cout << "could not open the file" << std::endl;
		return 1;
	}
	
	if (is_grid){
		std::cout << "Grid mode enabled" << std::endl;
	}
	else{
		std::cout << "Grid mode disabled: Appending to previous data" << std::endl;
	}
		
//set dimention to be paralelized
	const int paralel_dim = 0;
		
// do work sharing
	double glob_lower =  left_limits[paralel_dim], glob_upper = right_limits[paralel_dim];
	int glob_points = num_samples[paralel_dim] ;
	double steps = (glob_upper - glob_lower) / glob_points;
	int left_over_points = glob_points % world_size;
	
// assign to local nodes
	num_samples[paralel_dim] = int(glob_points / world_size); 	 
	left_limits[paralel_dim] = glob_lower + steps * num_samples[paralel_dim] * world_rank;
	right_limits[paralel_dim] = left_limits[paralel_dim] + steps * (num_samples[paralel_dim]);

	if (world_rank == (world_size-1)){
	
		num_samples[paralel_dim] += (left_over_points); 
		right_limits[paralel_dim] += steps * (left_over_points);		
		
	}
	
// cout work sharing
	MPI_Barrier (MPI_COMM_WORLD);
	int writing_rank = 0;
	while (writing_rank < world_size) {
		if (world_rank== writing_rank) {
			
			std::cout << "rank: " << world_rank
				  << " num_samples in shared dim: " << num_samples[paralel_dim]  
				  << " left_limits: " << left_limits[paralel_dim] 
				  << " right_limits: " << right_limits[paralel_dim] 
				  << std::endl ;
		}
		writing_rank ++;
		MPI_Barrier (MPI_COMM_WORLD);
	}
	

	
//work starts here
	MPI_Barrier (MPI_COMM_WORLD);
//create data structures to pass in	
	std::vector<double>* x_sampling;
	std::vector<double>* domain_mesh;
	
	int domain_size;
// generate domain as a 2D list of coordinates	 	

	if( generator(num_dimensions, is_grid, left_limits,  right_limits,  num_samples,
			 domain_size, x_sampling, domain_mesh)){
		 	
		std::cout << "problem with generator" << std::endl;
		return 1;
	}
// cout work sharing
	MPI_Barrier (MPI_COMM_WORLD);
	writing_rank = 0;
	while (writing_rank < world_size) {
		if (world_rank== writing_rank) 			
			std::cout << "rank: " << world_rank << std::endl << " domain_size: " <<domain_size<< std::endl;

		writing_rank ++;
		MPI_Barrier (MPI_COMM_WORLD);
	}
	
	

/*
// check genarated mesh

	MPI_Barrier (MPI_COMM_WORLD);	
	writing_rank = 0;
	while (writing_rank < world_size) {
		if (world_rank== writing_rank) {
		std::cout << std::endl << "=======" << std::endl;
		std::cout << "rank: " << world_rank << std::endl;	
		
			for (int i = 0; i< num_dimensions; i++){
				std::cout << "--- x_" << i << " ---" << std::endl;
				for (size_t ii = 0; ii< x_sampling[i].size(); ii++)
					std::cout << x_sampling[i][ii] << " ";
				std::cout <<  std::endl;
			}
		}
		writing_rank ++;
		MPI_Barrier (MPI_COMM_WORLD);
	}
	
	
	MPI_Barrier (MPI_COMM_WORLD);	
	writing_rank = 0;
	while (writing_rank < world_size) {
		if (world_rank== writing_rank) {
		std::cout << std::endl << "=======" << std::endl;
		std::cout << "rank: " << world_rank << std::endl;	
		
			for (int i = 0; i< domain_size; i++){
				std::cout << "point_" << i << ": " ;
				for (size_t ii = 0; ii< num_dimensions; ii++)
					std::cout << domain_mesh[i][ii] << " ";
				std::cout <<  std::endl;
			}
		}
		writing_rank ++;
		MPI_Barrier (MPI_COMM_WORLD);
	}
*/

//create data structures to pass in
	std::vector<double>* first_order_directional_derivative_tan;
	std::vector<double>* first_order_directional_derivative_adj;
	std::vector<double>* first_order_directional_derivative_xa;
	std::vector<double>* first_order_directional_derivative_yt;
	std::vector<double>* first_order_directional_derivative_FD;	
	std::vector<double>* second_order_directional_derivative_diagonals_adj;
	double* calculated_values; 
	
//derivatives are calculated inside	
	calculator( num_dimensions, domain_size, domain_mesh,
			 first_order_directional_derivative_tan,
			 first_order_directional_derivative_adj,
			 first_order_directional_derivative_xa,
			 first_order_directional_derivative_yt,
			 first_order_directional_derivative_FD,
			 second_order_directional_derivative_diagonals_adj, 
			 calculated_values);

// check derivatives
/*
	MPI_Barrier (MPI_COMM_WORLD);	 
	writing_rank = 0;
	while (writing_rank < world_size) {
		if (world_rank== writing_rank) {
		std::cout << std::endl << "=======" << std::endl;
		std::cout << "rank: " << world_rank << std::endl;	
		std::cout << "--- calculated values ---" << std::endl;
			for (int i = 0; i< domain_size; i++)
	  			std::cout << calculated_values[i] << " ";
			std::cout <<  std::endl;
			
		std::cout << "\n--- derivatives: tan ---" << std::endl;	
			for (int i = 0; i< domain_size; i++){
		     		std::cout << i << ": ";
				for (size_t ii = 0; ii< num_dimensions; ii++)
				  	std::cout <<  first_order_directional_derivative_tan[i][ii] << " ";
				std::cout <<  std::endl;
			}
			
		std::cout << "\n--- derivatives: adj ---" << std::endl;	
			for (int i = 0; i< domain_size; i++){
		     		std::cout << i << ": ";
				for (size_t ii = 0; ii< num_dimensions; ii++)
				  	std::cout <<  first_order_directional_derivative_adj[i][ii] << " ";
				std::cout <<  std::endl;
			}
						
		std::cout << "\n--- second - derivatives ---" << std::endl;	
			for (int i = 0; i< domain_size; i++){
		     		std::cout << i << ": ";
				for (size_t ii = 0; ii< num_dimensions; ii++)
				  	std::cout <<  second_order_directional_derivative_diagonals_adj[i][ii] << " ";
				std::cout <<  std::endl;
			}						
		std::cout << "\n--- derivatives: FD ---" << std::endl;	
			for (int i = 0; i< domain_size; i++){
		     		std::cout << i << ": ";
				for (size_t ii = 0; ii< num_dimensions; ii++)
				  	std::cout <<  first_order_directional_derivative_FD[i][ii] << " ";
				std::cout <<  std::endl;
			}
			
		}
		writing_rank ++;
		MPI_Barrier (MPI_COMM_WORLD);
	}
*/
// post process for gnu plot

	auto plot_data = new std::vector<double>[domain_size];
	
	double abs_first_der, abs_second_der;
	
	for (int i = 0; i< domain_size; i++){
	
		plot_data[i].reserve(num_dimensions*2);
		
		for (int j = 0; j< num_dimensions; j++){ 
				
			plot_data[i][j] = domain_mesh[i][j];
		}

		for (int j = 0; j< num_dimensions; j++){ 
			
			abs_first_der = std::abs(first_order_directional_derivative_tan[i][j]) ;
			//abs_first_der = abs_first_der * abs_first_der ;
			
			if(std::abs(second_order_directional_derivative_diagonals_adj[i][j]) < 0.05 )
				abs_second_der = 0.05;

			plot_data[i][j+num_dimensions] =   abs_first_der * abs_second_der ; //constant
		}
			
	}
//write to files
	
	MPI_Barrier (MPI_COMM_WORLD);
	
	struct Data_file{	
		std::ofstream file_obj;
		std::string name;
		std::vector<double>* data;
		double* data_1D;
		int num_rows;
		int num_cols;
	};
	const int num_files = 7;	
	Data_file data_file[num_files];
	
	data_file[0].name = "Mesh_points.dat";
	data_file[0].data = domain_mesh;
	data_file[0].num_rows = domain_size;
	data_file[0].num_cols = num_dimensions;
	
	data_file[1].name = "First_order_tangent_values.dat";
	data_file[1].data = first_order_directional_derivative_tan;
	data_file[1].num_rows = domain_size;
	data_file[1].num_cols = num_dimensions;
	
	data_file[2].name = "First_order_adjoint_values.dat";
	data_file[2].data = first_order_directional_derivative_adj;
	data_file[2].num_rows = domain_size;
	data_file[2].num_cols = num_dimensions;
	
	data_file[3].name = "First_order_FD_values.dat";
	data_file[3].data = first_order_directional_derivative_FD;
	data_file[3].num_rows = domain_size;
	data_file[3].num_cols = num_dimensions;
	
	data_file[4].name = "Second_order_adjoint_values.dat";
	data_file[4].data = second_order_directional_derivative_diagonals_adj;
	data_file[4].num_rows = domain_size;
	data_file[4].num_cols = num_dimensions;
	
	data_file[5].name = "Calculated_values.dat";
	data_file[5].data_1D = calculated_values;
	data_file[5].num_rows = domain_size;
	data_file[5].num_cols = 1;
	
	data_file[6].name = "First_order_tangent_plot.dat";
	data_file[6].data = plot_data;
	data_file[6].num_rows = domain_size;
	data_file[6].num_cols = num_dimensions*2;
	
	
	int writing_epocs = num_files + world_size -1;
	int epoc_counter = 0;

	while (epoc_counter < writing_epocs) {
		for (int k = 0; k < num_files; k++){		
			if (world_rank == epoc_counter - k) {
				usleep(world_rank*100000 +10);
				std::cout << "epoc: "<< epoc_counter << ", Rank " << world_rank << " is writing to file: " << k << ", " <<  data_file[k].name << std::endl; 	
				
				if ((world_rank == 0 ) &&  (is_grid)) data_file[k].file_obj.open (data_file[k].name);
				else data_file[k].file_obj.open (data_file[k].name, std::fstream::app);
				
				if (data_file[k].num_cols == 1)
					for(int i = 0; i<data_file[k].num_rows; i++){	
						data_file[k].file_obj << data_file[k].data_1D[i] << std::endl;
					}
				else	
					for(int i = 0; i<data_file[k].num_rows; i++){
						for(int j=0; j<data_file[k].num_cols;j++) 	
							data_file[k].file_obj << data_file[k].data[i][j]<<" ";
							
						data_file[k].file_obj << std::endl; 
					}
				data_file[k].file_obj.close();
				break;						
			}
		}
		epoc_counter ++;
		MPI_Barrier (MPI_COMM_WORLD);
	}

	std::cout << "Rank " << world_rank << " finished all." << std::endl; 	
	
	if (world_rank == 0){
		if (num_dimensions == 1)	
		std::cout << "Use:\n	python3 src/python/prepare_1D_pdf.py && pdflatex report.tex >/dev/null  && rm *.log *.out *.aux" << std::endl; 
		if (num_dimensions == 2)	
		std::cout << "Use:\n	python3 src/python/prepare_2D_pdf.py && pdflatex report.tex >/dev/null  && rm *.log *.out *.aux" << std::endl; 		
	}

	MPI_Finalize();
	
	return 0;
}
