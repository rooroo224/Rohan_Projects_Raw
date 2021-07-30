/** @file sampler_test.cpp
 *  @author alper
 *  @version 0.1
 *  @brief Test for uniform random sampler
 *  @date Jan, 26,2021
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include<vector>
#include <cassert>

void random_sampling(int num_points, double ul, double ur, std::vector<double>& x_sampling);

#define INPUTS 4

int main(int argc, char** argv){
	
	assert(argc == INPUTS && "Wrong number of arguments");
	
	double glob_lower =  std::atof(argv[1]), glob_upper = std::atof(argv[2]);
	int num_points = std::atof(argv[3]) ;
	
	std::vector<double> x_sampling;
	random_sampling( num_points, glob_lower, glob_upper, x_sampling);
	
	std::ofstream myfile;
	myfile.open ("test_sampler.dat");
	for(int i = 0; i < num_points; i++ ) 
		myfile << x_sampling[i] << std::endl;
	myfile.close();

	return 0;
}
