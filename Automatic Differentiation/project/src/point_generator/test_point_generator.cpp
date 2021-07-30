/** @file get_input.cpp
 *  @author alper
 *  @version 0.1
 *  @brief Read domin bounds from input file.
 *  @date Feb, 1,2021
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "generator.cpp"
#include "get_input.cpp"

int main(int argc, char** argv){


	std::vector<double> left_limits;
	std::vector<double> right_limits;
	std::vector<int> num_samples;
	
	int num_dimensions;

	if( get_input("./src/point_generator/test_input.dat", num_dimensions, left_limits,  right_limits,  num_samples)){
		std::cout << "could not open the file" << std::endl;
		return 1;
	}
		

	std::vector<double>* x_sampling;
	std::vector<double>* domain_mesh;
	int domain_size;
	 	
	
	if( generator(num_dimensions, left_limits,  right_limits,  num_samples,
			 domain_size, x_sampling, domain_mesh)){
		 	
		std::cout << "problem with generator" << std::endl;
		return 1;
	}
		
	
	std::ifstream input_file;
	std::string line, left, right, num_points;
	num_dimensions = int(domain_mesh[0].size());
		
  	std::cout << "num_dim:" << num_dimensions << std::endl;  		

   	
	for (int i = 0; i< num_dimensions; i++){
		for (size_t ii = 0; ii< x_sampling[i].size(); ii++)
		  std::cout << x_sampling[i][ii] << " ";
		std::cout <<  std::endl;
	}
	

	
	for (int i = 0; i< domain_size; i++){
		for (size_t ii = 0; ii< num_dimensions; ii++)
		  std::cout << i << ": "<<  domain_mesh[i][ii] << " ";
		std::cout <<  std::endl;
	}
	
	return 0;
}
