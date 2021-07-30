#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "test_drivers.cpp"
#include "calculator.cpp"

#include "../point_generator/generator.cpp"
#include "../point_generator/get_input.cpp"

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
		
		
  	std::cout << "num_dim:" << num_dimensions << std::endl;  		

   	
	for (int i = 0; i< num_dimensions; i++){
		for (size_t ii = 0; ii< x_sampling[i].size(); ii++)
		  std::cout << x_sampling[i][ii] << " ";
		std::cout <<  std::endl;
	}
	

	std::vector<double>* first_order_directional_derivative_tan;
	std::vector<double>* first_order_directional_derivative_adj;
	std::vector<double>* first_order_directional_derivative_FD;
	double* calculated_values; 
	
	
	calculator( num_dimensions, domain_size, domain_mesh,
			 first_order_directional_derivative_tan,
			 first_order_directional_derivative_adj,
			 first_order_directional_derivative_FD, 
			 calculated_values);
	
	for (int i = 0; i< domain_size; i++)
		  std::cout << calculated_values[i] << " ";
	std::cout <<  std::endl;
	
     	for (int i = 0; i< domain_size; i++){
     		std::cout << i << ": ";
		for (size_t ii = 0; ii< num_dimensions; ii++)
		  std::cout <<  first_order_directional_derivative_tan[i][ii] << " ";
		std::cout <<  std::endl;
	}
	
	return 0;
}
