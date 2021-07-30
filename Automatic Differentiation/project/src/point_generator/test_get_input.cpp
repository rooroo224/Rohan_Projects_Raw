#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "get_input.cpp"

int main(int argc, char** argv){

	std::vector<double> left_limits;
	std::vector<double> right_limits;
	std::vector<int> num_samples;
	
	int num_dimensions;
	bool is_grid;

	get_input("./input.in",  num_dimensions, is_grid, left_limits,  right_limits,  num_samples);

	std::cout <<"is grid "<< is_grid << std::endl;
	for(int i= 0; i < num_dimensions ;i++){
			std::cout << left_limits[i] << " "<<right_limits[i] << " "<< num_samples[i] << std::endl;
	}
	
	return 0;
}
