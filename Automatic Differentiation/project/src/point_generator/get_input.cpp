/** @file get_input.cpp
 *  @author alper
 *  @version 0.1
 *  @brief Read domin bounds from input file.
 *  @date Feb, 1,2021
 */

// #include <iostream>
#include <fstream>
#include <string>
#include <vector>



int get_input(std::string file_name,  int& num_dimensions, bool& is_grid, std::vector<double>& left_limits, std::vector<double>& right_limits, std::vector<int>& num_samples){

	std::ifstream input_file;
	std::string line, left, right, num_points;
		
	
	input_file.open(file_name);
  	if (input_file.is_open()){   //checking whether the file is open	   
  		input_file >> line;
  		num_dimensions = stoi(line);
  		input_file >> line;
  		is_grid = stoi(line);
  		is_grid = not is_grid;
  		
  		left_limits.reserve(num_dimensions);
  		right_limits.reserve(num_dimensions);
  		num_samples.reserve(num_dimensions);
  		
	      	for (int i = 0; i< num_dimensions; i++){  //read data from file object and put it into string.
	      		input_file >> left >> right >> num_points;
	      		
	      		left_limits[i] = std::stod(left);
	      		right_limits[i] = std::stod(right);
	      		num_samples[i] = std::stoi(num_points);
			}
		input_file.close();   //close the file object.		
   	}
   	else return 1; 
   	
	return 0;
}
