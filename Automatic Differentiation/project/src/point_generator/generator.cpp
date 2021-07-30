/** @file get_input.cpp
 *  @author alper
 *  @version 0.1
 *  @brief Read domain bounds from input file.
 *  @date Feb, 1,2021
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include "../sampler/sampler.cpp"

int generator(int num_dimensions, bool is_grid, std::vector<double>& left_limits,
					  std::vector<double>& right_limits,
					   std::vector<int>& num_samples,
		 int& domain_size,
		 std::vector<double>*& x_sampling, 
		 std::vector<double>*& domain_mesh){
	
	
	
	domain_size = 1 ; 
	for (int i = 0; i< num_dimensions; i++) domain_size *= num_samples[i];
	
	x_sampling = new std::vector<double>[num_dimensions];

	double delta;
	if (is_grid){
	      	for (int i = 0; i< num_dimensions; i++){  //read data from file object and put it into string.
	      	
	      		x_sampling[i].reserve(num_samples[i]);
	      		delta =  (right_limits[i] - left_limits[i]) / (num_samples[i] - 1);
	      		for (int ii = 0 ; ii < num_samples[i]; ii++)
	      			 x_sampling[i].emplace_back(left_limits[i] + ii*delta);
		}
	}
	else
	      	for (int i = 0; i< num_dimensions; i++){  //read data from file object and put it into string.
	      		random_sampling( domain_size,  left_limits[i],  right_limits[i], x_sampling[i], is_grid);
		}


	domain_mesh = new std::vector<double>[domain_size];
	
   	if (is_grid){		
		int index;
		int multiplier ;
		
		for (int i = 0; i< domain_size; i++){
		
			domain_mesh[i].reserve(num_dimensions);
		
			for (size_t ii = 0; ii< num_dimensions; ii++){
				multiplier =1;
				for (size_t k = ii+1; k< num_dimensions; k++)
					multiplier *= x_sampling[k].size();
					
				index = (i/multiplier)% x_sampling[ii].size();
				domain_mesh[i].emplace_back(x_sampling[ii][index]);
			}
		}
	}
	else{
		for (int i = 0; i< domain_size; i++){		
			domain_mesh[i].reserve(num_dimensions);		
			for (size_t ii = 0; ii< num_dimensions; ii++){
				domain_mesh[i].emplace_back(x_sampling[ii][i]);
			}
		}

	}
	return 0;
}
