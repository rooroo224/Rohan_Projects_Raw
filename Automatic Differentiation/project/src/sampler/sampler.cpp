#include<vector>
#include <random>
#include <algorithm>

void random_sampling(int num_points, double ul, double ur, std::vector<double>& x_sampling, bool is_grid){

	x_sampling.reserve(num_points);
        
	std::random_device rd;
	std::mt19937 gen(rd()); 
	std::uniform_real_distribution<> dis(ul,ur);

	for (int j = 0; j < num_points; j++){
		x_sampling.emplace_back(dis(gen));
	}
	if (is_grid)
		std::sort(x_sampling.begin(), x_sampling.end());
}


