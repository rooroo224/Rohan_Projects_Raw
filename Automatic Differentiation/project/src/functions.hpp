// target function
template<typename T>
void f(const std::vector<T>& x, T& y);

// Test function for mpi
void fun(double lower, double upper, double steps, int slices, double* data);

// !point_generator/

// /get_input.cpp
// read input data from text file
int get_input(std::string file_name, int& num_dimensions, bool& is_grid, std::vector<double>& left_limits, std::vector<double>& right_limits, std::vector<int>& num_samples);

// /generator.cpp
// sample from domain and transform n dimensional coordinates to 2D data
int generator(int num_dimensions, bool is_grid, std::vector<double>& left_limits,
				  std::vector<double>& right_limits,
				   std::vector<int>& num_samples,
		 int& domain_size,
		 std::vector<double>*& x_sampling, 
		 std::vector<double>*& domain_mesh);

// !calculator/

// calculate derivatives/calculator.cpp
void calculator(int& n, int& all_point,  std::vector<double>*& xv,
		std::vector<double>*& yt,
		std::vector<double>*& xa,
		std::vector<double>*& hessian_xa,
		std::vector<double>*& hessian_yt,
		std::vector<double>*& y_fd,
		std::vector<double>*& hessian_diag,
		double*& yv);

		
// !driver/drivers.cpp

// dcoc drivers
void dfdx_tangent(const std::vector<double>& xv, const std::vector<double>& xt, double& yv, std::vector<double>& yt);
void dfdx_adjoint(const std::vector<double>& xv, std::vector<double>& xa, double& yv);
void t2s_a1s_Hessian(const std::vector<double>& xv,double& yv, std::vector<double>& xa ,std::vector<double>& yt,std::vector<std::vector<double> >& hessian);
void dfdx_fd(const std::vector<double>& x, std::vector<double>& y_fd);




// !post_process/

// post_process_2D.cpp
int post_process_2D(int& num_dimensions, int& domain_size, std::vector<double>*& domain_mesh, std::vector<double>*& deriv_data);

// report_gen.cpp
int report_gen(int num_dim, std::string data_file);
