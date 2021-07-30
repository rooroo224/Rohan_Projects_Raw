#include <iostream>
#include <fstream>
#include <vector>
//using namespace std;

int post_process_2D(int& num_dimensions, int& domain_size, std::vector<double>*& domain_mesh, std::vector<double>*& deriv_data){

	std::ofstream myfile_out;

	myfile_out.open ("data.dat");
  	if (myfile_out.is_open()){   //checking whether the file is open	   

  		for (int i = 0; i< domain_size; i++){
		      	for (int j = 0; j< num_dimensions; j++){ 
				
				myfile_out <<	domain_mesh[i][j] << " ";
				}
				
		      	for (int j = 0; j< num_dimensions; j++){ 
				
				myfile_out <<	deriv_data[i][j] << " ";
				}
			myfile_out << "\n";
		}
			
		myfile_out.close();   //close the file object.		
   	}
   	else return 1; 
}
