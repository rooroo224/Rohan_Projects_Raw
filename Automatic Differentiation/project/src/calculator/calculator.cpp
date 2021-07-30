#include<iostream>
#include<vector>
#include<math.h>
#include<cmath>
#include <fstream>


#include "../functions.hpp"
	
void calculator(int& n, int& all_point,  std::vector<double>*& xv,
		std::vector<double>*& yt,
		std::vector<double>*& xa,
		std::vector<double>*& hessian_xa,
		std::vector<double>*& hessian_yt,
		std::vector<double>*& y_fd,
		std::vector<double>*& hessian_diag,
		double*& yv){

	std::vector<double> xt(n,0);  
	
	yt = new std::vector<double>[all_point];
	xa = new std::vector<double>[all_point];
	hessian_xa = new std::vector<double>[all_point];
	hessian_yt = new std::vector<double>[all_point];
	y_fd = new std::vector<double>[all_point];	
	
	std::vector<std::vector<double>> hessian( n , std::vector<double> (n, 0));	
	hessian_diag = new std::vector<double>[all_point];	
	
	yv = new double [all_point];
	auto yv_a = new double [all_point];
	
	int count=0;	
	while(count<all_point){

		yt[count].reserve(n);
		dfdx_tangent(xv[count],xt,yv[count],yt[count]);

		
		xa[count].reserve(n);
		for (int i=0;i<n;i++) xa[count][i]=1; 
		dfdx_adjoint(xv[count],xa[count],yv_a[count]);
		
		y_fd[count].reserve(n);
		dfdx_fd(xv[count],y_fd[count]);
		
		hessian_xa[count].reserve(n);
       		hessian_yt[count].reserve(n);
     		t2s_a1s_Hessian(xv[count], yv[count], hessian_xa[count], hessian_yt[count], hessian);
	// extract diagomal elements
		hessian_diag[count].reserve(n);
		for (int i=0;i<n;i++){
			hessian_diag[count][i] = hessian[i][i];
		}  
		
		
	count++;
	}
	
	delete[] yv_a;
}	

