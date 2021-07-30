#include<iostream>
#include<vector>
#include<math.h>
#include<cmath>

#include <string>
#include <fstream>

#define _USE_MATH_DEFINES
using namespace std;

#include "dco.hpp"
using namespace dco;


vector<vector<double>> first_tangent(const int n,const int m, double ul, double ur, int num_points)
{
	double yv;
	vector<double> xv(n);
	vector<double> xt(n);
	vector<double> yt(n);
	vector<vector<double>> hessian(n, vector<double>(n));
	double c;
	
	vector<double> x_sampling;
	random_sampling(num_points, ul, ur, x_sampling);

	vector<vector<double>> firstordertangent;

	for(int i=0;i<n;i++){xt[i]=0;} //initialize xt

	//go through all points
	int count=0; vector<int> index(n);
	while(count<pow(num_points,n)){
		for(int i=0;i<n;i++){
			index[i]=fmod(count/pow(num_points,i),num_points);
			xv[i]=x_sampling[index[i]];
		}

		dfdx_tangent(xv,xt,yv,yt);
		
		firstordertangent.push_back(yt); //INSTED OF FILE IT IS STORED IN A ARRAY
	
  		count++;
	}
return firstordertangent;
}



vector<vector<double>> second(const int n,const int m, double ul, double ur, int num_points)
{
	double yv;
	vector<double> xv(n);
	vector<double> xt(n);
	vector<double> yt(n);
	vector<vector<double>> hessian(n, vector<double>(n));
	double c;
	vector<double> x_sampling;
	random_sampling(num_points, ul, ur, x_sampling);

	vector<vector<double>> secondorder;

	for(int i=0;i<n;i++){xt[i]=0;} //initialize xt

	//go through all points
	int count=0; vector<int> index(n);
	while(count<pow(num_points,n)){
		for(int i=0;i<n;i++){
			index[i]=fmod(count/pow(num_points,i),num_points);
			xv[i]=x_sampling[index[i]];
		}
		t2s_a1s_Hessian(xv,yv,hessian); //INSTEAD OF FILE IT IS STORED IN A ARRAY
		
		vector<double> temp;	
	
		for(int p = 0; p<hessian.size(); p++){
			for(int q=0; q<hessian[p].size(); q++){
				temp.push_back(hessian[p][q]);
				}
			}
		secondorder.push_back(temp);
		temp.clear();
  		count++;
	}
return secondorder;
}

