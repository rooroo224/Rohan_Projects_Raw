/** @file fun.cpp
 *  @author alper
 *  @version 0.1
 *  @brief MAster C code to run and divide the work between nodes.
 *  @date Jan, 10, 2021
 */

#include <iostream>

void fun(double lower, double upper, double steps, int slices, double* data){


	for (int i = 0; i < slices; i++){
		double x =  (lower + steps * i);
		data[i] = x*x;
	}
}

