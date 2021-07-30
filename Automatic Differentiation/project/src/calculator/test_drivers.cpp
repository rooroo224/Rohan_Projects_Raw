

#define _USE_MATH_DEFINES
using namespace std;

#include "dco.hpp"
using namespace dco;


template<typename T>
void f(const vector<T>& x, T& y) {
  const size_t n=x.size();T temp = 0;
  for(size_t i = 0; i<n ; i++)
	temp += x[i];
	y = 1/exp(-temp);
}

template<typename T>
void dfdx_fd(const vector<T>& x, vector<T>& y_fd) {
	const size_t n=x.size();
	const double h = 1e-5;
	vector<double> x_l(n), x_r(n); x_l = x_r = x;
	double y_l, y_r;
	for(size_t i = 0; i < n; i++){
		x_l[i] -= h; x_r[i] += h;
		f(x_r,y_r); f(x_l,y_l);
		y_fd[i] = (y_r-y_l)/(2*h);
		x_l[i] += h; x_r[i] -= h;
		
	}
}


void t2s_a1s_Hessian(const vector<double>& xv,double& yv,vector<vector<double> >& hessian) {
  typedef gt1s<double> DCO_BASE_MODE;
  typedef DCO_BASE_MODE::type DCO_BASE_TYPE;
  typedef ga1s<DCO_BASE_TYPE> DCO_MODE;
  typedef DCO_MODE::type DCO_TYPE;
  typedef DCO_MODE::tape_t DCO_TAPE_TYPE;

  DCO_MODE::global_tape=DCO_TAPE_TYPE::create();

  const size_t n=xv.size();    // the length of the vector is returned with .size() method with returntype size_t

  vector<DCO_TYPE> x(n); DCO_TYPE y;   //declaring input vector and output scalar 

  for (size_t i=0;i<n;i++) 
	x[i]=xv[i];

  for (size_t i=0;i<n;i++) {

    for (size_t j=0;j<n;j++)
      DCO_MODE::global_tape->register_variable(x[j]);       // registering all the DCO input variables into the tape

      derivative(value(x[i])) = 1.0;    // seeding 

      f(x,y);     // y is modified by reference

      DCO_MODE::global_tape->register_output_variable(y);       // registering o/p variable into the tape
      derivative(y) = 1.0;

      DCO_MODE::global_tape->interpret_adjoint();

    for (size_t j = 0; j <= i; ++j)
      hessian[i][j] = derivative( derivative(x[j]) );    // this is desired second order derivative

      derivative(value(x[i])) = 0.0;
      DCO_MODE::global_tape->reset();
  }

  DCO_TAPE_TYPE::remove(DCO_MODE::global_tape);
  yv = passive_value(y);
}


void dfdx_tangent(const vector<double>& xv, const vector<double>& xt, double& yv, vector<double>& yt) {

	typedef double DCO_BASE_TYPE;
	typedef gt1s<DCO_BASE_TYPE> DCO_MODE; 
	typedef DCO_MODE::type DCO_TYPE;
	const size_t n=xv.size();

	vector<DCO_TYPE> x(n);  DCO_TYPE y;
	value(x)=xv;
	derivative(x)=xt;
	for (size_t i=0;i<n;i++) {
		value(x)=xv;
		derivative(x[i])=1;   // seed
		f(x,y);
		yt[i] = derivative(y); // harvest
		derivative(x[i])=0;   // reset seed	
	}	
	yv=value(y);
}

void dfdx_adjoint(const vector<double>& xv, vector<double>& xa, double& yv, double& ya) {

	typedef double DCO_BASE_TYPE;
	typedef ga1s<DCO_BASE_TYPE> DCO_MODE;
	typedef DCO_MODE::type DCO_TYPE;
	typedef DCO_MODE::tape_t DCO_TAPE_TYPE;
	size_t n=xv.size();

	DCO_MODE::global_tape=DCO_TAPE_TYPE::create();
	
	vector<DCO_TYPE> x(n); DCO_TYPE y;

	for (size_t i=0;i<n;i++) {
		x[i]=xv[i];
		DCO_MODE::global_tape->register_variable(x[i]);
	}
	f(x,y);
	yv=value(y); 
	derivative(y)=ya;
	DCO_MODE::global_tape->interpret_adjoint();
	for (size_t i=0;i<n;i++) xa[i]=derivative(x[i]);
	ya=derivative(y);
	DCO_TAPE_TYPE::remove(DCO_MODE::global_tape);
}

