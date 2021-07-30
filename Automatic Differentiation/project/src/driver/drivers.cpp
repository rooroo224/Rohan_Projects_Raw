#define _USE_MATH_DEFINES
using namespace std;

#include "dco.hpp"
using namespace dco;

#include "../../f.cpp"

void dfdx_tangent(const vector<double>& xv, const vector<double>& xt, double& yv, vector<double>& yt) {

    typedef double DCO_BASE_TYPE;
    typedef gt1s<DCO_BASE_TYPE> DCO_MODE; // tangent mode
    typedef DCO_MODE::type DCO_TYPE; // tangent type
    const size_t n=xv.size(); // size of gradient

    vector<DCO_TYPE> x(n);  DCO_TYPE y; // activate
    value(x)=xv;
    derivative(x)=xt;
    for (size_t i=0;i<n;i++) {
        derivative(x[i])=1;   // seed x^(1)=1
        f(x,y);
        yt[i] = derivative(y); // harvest
        derivative(x[i])=0;   // reset seed x^(1)=0
    }    
    yv=value(y);
}

void dfdx_adjoint(const vector<double>& xv, vector<double>& xa, double& yv) { 

    typedef double DCO_BASE_TYPE;
    typedef ga1s<DCO_BASE_TYPE> DCO_MODE; // adjoint mode
    typedef DCO_MODE::type DCO_TYPE; // adjoint type
    typedef DCO_MODE::tape_t DCO_TAPE_TYPE; // tape type
    size_t n=xv.size(); // size of gradient

    DCO_MODE::global_tape=DCO_TAPE_TYPE::create(); // touch type
    
    vector<DCO_TYPE> x(n); DCO_TYPE y; //activate

    for (size_t i=0;i<n;i++) {
        x[i]=xv[i];
        DCO_MODE::global_tape->register_variable(x[i]); // record active input
    }
    f(x,y);
    yv=value(y); // harvest
    derivative(y)=1.0; //seed y_(1)=1
    DCO_MODE::global_tape->interpret_adjoint(); // interpret tape
    for (size_t i=0;i<n;i++) xa[i]=derivative(x[i]); // harvest
    DCO_TAPE_TYPE::remove(DCO_MODE::global_tape); // deallocates global tape 
}

void t2s_a1s_Hessian(const vector<double>& xv,double& yv, vector<double>& xa ,vector<double>& yt,vector<vector<double> >& hessian) {
    typedef gt1s<double> DCO_BASE_MODE; //tangent (based) mode
    typedef DCO_BASE_MODE::type DCO_BASE_TYPE; //tangent type
    typedef ga1s<DCO_BASE_TYPE> DCO_MODE; //base type of adjoint mode
    typedef DCO_MODE::type DCO_TYPE; //adjoint type
    typedef DCO_MODE::tape_t DCO_TAPE_TYPE;// tape type
    DCO_MODE::global_tape=DCO_TAPE_TYPE::create(); // touch tape
    const size_t n=xv.size(); //size of gradient

    vector<DCO_TYPE> x(n); DCO_TYPE y;   //declaring input vector and output scalar 
    for (size_t i=0;i<n;i++) x[i]=xv[i]; //set
    for (size_t i=0;i<n;i++) {
        for (size_t j=0;j<n;j++) DCO_MODE::global_tape->register_variable(x[j]); //registering all the DCO input variables into the tape
        derivative(value(x[i])) = 1.0; // seed tangent x^(2)=1
        f(x,y);     
        DCO_MODE::global_tape->register_output_variable(y);// registering o/p variable into the tape
        derivative(y) = 1.0; //seed adjoint y_(1)=1
        DCO_MODE::global_tape->interpret_adjoint(); //interpret
        for (size_t j = 0; j <= i; ++j) hessian[i][j] = derivative( derivative(x[j]) );    //harvest 2nd-order derivative
        xa[i] = value(derivative(x[i]));// harvest gradient (adjoint mode)
        yt[i] = derivative(value(y)); // harvest grdient (tangent mode) 
        derivative(value(x[i])) = 0.0; //reset
        DCO_MODE::global_tape->reset();
    }
    DCO_TAPE_TYPE::remove(DCO_MODE::global_tape);//deallocate tape
    yv = passive_value(y);
}



void dfdx_fd(const vector<double>& x, vector<double>& y_fd) {
    const size_t n=x.size(); //size of gradient
    const double h = 1e-5; //perturbation
    vector<double> x_l(n), x_r(n); x_l = x_r = x;
    double y_l, y_r;
    for(size_t i = 0; i < n; i++){
        x_l[i] -= h; x_r[i] += h; 
        f(x_r,y_r); f(x_l,y_l); //perturb function to left and right
        y_fd[i] = (y_r-y_l)/(2*h); //finite difference quotient
        x_l[i] += h; x_r[i] -= h; //remove perturbation
        
    }
}

template<typename T>
void ddfdx_fd(const vector<T>& x, vector<vector<T> >& sec_fd) {
    const int n=x.size();//size of gradient
    const double h = 1e-5;//perturbation
    vector<double> x_l(n), x_r(n); x_l = x_r = x;
    vector<double> dfdx_l(n), dfdx_r(n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j <= i; j++){
            x_l[j] -= h; x_r[j] += h;
            dfdx_fd(x_r,dfdx_r); //perturb function to right 
            dfdx_fd(x_l,dfdx_l); //perturb function to left
            sec_fd[i][j] = (dfdx_r[i]-dfdx_l[i])/(2*h);//finite difference quotient
            x_l[i] += h; x_r[j] -= h; //remove perturbation
        }
    }
}




