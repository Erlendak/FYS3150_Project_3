#ifndef UNIT_TESTS_H
#define UNIT_TESTS_H
#include <methods.h>
#include <integrand.h>

double func(double x1){
    return x1*x1;
}


void test_leg(){
    /*test_leg:
    *Unit test for gaussleg. Set*/
    int N = 26;
    int a = -3;
    int b = 3;
    int tol  = 1E-8;
    double *x = new double [N];
    double *w = new double [N];
//  Function from lib, sets up mesh grid and weights.
    gauleg(a,b,x,w, N);

//   Gauss-Legendre method
    double gauss = 0.0;
    for (int i=0;i<N;i++){
       gauss+=w[i]*func(x[i]);
           }
    try{
        if (gauss-18 > tol){
        throw "Warning: The Jacobi rotation method is not accurate and is not giving excpected values, there may be a serious problem";
    }
    }
        catch (const char* msg){
            cerr << msg <<endl;
            }
};

void test_ans(){
    int tol = 1E-6;
    try{
        if (abs(0.192765-task3a()) > tol){
        throw "Warning: Bot even close to the analytical value";
    }
    }
        catch (const char* msg){
            cerr << msg <<endl;
            }


};

#endif // UNIT_TESTS_H
