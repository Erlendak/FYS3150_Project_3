#ifndef UNIT_TESTS_H
#define UNIT_TESTS_H
#include <methods.h>
#include <integrand.h>
#include <monte_carlo.h>

double func(double x1){
    return x1*x1;
}

void test_leg(){
    /*test_leg():
    *Unit test for gaussleg. We have found an
    analytical solution to the integral and compare our method in
    one dimention tot he analytical answer.*/

    int N = 27;
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

void test_laguerre(){
    /*
    Unit test for our Gauss Laguerre. We test if our implementiaion off Gauss Laguere will calculate the integral
    with 20 integration points to our expected solution 0.19478557.
    */

    double expected = 0.19478557;
    double t;
    double tested = task3b(20,t);
    double tol = 0.00001;
    try{
        if (abs(abs(tested)-abs(expected)) > tol){
        throw "Warning: Integration using Gauss Laguerre method is not as accurate as expected, there may be a serious problem that needs to be checked up upon.";
    }
    }
        catch (const char* msg){
            cerr << msg <<endl;
            }
};

void test_monte_carlo(){
    /*
    Unit test for our all Monte Carlo methods. We are testing that Monte Carlo methods gets
    sufficient varriance to make a accurate integral.
    */
    double t;
    double vari;
    double tol = 0.1;

    Brute_MonteCarlo(10000000,t,vari);


    try{
        if (vari > tol){
        throw "Warning: Integration using Brute Force Monte Carlo method does not give a satisfying variance, there may be a serious problem that needs to be checked up upon.";
    }
    }
        catch (const char* msg){
            cerr << msg <<endl;
            }

    Importance_MonteCarlo(10000000,t,vari);


    try{
        if (vari > tol){
        throw "Warning: Integration using Importance Sampling method method does not give a satisfying variance, there may be a serious problem that needs to be checked up upon.";
    }
    }
        catch (const char* msg){
            cerr << msg <<endl;
     }

    Parallel_Brute_MonteCarlo(10000000,t,vari);
    try{
        if (vari > tol){
        throw "Warning: Integration using Parallell Brute Force Monte Carlo method does not give a satisfying variance, there may be a serious problem that needs to be checked up upon.";
    }
    }
        catch (const char* msg){
            cerr << msg <<endl;
            }

    Parallel_Importance_MonteCarlo(10000000,t,vari);


    try{
        if (vari > tol){
        throw "Warning: Integration using Parallell Importance Sampling method method does not give a satisfying variance, there may be a serious problem that needs to be checked up upon.";
    }
    }
        catch (const char* msg){
            cerr << msg <<endl;
            }
};
void tests(){
/*
Calls all the above tests at once. For a better looking main.
*/
    test_leg();
    test_laguerre();
    test_monte_carlo();
}

#endif // UNIT_TESTS_H
