#include <iostream>
#include <integrand.h>
#include <methods.h>
#include <cmath>
#include <armadillo>
#include <monte_carlo.h>
#include <unit_tests.h>
//#include <omp.h> //The Open MPI package


using namespace std;
using namespace arma;

    int main(){

    //task3a();
    //task3b();
    //test_ans();
    double integralening ;
    double _std;
    Importance_MonteCarlo(10000000, -2,2, integralening,_std);
    cout <<integralening <<endl;
    return 0;
    }
