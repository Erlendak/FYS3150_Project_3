#include <iostream>
#include <integrand.h>
#include <methods.h>
#include <cmath>
#include <armadillo>
#include <monte_carlo.h>


using namespace std;
using namespace arma;

    int main(){
    //cout << func_polar_lag(-1,2,1, 0.5,2 ,-1)<<endl:
    //task3a();
    //task3b();
    double integralening ;
    double _std;
    Importance_MonteCarlo(10000000, -3,3, integralening,_std);
    cout <<integralening <<endl;
    return 0;
    }
