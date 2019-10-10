#include <iostream>
#include <gauss_legendre.h>
#include <gauss_lag.h>

using namespace std;

int main()
{
//task3a();

int n = 4;
double a = 0;
double b= 1;
//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method
// Gauss-Laguerre is old-fashioned translation of F77 --> C++
// arrays start at 1 and end at n
double *xgl = new double [n+1];
double *wgl = new double [n+1];
        // These arrays are used for improved Gauss-Legendre, mapping of
        // x \in [-1,1] to x \in [0, infinity)

   //   set up the mesh points and weights
   //   set up the mesh points and weights and the power of x^alf
double alf = 2.0;
gauss_laguerre(xgl,wgl, n, alf);

//   evaluate the integral with the Gauss-Laguerre method
//   Note that we initialize the sum
double int_gausslag = 0.;
for ( int i = 1;  i <= n; i++){
    int_gausslag += wgl[i]*sin(xgl[i]);
}
cout<< int_gausslag<<endl;
}
