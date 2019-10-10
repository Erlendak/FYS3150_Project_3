#include <iostream>
#include <gauss_legendre.h>
#include <gauss_lag.h>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
{
//task3a();
double pi = 3.14159265;

int N = 20;
//double dpi = pi/N;
//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method
// Gauss-Laguerre is old-fashioned translation of F77 --> C++
// arrays start at 1 and end at n
double *xgl = new double [N+1];
double *wgl = new double [N+1];
        // These arrays are used for improved Gauss-Legendre, mapping of
        // x \in [-1,1] to x \in [0, infinity)

   //   set up the mesh points and weights
   //   set up the mesh points and weights and the power of x^alf
double alf = 2.0;
gauss_laguerre(xgl,wgl, N, alf);

//   evaluate the integral with the Gauss-Laguerre method
//   Note that we initialize the sum

vec _theta = linspace(0,pi,N+1);
vec _phi =linspace(0,2*pi,N+1);
double int_gausslag = 0.0;
for (int i=1;i<=N;i++){
    cout << i << endl;
    for (int j = 1; j<=N; j++){
    for (int k = 1; k<=N; k++){
    for (int l = 1; l<=N; l++){
    for (int m = 1; m<=N; m++){
    for (int n = 1; n<=N; n++){

    int_gausslag += wgl[i]*wgl[k]*wgl[m]*wgl[j]*wgl[l]*wgl[n]
            *func_spherical(i,_theta(k),_phi(m),j,_theta(l),_phi(n));

//cout<< int_gausslag<<endl;
}}}}}}
delete []xgl;
delete []wgl;
cout<< int_gausslag<<endl;
return 0;
}
