#include <iostream>
#include <gauss_legendre.h>
#include <gauss_lag.h>
#include <cmath>

using namespace std;

int main()
{
//task3a();
double pi = 3.14159265;
int N = 10;
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
double int_gausslag = 0.0;
for (int i=0;i<N;i++){
    cout << i << endl;
    for (int j = 0; j<N; j++){
    for (int k = 0; k<pi; k++){
    for (int l = 0; l<pi; l++){
    for (int m = 0; m< 2*pi; m++){
    for (int n = 0; n < 2*pi; n++){

    int_gausslag += wgl[i]*wgl[k]*wgl[m]*wgl[j]*wgl[l]*wgl[n]
            *func_spherical(xgl[i],xgl[k],xgl[m],xgl[j],xgl[l],xgl[n]);

//cout<< int_gausslag<<endl;
}}}}}}
delete []xgl;
delete []wgl;
cout<< int_gausslag<<endl;
return 0;
}
