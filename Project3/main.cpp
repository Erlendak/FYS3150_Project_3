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

        // These arrays are used for improved Gauss-Legendre, mapping of
        // x \in [-1,1] to x \in [0, infinity)

   //   set up the mesh points and weights
   //   set up the mesh points and weights and the power of x^alf

double *x1 = new double [N];//Meshgrid for theta
double *w1 = new double [N];// vekttall for theta
double *x2 = new double [N];//Meshgrid for phi
double *w2 = new double [N];//vekttall for phi
int a = 0;
double b = pi;
gauleg(a,b,x1,w1, N);// For theta
gauleg(a,2*b,x2,w2, N); // For phi
double *xgl = new double [N+1]; // radiusens meshgrid
double *wgl = new double [N+1]; // radiusen vekttall
double alf = 2.0;
gauss_laguerre(xgl,wgl, N, alf);//For radiusen

//vec _theta = linspace(0,pi,N+1);
//vec _phi =linspace(0,2*pi,N+1);
double int_gausslag = 0.0;
for (int i=1;i<=N;i++){ //Radius 1
    cout << i << endl;
    for (int j = 1; j<=N; j++){//Radius 2

    for (int k = 0; k<N; k++){//Theta 1
    for (int l = 0; l<N; l++){// Theta 2

    for (int m = 0; m<N; m++){// Phi 1
    for (int n = 0; n<N; n++){// Phi 2

        int_gausslag += wgl[i]*wgl[j]*w1[k]*w1[l]*w2[m]*w2[n]
            * func_spherical(xgl[i], x1[k], x2[m],xgl[j],x1[l],x2[n]);

//cout<< int_gausslag<<endl;
}}}}}}
delete []xgl;
delete []wgl;
cout<< int_gausslag<<endl;
return 0;
}
