#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include <gauss_legendre.h>
#include <gauss_lag.h>
double func(double x);

//     Main function begins here
int motecarlo(int N)
{
     double MCint, MCintsqr2, fx, Variance;
     MCint = MCintsqr2=0.;
     double invers_period = 1./RAND_MAX; // initialise the random number generator
     srand(time(NULL));  // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method
     for ( int i = 1;  i <= N; i++){
             for (int j = 1; j<=N; j++){//Radius 2
                for (int k = 0; k<N; k++){//Theta 1
                     for (int l = 0; l<N; l++){// Theta 2
                           for (int m = 0; m<N; m++){// Phi 1
                                for (int n = 0; n<N; n++){// Phi 2
  // obtain a floating number x in [0,1]
           double x = double(rand())*invers_period;
           fx = func(x);
           MCint += fx;
           MCintsqr2 += fx*fx;
     }
     MCint = MCint/((double) N );
     MCintsqr2 = MCintsqr2/((double) N );
     double variance=MCintsqr2-MCint*MCint;
//   final output
     cout << " variance= " << variance << " Integral = " << MCint << " Exact= " << M_PI << endl;
}  // end of main program
// this function defines the function to integrate
double func(double x)
{
  double value;
  value = 4/(1.+x*x);
  return value;
}
#endif // MONTE_CARLO_H
