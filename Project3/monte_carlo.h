#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <iostream>
#include <integrand.h>
#include <methods.h>
#include <cmath>
#include <armadillo>

using namespace std;

double brute_force_MC(double *);
//     Main function begins here

void brute(int n){

     double a = 0; // Kunne brukt indekser for forskjellige grense verdier.
     double b = 3;
     double jacobidet = 1;
     for ( int i = 0;  i < 6; i++){
         jacobidet = jacobidet * (b- a) ; // \prod_{i=0}^{8} b_i - a_i
     }

     double x[6], y, fx;
     double int_mc = 1.;  double variance = 0.;
     double sum_sigma= 0. ; long idum=-1 ;
     double length = 5.; // Gjelder for eksempel

     for ( int i = 1;  i <= n; i++){ // Monte Carlo integrasjons loopen
       for (int j = 0; j< 6; j++) { // Setter opp virkårlige tall for dimensjonene
           x[j]=-length+2*length*ran0(&idum); // Må skrives om for å passe hvor integrasjon?
       }
       fx= func_cartesian(x[0],x[1],x[2],x[3],x[4],x[5]); // Kan vi bruke integranten fra tidligere?
       int_mc += fx;
       sum_sigma += fx*fx;
     }
     int_mc = int_mc/((double) n );
     sum_sigma = sum_sigma/((double) n );
     variance=sum_sigma-int_mc*int_mc;


     cout <<  jacobidet*int_mc<< endl;
     cout <<  sum_sigma<< endl;
     }
  // end of main program

// this function defines the integrand to integrate


double  brute_force_MC(double *x)
{
// evaluate the different terms of the exponential
   double xx=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
   double yy=x[3]*x[3]+x[4]*x[4]+x[5]*x[5];
   double xy=pow((x[0]-x[3]),2)+pow((x[1]-x[4]),2)+pow((x[2]-x[5]),2);
   return exp(-xx-yy)*xy;
} // end function for the integrand

#endif // MONTE_CARLO_H
