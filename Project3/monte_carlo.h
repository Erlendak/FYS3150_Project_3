#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <iostream>
#include <integrand.h>
#include <methods.h>
#include <cmath>
#include <armadillo>

using namespace std;

void Brute_MonteCarlo(int n, double a, double b, double  &integral, double  &std){
        double * x = new double [n];
        double x1, x2, y1, y2, z1, z2, f;
        double mc = 0.0;
        double sigma = 0.0;
        int i;
        double jacob = pow((b-a),6);

        #pragma omp parallel for reduction(+:mc)  private (i, x1, x2, y1, y2, z1, z2, f)
        for (i = 0; i < n; i++){
                x1 = rand()*(b-a)+a;
                x2 = rand()*(b-a)+a;
                y1 = rand()*(b-a)+a;
                y2 = rand()*(b-a)+a;
                z1 = rand()*(b-a)+a;
                z2 = rand()*(b-a)+a;
                f = func_cartesian(x1, y1, z1, x2, y2, z2);
                mc += f;
                x[i] = f;
        }
        mc = mc/( n );
        #pragma omp parallel for reduction(+:sigma)  private (i)
        for (i = 0; i < n; i++){
                sigma += (x[i] - mc)*(x[i] - mc);
        }
        double _n = n;
        sigma = sigma*jacob/(_n );
        std = sqrt(sigma)/sqrt(_n);
        integral = mc*jacob;
        cout<< integral<<endl;
        delete [] x;
} // end function for the integrand

#endif // MONTE_CARLO_H

