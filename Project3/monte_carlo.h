#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <iostream>
#include <integrand.h>
#include <methods.h>
#include <cmath>
#include <armadillo>
#include <random>

using namespace std;

void Brute_MonteCarlo(int n, double a, double b, double  &integral, double  &std){
        random_device ran;
        mt19937_64 gen(ran());

        uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
        double * x = new double [n];
        double x1, x2, y1, y2, z1, z2, f;
        double mc = 0.0;
        double sigma = 0.0;
        int i;
        double jacob = pow((b-a),6);


        #pragma omp parallel for reduction(+:mc)  private (i, x1, x2, y1, y2, z1, z2, f)
        for (i = 0; i < n; i++){
                x1 = RandomNumberGenerator(gen)*(b-a) + a;
                x2 = RandomNumberGenerator(gen)*(b-a) +a;
                y1 = RandomNumberGenerator(gen)*(b-a) +a;
                y2 = RandomNumberGenerator(gen)*(b-a) +a;
                z1 = RandomNumberGenerator(gen)*(b-a) +a;
                z2 = RandomNumberGenerator(gen)*(b-a) +a;
                //cout <<  x1 << " " <<  x2 << " "<<  y1 << " "<<  y2 << " "<<  z1 << " "<<  z2 << endl;
                f = func_cartesian(x1, y1, z1, x2, y2, z2);
                mc += f;
                //cout << f << endl;
                x[i] = f;
        }
        mc = mc/( (double) n );
        #pragma omp parallel for reduction(+:sigma)  private (i)
        for (i = 0; i < n; i++){
                sigma += (x[i] - mc)*(x[i] - mc);
        }
        double _n = n;
        sigma = sigma*jacob/((double)_n );
        std = sqrt(sigma)/sqrt((double)_n);
        integral = mc*jacob;
        cout<< integral<<endl;
        delete [] x;
} // end function for the integrand

double guass_deviant(long *idum){
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if ( idum < 0) iset =0;
  if (iset == 0) {
    do {
      v1 = 2.*ran0(idum) -1.0;
      v2 = 2.*ran0(idum) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset =0;
    return gset;
  }
}


void Importance_MonteCarlo(int n, double a, double b, double  &integral, double  &std){
        random_device ran;
        mt19937_64 gen(ran());

        uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
        double * x = new double [n];
        double x1, x2, y1, y2, z1, z2, f;
        double mc = 0.0;
        double variance = 0.0 ;
        long idum = - 1;
        double sigma = 0.0;
        int i;
        double jacob = pow((b-a),6);


        #pragma omp parallel for reduction(+:mc)  private (i, x1, x2, y1, y2, z1, z2, f)
        for (i = 0; i < n; i++){
                x1 = guass_deviant(idum)*((b-a) + a); //*RandomNumberGenerator(gen);
                x2 = RandomNumberGenerator(gen)*(b-a) +a;
                y1 = RandomNumberGenerator(gen)*(b-a) +a;
                y2 = RandomNumberGenerator(gen)*(b-a) +a;
                z1 = RandomNumberGenerator(gen)*(b-a) +a;
                z2 = RandomNumberGenerator(gen)*(b-a) +a;
                //cout <<  x1 << " " <<  x2 << " "<<  y1 << " "<<  y2 << " "<<  z1 << " "<<  z2 << endl;
                f = func_cartesian(x1, y1, z1, x2, y2, z2);
                mc += f;
                //cout << f << endl;
                x[i] = f;
        }
        mc = mc/( (double) n );
        #pragma omp parallel for reduction(+:sigma)  private (i)
        for (i = 0; i < n; i++){
                sigma += (x[i] - mc)*(x[i] - mc);
        }
        double _n = n;
        sigma = sigma*jacob/((double)_n );
        std = sqrt(sigma)/sqrt((double)_n);
        integral = mc*jacob;
        cout<< integral<<endl;
        delete [] x;
}
#endif // MONTE_CARLO_H

