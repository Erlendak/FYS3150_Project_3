#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <iostream>
#include <integrand.h>
#include <methods.h>
#include <cmath>
#include <armadillo>
#include <random>
#include <omp.h>
#include <cstdio>
using namespace std;

double Brute_MonteCarlo(int n, double  &_t, double  &_sig){
    /* Brute_MonteCarlo(int n, double  &_t, double  &_sig):
       In the Brute Force Monte Carlo method we use the same integrand as we did in Gauss-Legendre,
       and we use it in cartesian coordiantes. We also use uniform distribution and a random
       number generator to find random numbers within our interval. We put the random numbers
       into our integrand. We do this several times and the sum times the Jacobi determinant
       is the answer to our integral.*/

        clock_t start, finish;
        //Integration intervall
        double a = -3;
        double b = 3;
        //Random generators
        random_device ran;
        mt19937_64 gen1(ran());
        mt19937_64 gen2(ran());
        mt19937_64 gen3(ran());
        mt19937_64 gen4(ran());
        mt19937_64 gen5(ran());
        mt19937_64 gen6(ran());
        //Distribution
        uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0); // (a,b)
        double * x = new double [n];
        double x1, x2, y1, y2, z1, z2, f;
        double mc = 0.0;
        double sigma = 0.0;
        int i;
        //Jacobi determinant
        double jacob = pow((b-a),6);

        start = clock();

        for (i = 0; i < n; i++){
                x1 = RandomNumberGenerator(gen1)*(b-a) + a; //RandomNumberGenerator(gen)
                x2 = RandomNumberGenerator(gen2)*(b-a) +a;
                y1 = RandomNumberGenerator(gen3)*(b-a) +a;
                y2 = RandomNumberGenerator(gen4)*(b-a) +a;
                z1 = RandomNumberGenerator(gen5)*(b-a) +a;
                z2 = RandomNumberGenerator(gen6)*(b-a) +a;
                f = func_cartesian(x1, y1, z1, x2, y2, z2);
                mc += f;
                x[i] = f;
        }
        finish = clock();
        _t = ((((double)finish - (double)start)/CLOCKS_PER_SEC));
        mc = mc/( (double) n );

        for (i = 0; i < n; i++){
                sigma += (x[i] - mc)*(x[i] - mc);
        }
        double _n = n;
        sigma = sigma*jacob/((double)_n ); //Variance
        _sig = sigma;
        //std = sqrt(sigma)/sqrt((double)_n);
        double integral = mc*jacob;
        delete [] x;
        return( integral);
}



double Importance_MonteCarlo(int n, double  &_t, double & _sig){
    /* Importance_MonteCarlo(int n, double  &_t, double & _sig):
       In this improved version of Monte Carlo we use spherical coordinates instead of cartesian coordinates.
       We also use exponential distribution to find the numbers our interval for the radius, and
       uniform distribution for the angles. On top of this, because we're using exponential distribution
       we can simplify parts of our integrand because it gets absorbed.*/

        double pi = 3.14159265;
        clock_t start, finish;
        //Random generators
        random_device ran;
        mt19937_64 gen1(ran());
        mt19937_64 gen2(ran());
        mt19937_64 gen3(ran());
        mt19937_64 gen4(ran());
        mt19937_64 gen5(ran());
        mt19937_64 gen6(ran());
        //Distributions
        exponential_distribution<double> Exponential_R(4);
        uniform_real_distribution<double> UniformTheta(0,pi);
        uniform_real_distribution<double> UniformPhi(0,2*pi);
        double * x = new double [n];
        double r1, r2, t1, t2, p1, p2, f;
        double mc = 0.0;

        double sigma = 0.0;

        //Jacobi determinant
        double jacob = 4*pi*pi*pi*pi /16;

        int i;
        start = clock();

        for (i = 0; i < n; i++){
                //Define r with exponential distribution
                r1 = Exponential_R(gen1);
                r2 = Exponential_R(gen2);
                //Angles with uniform distribution
                t1 = UniformTheta(gen3);
                t2 = UniformTheta(gen4);
                p1 = UniformPhi(gen5);
                p2 = UniformPhi(gen6);

                f = func_polar_mc(r1, t1, p1, r2, t2, p2);
                mc += f;
                x[i] = f;
        }
        finish = clock();
        _t = ((((double)finish - (double)start)/CLOCKS_PER_SEC));
        mc = mc/( (double) n );

        for (i = 0; i < n; i++){
                sigma += (x[i] - mc)*(x[i] - mc);
        }
        sigma = sigma*jacob/((double)n ); //Variance
        _sig = sigma;
        //std = sqrt(sigma)/sqrt((double)n);
        double integral = mc*jacob;
        delete [] x;
        return( integral);
}

double Parallel_Brute_MonteCarlo(int n, double  &_t, double  &_sig){
    /* Parallel_Brute_MonteCarlo(int n, double  &_t, double  &_sig):
        We parallellize the Brute Force Monte Carlo method with OpenMP
        so we can use more cores and threads when we run the code.
        This should give us a speed-up.*/
        clock_t start, finish;
        //Integration interval
        double a = -3;
        double b = 3;
        //Random generators
        random_device ran;
        mt19937_64 gen1(ran()+omp_get_thread_num());
        mt19937_64 gen2(ran()+omp_get_thread_num());
        mt19937_64 gen3(ran()+omp_get_thread_num());
        mt19937_64 gen4(ran()+omp_get_thread_num());
        mt19937_64 gen5(ran()+omp_get_thread_num());
        mt19937_64 gen6(ran()+omp_get_thread_num());
        //Distributions
        uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0); // (a,b)
        double * x = new double [n];
        double x1, x2, y1, y2, z1, z2, f;
        double mc = 0.0;
        double sigma = 0.0;
        int i;
        //Jacobi determinant
        double jacob = pow((b-a),6);

        start = clock();
        #pragma omp parallel for reduction(+:mc) private(i, x1, x2, y1,y2, z1, z2, f)
        for (i = 0; i < n; i++){
                x1 = RandomNumberGenerator(gen1)*(b-a) + a; //RandomNumberGenerator(gen)
                x2 = RandomNumberGenerator(gen2)*(b-a) +a;
                y1 = RandomNumberGenerator(gen3)*(b-a) +a;
                y2 = RandomNumberGenerator(gen4)*(b-a) +a;
                z1 = RandomNumberGenerator(gen5)*(b-a) +a;
                z2 = RandomNumberGenerator(gen6)*(b-a) +a;
                f = func_cartesian(x1, y1, z1, x2, y2, z2);
                mc += f;
                x[i] = f;

        }
        finish = clock();
        _t = ((((double)finish - (double)start)/CLOCKS_PER_SEC));
        mc = mc/( (double) n );
        #pragma omp parallel for reduction(+:sigma)  private (i)
        for (i = 0; i < n; i++){
                sigma += (x[i] - mc)*(x[i] - mc);
        }
        double _n = n;
        sigma = sigma*jacob/((double)_n ); //Variance
        _sig = sigma;
        //std = sqrt(sigma)/sqrt((double)_n);
        double integral = mc*jacob;
        delete [] x;
        return( integral);
}



double Parallel_Importance_MonteCarlo(int n, double  &_t, double & _sig){
    /*Parallel_Importance_MonteCarlo(int n, double  &_t, double & _sig)
      We parallellize the Importance Monte Carlo method with OpenMP
      so we can use more cores and threads when we run the code.
      This should give us a speed-up.*/
        double pi = 3.14159265;
        clock_t start, finish;
        //Random generator
        random_device ran;
        mt19937_64 gen1(ran()+omp_get_thread_num());
        mt19937_64 gen2(ran()+omp_get_thread_num());
        mt19937_64 gen3(ran()+omp_get_thread_num());
        mt19937_64 gen4(ran()+omp_get_thread_num());
        mt19937_64 gen5(ran()+omp_get_thread_num());
        mt19937_64 gen6(ran()+omp_get_thread_num());
        //Distribution
        exponential_distribution<double> Exponential_R(4);
        uniform_real_distribution<double> UniformTheta(0,pi);
        uniform_real_distribution<double> UniformPhi(0,2*pi);
        double * x = new double [n];
        double r1, r2, t1, t2, p1, p2, f;
        double mc = 0.0;

        double sigma = 0.0;

        double jacob = 4*pi*pi*pi*pi /16;


        int i;
        start = clock();

         #pragma omp parallel for reduction(+:mc) private(i, r1, r2, t1, t2, p1, p2, f)
        for (i = 0; i < n; i++){
                //Find radius with exponential distribution
                r1 = Exponential_R(gen1);
                r2 = Exponential_R(gen2);
                //Find angles with uniform distribution
                t1 = UniformTheta(gen3);
                t2 = UniformTheta(gen4);
                p1 = UniformPhi(gen5);
                p2 = UniformPhi(gen6);

                f = func_polar_mc(r1, t1, p1, r2, t2, p2);
                mc += f;
                x[i] = f;
        }
        finish = clock();
        _t = ((((double)finish - (double)start)/CLOCKS_PER_SEC));
        mc = mc/( (double) n );
        #pragma omp parallel for reduction(+:sigma)  private (i)
        for (i = 0; i < n; i++){
                sigma += (x[i] - mc)*(x[i] - mc);
        }
        sigma = sigma*jacob/((double)n);//Variance
        _sig = sigma;
        //std = sqrt(sigma)/sqrt((double)n);
        double integral = mc*jacob;
        delete [] x;
        return( integral);
}



#endif // MONTE_CARLO_H

