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

void Brute_MonteCarlo(int n, double a, double b, double  &integral, double  &std){
        /*Brute_Montecarlo(int n, double a, double b, double  &integral, double  &std):
         * Brute force way of implementing the Monte Carlo integration method. We have
         * a random number generation making numbers for us in the interval [a,b] and
         * then we use them in our function. We do this several times, and the
         * sum of the function should be an approximation to our integral. Here we use
         * uniform distribution and cartesian coordinates in our function.
        */
        //Random number generation
        random_device ran;
        mt19937_64 gen(ran());

        clock_t start, finish;

        //Uniform distribution
        uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0); // (a,b)
        double *x = new double [n];
        double x1, x2, y1, y2, z1, z2, f;
        double mc = 0.0;
        double sigma = 0.0;
        int i;
        //Jacobi determinant for uniform distribution
        double jacob = pow((b-a),6);

        start = clock();
        //#pragma omp parallel for reduction(+:mc)  private (i)
        for (i = 0; i < n; i++){
                x1 = RandomNumberGenerator(gen)*(b-a) + a; //RandomNumberGenerator(gen)
                x2 = RandomNumberGenerator(gen)*(b-a) +a;
                y1 = RandomNumberGenerator(gen)*(b-a) +a;
                y2 = RandomNumberGenerator(gen)*(b-a) +a;
                z1 = RandomNumberGenerator(gen)*(b-a) +a;
                z2 = RandomNumberGenerator(gen)*(b-a) +a;
                f = func_cartesian(x1, y1, z1, x2, y2, z2);
                mc += f;
                x[i] = f;
        }
        mc = mc/( (double) n );
        //#pragma omp parallel for reduction(+:sigma)  private (i)
        for (i = 0; i < n; i++){
                sigma += (x[i] - mc)*(x[i] - mc);
        }
        finish = clock();
        cout << "Algorithm time:" <<((((double)finish - (double)start)/CLOCKS_PER_SEC)) << "s"<<endl;
        double _n = n;
        sigma = sigma*jacob/((double)_n ); //Variance
        std = sqrt(sigma)/sqrt((double)_n);//Standard Deviation
        integral = mc*jacob;
        cout<< integral<<endl;
        delete [] x;
}




void Importance_MonteCarlo(int n, double a, double b, double  &integral, double  &std){
        /* Importance_MonteCarlo(int n, double a, double b, double  &integral, double  &std):
           In this improved version of Monte Carlo we use spherical coordinates instead of cartesian coordinates.
           We also use exponential distribution to find the numbers our interval for the radius, and
           uniform distribution for the angles. On top of this, because we're using exponential distribution
           we can simplify parts of our integrand because it gets absorbed.*/
        clock_t start, finish;
        double pi = 3.14159265;
        //Random number generation
        random_device ran;
        mt19937_64 gen(ran());

        //Exponential distribution for r
        exponential_distribution<double> Exponential_R(-4);

        //Uniform distribution for angles
        uniform_real_distribution<double> UniformTheta(0,pi);
        uniform_real_distribution<double> UniformPhi(0,2*pi);

        double * x = new double [n];
        double r1, r2, t1, t2, p1, p2, f;
        double mc = 0.0;
        double sigma = 0.0;

        //Jacobi determinant for exponential distribution
        double jacob = 4*pi*pi*pi*pi/16;

        int i;
        start = clock();
        //omp_set_num_threads(2);
        //#pragma omp parallel for reduction(+:mc)  private (i)
        for (i = 0; i < n; i++){
                r1 = Exponential_R(gen);
                r2 = Exponential_R(gen);
                t1 = UniformTheta(gen);
                t2 = UniformTheta(gen);
                p1 = UniformPhi(gen);
                p2 = UniformPhi(gen);

                f = func_spherical_mc(r1, t1, p1, r2, t2, p2);
                mc += f;
                x[i] = f;
                //printf("Using %d threads\n",omp_get_num_threads());
        }
        //int nthreads = omp_get_num_threads();
        //printf("Using %d threads\n",nthreads);
        //int mthreads = omp_get_max_threads();
        //printf("There are %d threads available at a time\n",mthreads);
        mc = mc/( (double) n );
        //#pragma omp parallel for reduction(+:sigma)  private (i)
        for (i = 0; i < n; i++){
                sigma += (x[i] - mc)*(x[i] - mc);
        }
        finish = clock();
        cout << "Algorithm time:" <<((((double)finish - (double)start)/CLOCKS_PER_SEC)) << "s"<<endl;
        double _n = n;
        sigma = sigma*jacob/((double)_n ); //Variance
        std = sqrt(sigma)/sqrt((double)_n);//Standard deviation
        integral = mc*jacob;
        cout<< integral<<endl;
        delete [] x;
}

#endif // MONTE_CARLO_H
