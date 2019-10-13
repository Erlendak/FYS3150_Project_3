#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <time.h>
#include <stdlib.h>
#include "armadillo"
#include <iostream>
#include <lib.h>
#include <methods.h>

using namespace arma;
using namespace std;

double func_cartesian(double x1,double y1, double z1, double x2, double y2, double z2){
    /*func(double x1,double y1, double z1, double x2, double y2, double z2):
     * Function that sets up the integral we're going to integrate.
      This is a six-dimensional integral.*/

    double alpha = 2.0;
    double tol = 1E-10; //Tolerance: |r1-r2|=0 can cause problems, so have to check.
    double deno = sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));

    //We make an if test, because we have to account for when deno goes towards zero,
    //if we don't, our expression will go against infinity.
    if(deno < tol){
        return 0;
    }else{
        double exp1 = -2*alpha*sqrt((x1*x1)+(y1*y1)+(z1*z1));
        double exp2 = -2*alpha*sqrt((x2*x2)+(y2*y2)+(z2*z2));
        return exp(exp1+exp2)/deno;
    }
}

double func_spherical(double r1,double theta1, double phi1, double r2, double theta2, double phi2){
    double alpha = 2.0;
    double tol = 1E-10; //Tolerance: Too small denominator can cause problems, so have to check

    //We make an if test, because we have to account for when deno goes towards zero,
    //if we don't, our expression will go against infinity.
    double beta = cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2);
    double deno = sqrt(abs(pow(r1,2) + pow(r2,2) - (2*r1*r2*beta)));
    if(deno < tol){
        return 0;
    }
    else{
    double exp1 = -2*alpha*r1;
    double exp2 = -2*alpha*r2;
    return (sin(theta1)*sin(theta2)*pow(r1,2)*pow(r2,2)*exp(exp1+exp2))/deno;
    }
}

void task3a() {
    int N = 27;
    int a = -3;
    int b = 3;
    double *x = new double [N];
    double *w = new double [N];
//  Function from lib, sets up mesh grid and weights.
    gauleg(a,b,x,w, N);

//   Gauss-Legendre method
    double gauss = 0.0;
    for (int i=0;i<N;i++){
        cout << i << endl;
        for (int j = 0;j<N;j++){
        for (int k = 0;k<N;k++){
        for (int l = 0;l<N;l++){
        for (int m = 0;m<N;m++){
        for (int n = 0;n<N;n++){
       gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]*func_cartesian(x[i],x[j],x[k],x[l],x[m],x[n]);
           }}}}}
    }
    cout << gauss<< endl;
    delete [] x;
    delete [] w;
}

void task3b(){
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
}



#endif // GAUSS_LEGENDRE_H
