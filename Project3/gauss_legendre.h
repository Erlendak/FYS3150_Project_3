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
    };
};

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





#endif // GAUSS_LEGENDRE_H
