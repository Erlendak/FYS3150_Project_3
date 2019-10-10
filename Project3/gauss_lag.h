#ifndef GAUSS_LAG_H
#define GAUSS_LAG_H


#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10

using namespace std;


double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

//  Note that you need to call it with a given value of alpha,
// called alf here. This comes from x^{alpha} exp(-x)

void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}
// end function gaulag

// end function gammln
#undef EPS
#undef MAXIT



double func_spherical(double r1,double theta1, double phi1, double r2, double theta2, double phi2){
    /*func(double x1,double y1, double z1, double x2, double y2, double z2):
     * Function that sets up the integral we're going to integrate.
      This is a six-dimensional integral.*/

    double alpha = 2.0;
    double tol = 1E-10; //Tolerance: |r1-r2|=0 can cause problems, so have to check.
    //double deno = sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));

    //We make an if test, because we have to account for when deno goes towards zero,
    //if we don't, our expression will go against infinity.
    double beta = cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1-phi2);
    double deno = sqrt((r1*r1) + (r2*r2) - 2*r1*r2*beta);
    if(deno < tol){
        return 0;
    }
    else{
    double exp1 = -2*alpha*sqrt(r1);
    double exp2 = -2*alpha*sqrt(r2);
    return exp(exp1+exp2)/deno;
    }
}



#endif // GAUSS_LAG_H
