#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H

#endif // GAUSS_LEGENDRE_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <time.h>
#include <stdlib.h>
#include "armadillo"

using namespace arma;
using namespace std;

double func(double x1,double y1, double z1, double x2, double y2, double z2){
    double alpha = 2.0;
    double tol = 1E-10;
    double deno = sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));

    if(deno < tol){
        return 0;
    }else{
        double exp1 = -2*alpha*sqrt((x1*x1)+(y1*y1)+(z1*z1));
        double exp2 = -2*alpha*sqrt((x2*x2)+(y2*y2)+(z2*z2));
        return exp(exp1+exp2)/deno;
    }


}

