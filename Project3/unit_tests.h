#ifndef UNIT_TESTS_H
#define UNIT_TESTS_H
#include <methods.h>
#include <integrand.h>

double func(double x1){
    return x1*x1;
}

double func_2(double r,double t, double p){
    //sqrt(pow(x,2)+pow(y,2)+pow(z^2))
    p;
}

void test_leg(){
    /*test_leg():
    *Unit test for gaussleg. We have found an
    analytical solution to the integral and compare our method in
    one dimention tot he analytical answer.*/
    int N = 15;
    int a = -3;
    int b = 3;
    int tol  = 1E-8;
    double *x = new double [N];
    double *w = new double [N];

    gauleg(a,b,x,w, N);

    double gauss = 0.0;
    for (int i=0;i<N;i++){
       gauss+=w[i]*func(x[i]);
           }
    try{
        if (gauss-18 > tol){
        throw "Warning: The Jacobi rotation method is not accurate and is not giving excpected values, there may be a serious problem";
    }
    }
        catch (const char* msg){
            cerr << msg <<endl;
            }
    delete[] x; delete[] w;
};



void test_lag(){
    /*task3b():
     * Improved Gauss-Legendre/ Gauss-Laguerre.
       We set up weights and grid similar to how we did in task3a(),
       but because we're in polar coordinates so we have to make adjustments.
       We use Gauss-Legendre to find the weights and grid of phi and theta,
       and then we use Gauss-Laguerre to find the weights and grid of r.
       Gauss-Laguerre is defined from [0,inf).*/

    double pi = 3.14159265;
    int n_lag = 15;
    int n_leg = 15;
    double *x1 = new double [n_lag];//Meshgrid for theta
    double *w1 = new double [n_lag];// Weight for theta
    double *x2 = new double [n_leg];//Meshgrid for phi
    double *w2 = new double [n_leg];//Weightfor phi

    int a = 0;
    double b = pi;
    gauleg(a,b,x1,w1, n_leg);// For theta
    gauleg(a,2*b,x2,w2, n_leg); // For phi

    double *xgl = new double [n_lag+1]; // Radius meshgrid
    double *wgl = new double [n_lag+1]; // Radius vekttall


    double alf = 0;
    gauss_laguerre(xgl,wgl, n_lag, alf);//For radiusen

    double int_gausslag = 0.0;
   // cout<< xgl[1]<<endl;
    for (int i=1;i<=n_lag;i++){ //Radius 1
        cout << i << endl;
        for (int j = 1; j<=n_lag; j++){//Radius 2

        for (int k = 0; k<n_leg; k++){//Theta 1
        for (int l = 0; l<n_leg; l++){// Theta 2

        for (int m = 0; m<n_leg; m++){// Phi 1
        for (int n = 0; n<n_leg; n++){// Phi 2

            int_gausslag += wgl[i]
                * func(xgl[i]);
}}}}}}

    cout<< int_gausslag<<endl;
    delete []xgl;
    delete []wgl;
    delete []x1;
    delete []w1;
    delete []x2;
    delete []w2;

}


void test_ans(){
    int tol = 1E-2;
    cout << task3a() << endl;
    try{
        if (abs(0.192765-task3a()) > tol){
        throw "Warning: Not even close to the analytical value";
    }
    }
        catch (const char* msg){
            cerr << msg <<endl;
            }


};

#endif // UNIT_TESTS_H
