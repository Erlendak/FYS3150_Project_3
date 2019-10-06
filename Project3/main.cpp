#include <iostream>
#include <gauss_legendre.h>
#include <lib.h>

using namespace std;

int main()
{
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
    return 0;
}
