#include <iostream>
#include <integrand.h>
#include <methods.h>
#include <cmath>
#include <armadillo>
#include <monte_carlo.h>
#include <unit_tests.h>
#include <omp.h>


using namespace std;
using namespace arma;

    int main(){
        int nthreads = omp_get_num_threads();
        printf("Using %d threads outside parallel loop\n",nthreads);
        int mthreads = omp_get_max_threads();
        printf("There are %d threads available at a time\n",mthreads);
        int M = 9;
        double std;
        int i = 2;
        double _t ;
        double t_max = 5;
        int _x = 10;
        double _sig;
        //test_Laguerre();
/*
    ofstream afile;
    string afilename = "gauss_legendre.dat";
    afile.open(afilename);
    afile << setiosflags(ios::showpoint | ios::uppercase);
         while (_t<t_max){
         i ++;
         cout <<  i <<endl;
         afile << setw(15) << setprecision(8) <<  i;
         afile << setw(20) << setprecision(8) <<  task3a(i,_t);
         afile << setw(20) << setprecision(8) <<  _t<<endl;
    }
     afile.close();
     _t = 0;
     i  = 2;

     ofstream bfile;
     string bfilename = "gauss_laguerre.dat";
     bfile.open(bfilename);
     bfile << setiosflags(ios::showpoint | ios::uppercase);
     while (_t<t_max){
        i++;
        cout <<  i <<endl;
        bfile << setw(15) << setprecision(8) <<  i;
        bfile << setw(20) << setprecision(8) <<  task3b(i,_t);
        bfile << setw(20) << setprecision(8) <<  _t<<endl;
   }

     bfile.close();
     _t = 0;
     //i  = 3;
*/
     ofstream cfile;
     string cfilename = "brute_force_monte_carlo.dat";
     cfile.open(cfilename);
     cfile << setiosflags(ios::showpoint | ios::uppercase);
     //while (_t<t_max){
       // i++;

     for(int i = 4; i<M; i++){
        cout <<  pow(_x,i) <<endl;
        cfile << setw(15) << setprecision(8) <<  pow(_x,i);
        cfile << setw(20) << setprecision(8) <<  Brute_MonteCarlo(pow(_x,i),_t,_sig);
        cfile << setw(20) << setprecision(8) <<  _t;
        cfile << setw(20) << setprecision(8) <<  _sig<<endl;
   }
     cfile.close();
    _t = 0;
    //i  = 3;

     ofstream dfile;
     string dfilename = "importance_sampling_monte_carlo.dat";
     dfile.open(dfilename);
     dfile << setiosflags(ios::showpoint | ios::uppercase);

     //while (_t<t_max) {
        //i++;
       for (int i=4;i<M;i++){
        cout <<  pow(_x,i) <<endl;
        dfile << setw(15) << setprecision(8) <<  pow(_x,i);
        dfile << setw(20) << setprecision(8) <<  Importance_MonteCarlo(pow(_x,i),_t, _sig);
        dfile << setw(20) << setprecision(8) <<  _t;
        dfile << setw(20) << setprecision(8) <<  _sig<<endl;
   }

     dfile.close();
     _t = 0;
     //i  = 3;

     ofstream cpfile;
     string cpfilename = "parallel_brute_force_monte_carlo.dat";
     cpfile.open(cpfilename);
     cpfile << setiosflags(ios::showpoint | ios::uppercase);
     //while (_t<t_max){
        //i++;
       for (int i=4;i<M;i++){
        cout <<  pow(_x,i) <<endl;
        cpfile << setw(15) << setprecision(8) <<  pow(_x,i);
        cpfile << setw(20) << setprecision(8) <<  Parallel_Brute_MonteCarlo(pow(_x,i),_t,_sig);
        cpfile << setw(20) << setprecision(8) <<  _t;
        cpfile << setw(20) << setprecision(8) <<  _sig <<endl;
       }

     cpfile.close();
    _t = 0;
    //i  = 3;

     ofstream dpfile;
     string dpfilename = "parallel_importance_sampling_monte_carlo.dat";
     dpfile.open(dpfilename);
     dpfile << setiosflags(ios::showpoint | ios::uppercase);

     //while (_t<t_max) {
       // i++;
       for (int i=4;i<M;i++){
        cout <<  pow(_x,i) <<endl;
        dpfile << setw(15) << setprecision(8) <<  pow(_x,i);
        dpfile << setw(20) << setprecision(8) <<  Parallel_Importance_MonteCarlo(pow(_x,i),_t, _sig);
        dpfile << setw(20) << setprecision(8) <<  _t;
        dpfile << setw(20) << setprecision(8) <<  _sig<<endl;
   }

     dpfile.close();


    return 0;
    }
