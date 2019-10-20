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

        double std;
        int i = 2;
        double _t ;
        double t_max = 120;


    ofstream afile;
    string afilename = "gauss_legendre.dat";
    afile.open(afilename);
    afile << setiosflags(ios::showpoint | ios::uppercase);
         while (_t<t_max){
         i ++;
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
        bfile << setw(15) << setprecision(8) <<  i;
        bfile << setw(20) << setprecision(8) <<  task3b(i,_t);
        bfile << setw(20) << setprecision(8) <<  _t<<endl;
   }

     bfile.close();
     _t = 0;
     i  = 3;

     ofstream cfile;
     string cfilename = "brute_force_monte_carlo.dat";
     cfile.open(cfilename);
     cfile << setiosflags(ios::showpoint | ios::uppercase);
     while (_t<t_max){
        i++;
        cfile << setw(15) << setprecision(8) <<  pow(2,i);
        cfile << setw(20) << setprecision(8) <<  Brute_MonteCarlo(pow(2,i),_t,std);
        cfile << setw(20) << setprecision(8) <<  _t<<endl;
   }
     cfile.close();
    _t = 0;
    i  = 3;

     ofstream dfile;
     string dfilename = "importance_sampling_monte_carlo.dat";
     dfile.open(dfilename);
     dfile << setiosflags(ios::showpoint | ios::uppercase);

     while (_t<t_max) {
        i++;
        dfile << setw(15) << setprecision(8) <<  pow(2,i);
        dfile << setw(20) << setprecision(8) <<  Importance_MonteCarlo(pow(2,i),_t, std);
        dfile << setw(20) << setprecision(8) <<  _t<<endl;
   }

     dfile.close();
     _t = 0;
     i  = 3;

     ofstream cpfile;
     string cpfilename = "parallel_brute_force_monte_carlo.dat";
     cpfile.open(cpfilename);
     cpfile << setiosflags(ios::showpoint | ios::uppercase);
     while (_t<t_max){
        i++;
        cpfile << setw(15) << setprecision(8) <<  pow(2,i);
        cpfile << setw(20) << setprecision(8) <<  Parallel_Brute_MonteCarlo(pow(2,i),_t,std);
        cpfile << setw(20) << setprecision(8) <<  _t<<endl;
   }
     cpfile.close();
    _t = 0;
    i  = 3;

     ofstream dpfile;
     string dpfilename = "parallel_importance_sampling_monte_carlo.dat";
     dpfile.open(dpfilename);
     dpfile << setiosflags(ios::showpoint | ios::uppercase);

     while (_t<t_max) {
        i++;
        dpfile << setw(15) << setprecision(8) <<  pow(2,i);
        dpfile << setw(20) << setprecision(8) <<  Parallel_Importance_MonteCarlo(pow(2,i),_t, std);
        dpfile << setw(20) << setprecision(8) <<  _t<<endl;
   }

     dpfile.close();


    //task3a();
    //task3b();
    //test_ans();
    //test_lag();
    double integralening ;
    double _std;
    Brute_MonteCarlo(10000000, -3,3, integralening,_std);
    cout <<integralening <<endl;
    return 0;
    }
