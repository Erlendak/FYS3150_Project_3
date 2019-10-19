#include <iostream>
#include <integrand.h>
#include <methods.h>
#include <cmath>
#include <armadillo>
#include <monte_carlo.h>
#include <unit_tests.h>
#include <omp.h> //The Open MPI package


using namespace std;
using namespace arma;

    int main(){
        int nthreads = omp_get_num_threads();
        printf("Using %d threads outside parallel loop\n",nthreads);
        int mthreads = omp_get_max_threads();
        printf("There are %d threads available at a time\n",mthreads);
        double std;
        int M = 9;
        double _t ;

        int _N = 25;

/*
    ofstream afile;
    string afilename = "gauss_legendre.dat";
    afile.open(afilename);
    afile << setiosflags(ios::showpoint | ios::uppercase);
      for (int i=2;i<= _N;i++){

         afile << setw(15) << setprecision(8) <<  i;
         afile << setw(15) << setprecision(8) <<  task3a(i,_t);
         afile << setw(15) << setprecision(8) <<  _t<<endl;
    }
     afile.close();

     ofstream bfile;
     string bfilename = "gauss_laguerre.dat";
     bfile.open(bfilename);
     bfile << setiosflags(ios::showpoint | ios::uppercase);
     for (int i=2;i<= _N;i++){

        bfile << setw(15) << setprecision(8) <<  i;
        bfile << setw(15) << setprecision(8) <<  task3b(i,_t);
        bfile << setw(15) << setprecision(8) <<  _t<<endl;
   }

     bfile.close();
*/
     ofstream cfile;
     string cfilename = "brute_force_monte_carlo.dat";
     cfile.open(cfilename);
     cfile << setiosflags(ios::showpoint | ios::uppercase);
     for (int i=0;i<= M;i++){

        cfile << setw(15) << setprecision(8) <<  pow(10,i);
        cfile << setw(15) << setprecision(8) <<  Brute_MonteCarlo(pow(10,i),_t,std);
        cfile << setw(15) << setprecision(8) <<  _t<<endl;
   }


     cfile.close();
/*
     ofstream dfile;
     string dfilename = "importance_sampling_monte_carlo.dat";
     dfile.open(dfilename);
     dfile << setiosflags(ios::showpoint | ios::uppercase);

     for (int i=0;i<= M;i++){

        dfile << setw(15) << setprecision(8) <<  pow(10,i);
        dfile << setw(15) << setprecision(8) <<  Importance_MonteCarlo(pow(10,i),_t, std);
        dfile << setw(15) << setprecision(8) <<  _t<<endl;
   }

     dfile.close();
/*
     ofstream fdfile;
     string fdfilename = "importance_sampling_monte_carlo.dat";
     fdfile.open(fdfilename);
     fdfile << setiosflags(ios::showpoint | ios::uppercase);
     for (int i=2;i<= _N;i++){

        afile << setw(15) << setprecision(8) <<  i;
        afile << setw(15) << setprecision(8) <<  task3a(i,_t);
        afile << setw(15) << setprecision(8) <<  _t<<endl;
   }
     fdfile.close();

*/
    return 0;
    }
