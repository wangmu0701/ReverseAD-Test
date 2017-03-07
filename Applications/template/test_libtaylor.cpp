/*
  This should be a unified test file for some problems, with two rules:
  1. ALL application data are passed as arguments to the program
  2. ALL method data are controled by macros.
  NOTE : The compile may take a long time. Deal with it.

  INTERFACE:
 
  void set_up(argc, argv, vector<double>& ind, &n, &m); // number of independent and dependent variables

  void func_eval<T>(n, xad, m, yad);

  void tear_down();   // may not be necessary
*/
#include <memory>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "sys/time.h"

#include "func_eval.hpp"

#include "ctaylor.hpp"

//#define DUMMY_SCALAR
//#define PRINT_RESULT

//#define ORDER 2

// 0 : W/O sparsity, 1: need sparsity.txt 
//#define METHOD 0

double k_getTime() {
  struct timeval v;
  struct timezone z;
  gettimeofday(&v, &z);
  return ((double)v.tv_sec)+((double)v.tv_usec)/1000000;
}


double evaluate_derivatives(int n, int m, double* x) {
  typedef ctaylor<double, ORDER> ttype;
  ttype* xad = new ttype[n];
  ttype* yad = new ttype[m];
  for (int i = 0; i < n; i++) {
    xad[i] = x[i];
  }

  double t1 = k_getTime();
  if (METHOD == 0) {
    if (ORDER == 1) {
      for (int i = 0; i < n; i++) {
        xad[i].set(VAR0, 1);
        func_eval<ttype>(n, xad, m, yad);
        xad[i].set(VAR0, 0);
      }
    } else if (ORDER == 2) {
      for (int i = 0; i < n; i++) {
        xad[i].set(VAR0, 1);
        for (int j = 0; j <= i; j++) {
          xad[j].set(VAR1, 1);
          func_eval<ttype>(n, xad, m, yad);
          xad[j].set(VAR1, 0);
        }
        xad[i].set(VAR0, 0);
      }
    } else if (ORDER == 3) {
      for (int i = 0; i < n; i++) {
        xad[i].set(VAR0, 1);
        for (int j = 0; j <= i; j++) {
          xad[j].set(VAR1, 1);
          for (int k = 0; j <= k; k++) {
            xad[k].set(VAR2, 1);
            func_eval<ttype>(n, xad, m, yad);
            xad[k].set(VAR2, 0);
          }
          xad[j].set(VAR1, 0);
        }
        xad[i].set(VAR0, 0);
      }
    } else {
      std::cout << "Up to third order for libtaylor" << std::endl;
    }
  } else {
    std::ifstream inf("sparsity.txt", std::ifstream::in);
    size_t order;
    size_t size;
    inf >> order >> size;
    assert(order == ORDER);
    std::cout << "SpTaylor with " << size << " nonzeros." << std::endl;
    if (ORDER == 1) {
    } else if (ORDER == 2) {
      size_t tind[2];
      for (int i = 0; i < size; i++) {
        inf >> tind[0] >> tind[1];
        xad[tind[0]].set(VAR0, 1);
        xad[tind[1]].set(VAR1, 1);
        func_eval<ttype>(n, xad, m, yad);
        xad[tind[0]] = x[tind[0]];
        xad[tind[1]] = x[tind[1]];
      }
    } else if (ORDER == 3) {
      size_t tind[3];
      for (int i = 0; i < size; i++) {
        inf >> tind[0] >> tind[1] >> tind[2];
        xad[tind[0]].set(VAR0, 1);
        xad[tind[1]].set(VAR1, 1);
        xad[tind[2]].set(VAR2, 1);
        func_eval<ttype>(n, xad, m, yad);
        xad[tind[0]] = x[tind[0]];
        xad[tind[1]] = x[tind[1]];
        xad[tind[2]] = x[tind[2]];
      }
    }
  }

  double time_elapsed = k_getTime() - t1;
  printf("LibTaylor method[%d] order[%d] timing = %.6f\n", METHOD, ORDER, time_elapsed);
  return time_elapsed;
}

int main(int argc, char* argv[]) {
  vector<double> ind;
  int n = 0;
  int m = 0;
  set_up(argc, argv, ind, n, m);
  std::cout << "num_dep = " << m << ", num_ind = " << n << std::endl;

  // Evaluate derivatives;
  double t = evaluate_derivatives(n, m, ind.data());
  tear_down();

}


