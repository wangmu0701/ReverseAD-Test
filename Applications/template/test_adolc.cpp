/*
  This should be a unified test file for some problems, with two rules:
  1. ALL application data are passed as arguments to the program
  2. ALL method data are controled by macros.
  NOTE : The compile may take a long time. Deal with it.

  INTERFACE:
 
  void set_up(argc, argv, vector<double>& ind, &n, &m); // number of independent and dependent variables

  void func_eval(n, xad, m, yad);

  void tear_down();   // may not be necessary
*/
#include <memory>
#include <cstdlib>
#include <cassert>
#include "sys/time.h"

#include "func_eval.hpp"

#include "adolc/adolc.h"
#include "adolc/adolc_sparse.h"

#define TAG 1

//#define DUMMY_SCALAR
//#define PRINT_RESULT

//#define ORDER 2
//#define METHOD 4

double k_getTime() {
  struct timeval v;
  struct timezone z;
  gettimeofday(&v, &z);
  return ((double)v.tv_sec)+((double)v.tv_usec)/1000000;
}

/*
   options[0] : the order to be evaluated
   options[1] : 0 : Tensor
                1 : Direct
                2 : Indirect
                3 : FullHess
                4 : Single HV
                5 : Dense SecRev
                6 : Sparse SecRev
                7 : Hessian-Matrix                
*/

double evaluate_derivatives(int n, int m, double* x, int* options) {
  int order = options[0];
  int nnz;
  double t1 = k_getTime();
  if (options[1] == 0) { // Teed = new double*[n];
    assert(m == 1);
    double** seed = new double*[n];
    for (int i = 0; i < n; i++) {
      seed[i] = new double[n];
      for (int j = 0; j < n; j++) {
        seed[i][j] = ((i==j)?1.0:0.0);
      }
    }
    int dim = binomi(n+order, order);
    double** tensorhelp = myalloc2(1, dim);
    tensor_eval(TAG, 1, n, order, n, x, tensorhelp, seed);
    for (int i = 0; i < n; i++) {
      delete[] seed[i];
    }
    delete[] seed; 
    myfree2(tensorhelp);
  } else {
    if (order == 2)  { // Hessian
      assert(m == 1);
      if (options[1] == 1 || options[1] == 2) { // Direct or Indirect
        int opt[2] = {0, 0}; // default is indirect;
        if (options[1] == 1) {opt[0] = 1;} // set direct;
        unsigned int * rind = NULL;
        unsigned int * cind = NULL;
        double * values = NULL;
        sparse_hess(TAG, n, 0, x, &nnz, &rind, &cind, &values, opt);
#ifdef PRINT_RESULT
        for (int i = 0; i < nnz; i++) {
          printf("H[%d, %d] = %.6f\n", rind[i], cind[i], values[i]);
        }
#endif
        free(rind);
        free(cind);
        free(values);
      } else if (options[1] == 3) { // FullHess
        double** H = new double*[n];
        for (int i = 0; i < n; i++) {
          H[i] = new double[n];
        }
        hessian(TAG, n, x, H);
        nnz = n*n;
#ifdef PRINT_RESULT
        for (int i = 0; i < n; i++) {
          for (int j = 0; j <= i; j++) {
            printf("H[%d, %d] = %.6f\n", i, j, H[i][j]);
          }
        }
#endif
        for (int i = 0; i < n; i++) {
          delete[] H[i];
        }
        delete[] H;
      } else if (options[1] == 4) { // Single Hv
        double v[n];
        double Hv[n];
        for (int i = 0; i < n; i++) {
          v[i] = 1.0;
          Hv[i] = 0.0;
        }
        hess_vec(TAG, n, x, v, Hv);
        nnz = n;
      } else if (options[1] == 5) { // dense second order reverse
        double** H = new double*[n];
        for (int i = 0; i < n; i++) {
          H[i] = new double[n];
        }
        hessian_dense(TAG, n, x, H);
        nnz = n*n;
#ifdef PRINT_RESULT
        for (int i = 0; i < n; i++) {
          for (int j = 0; j <= i; j++) {
            printf("H[%d, %d] = %.6f\n", i, j, H[i][j]);
          }
        }
#endif
        for (int i = 0; i < n; i++) {
          delete[] H[i];
        }
        delete[] H;
      } else if (options[1] == 6){ // sparse second order reverse
        unsigned int * rind = NULL;
        unsigned int * cind = NULL;
        double * values = NULL;
        hessian_sparse(TAG, n, x, &nnz, &rind, &cind, &values);
#ifdef PRINT_RESULT
        for (int i = 0; i < nnz; i++) {
          printf("H[%d, %d] = %.6f\n", rind[i], cind[i], values[i]);
        }
#endif
        free(rind);
        free(cind);
        free(values);
      } else if (options[1] == 7) { // Hess-matrix options
        double** H  = myalloc2(n, n);
        double y;
        double*** Xppp = myalloc3(n, n, 1);
        double*** Yppp = myalloc3(1, n, 1);
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
             Xppp[i][j][0] = 0;
          }
          Xppp[i][i][0] = 1.0;
        }
        double** Upp = myalloc2(1,2);
        Upp[0][0] = 1; Upp[0][1] = 0;
        double*** Zppp = myalloc3(n, n, 2);
        int ret_val = hov_wk_forward(TAG,1,n,1,2,n,x,Xppp,&y,Yppp);
        ret_val = hos_ov_reverse(TAG,1,n,1,n,Upp,Zppp);
        for (int i = 0; i < n; ++i) {
          for (int l = 0; l < n; ++l) {
            H[l][i] = Zppp[i][l][1];
          }
        }
#ifdef PRINT_RESULT
        for (int i = 0; i < n; i++) {
          for (int j = 0; j <= i; j++) {
            printf("H[%d, %d] = %.6f\n", i, j, H[i][j]);
          }
        }
#endif
        myfree2(H);
        myfree3(Xppp);
        myfree3(Yppp);
        myfree2(Upp);
        myfree3(Zppp);
      }
    } else if (order == 1) { // Gradient or Jacobian
      if (m == 1) { // gradient
        double g[n];
        gradient(TAG, n, x, g);
#ifdef PRINT_RESULT
        for (int i = 0; i < n; i++) {
          printf("g[%d] = %.6f\n", i, g[i]);
        }
#endif
      } else { // jacobian
        double** J = new double*[m];
        for (int i = 0; i < m; i++) {
          J[i] = new double[n];
        }
        jacobian(TAG, m, n, x, J);
#ifdef PRINT_RESULT
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < n; j++) {
            printf("J[%d][%d] = %.6f\n", i, j, J[i][j]);
          }
        }
#endif
        for (int i = 0; i < m; i++) {
          delete[] J[i];
        }
        delete[] J;
      }
      nnz = n*m;
    }
  }


  double time_elapsed = k_getTime() - t1;
  size_t size;
  size_t** tind;
  double* values;
  printf("ADOLC nnz[%d] method[%d] order[%d] timing = %.6f\n", nnz, options[1], options[0], time_elapsed);
  return time_elapsed;
}

int main(int argc, char* argv[]) {
  vector<double> ind;
  int n = 0;
  int m = 0;
  set_up(argc, argv, ind, n, m);
  std::cout << "num_dep = " << m << ", num_ind = " << n << std::endl;
  adouble* xad = new adouble[n];
  adouble* yad = new adouble[m];
  double* y = new double[m];
  double* x = new double[n];

  for (int i = 0; i < n; i++) {
    x[i] = ind[i];
  }
  trace_on(TAG);
  for (int i = 0; i < n; i++) {
    xad[i] <<= x[i];
  }
  func_eval<adouble>(n, xad, m, yad);
#ifdef DUMMY_SCALAR
  adouble zad = 0;
  double z;
  for (int i = 0; i < m; i++) {
    zad = zad + yad[i];
  }
  zad >>= z;
  std::cout << "func_eval = " << z << std::endl;
  m = 1;
#else
  for (int i = 0; i < m; i++) {
    yad[i] >>= y[i];
    //std::cout << "y["<<i<<"] = " << y[i] << std::endl;
  }
#endif
  trace_off(); 

  tear_down();
  // Evaluate derivatives;
  int options[2] = {ORDER, METHOD};
  double t;
  std::cout << "m = " << m << " n = " << n << std::endl;
  t = evaluate_derivatives(n, m, x, options);

  delete[] xad;
  delete[] yad;

}


