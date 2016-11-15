#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <sys/time.h>

#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#include <adolc/hypertensor/hyper_main.h>

//#include "./../function.h"
//#include "./../physics.hpp"
#include "./../synthetic.hpp"

#define def_tol 1.0e-12

#define tag 1

void generate_multi(int curr_level, int max_level, int d, int n,
                    int prev_value, int* multi, double** tensorhelp) {
  int i;
  if (curr_level > max_level) {
    int add;
    double ret;
    add = tensor_address(d, multi);
    ret = tensorhelp[0][add];
    if (fabs(ret) > def_tol) {
      std::cout << "T[ ";
      for(i=0; i<= max_level; i++) {
        std::cout << multi[i]-1 << " ";
      }
      std::cout << "] = " << ret << std::endl;
    }
  } else {
    for(i = prev_value; i >0 ; i--) {
      multi[curr_level] = i;
      generate_multi(curr_level+1, max_level, d, n, i, multi, tensorhelp);
      multi[curr_level] = 0;
    }
  }
}
int main(int argc, char *argv[]) {
  int n = NUM_IND;
  int i,j,k;
  adouble *xad;
  adouble fad;
  double f;
  double *x;
  x=new double[n];
  xad=new adouble[n];
  struct timeval tv1, tv2;
  double time_elapsed = 0;

  get_initials(x, n);

  printf("evaluating the function...");
  gettimeofday(&tv1, NULL);
  trace_on(tag);
  for(i=0;i<n;i++)
  {
    xad[i] <<= x[i];  
  }
// function evaluation
  fad = eval_func<adouble>(xad, n);

  fad >>= f;
  trace_off();
  printf("done!\n");
  printf(" y = %.6lf\n", f);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) + (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  double func_time = time_elapsed;


  int d = DERIVATIVE_ORDER;
  double** S = new double*[n];
  for (i = 0; i < n; ++i) {
    S[i] = new double[n];
    for (j = 0; j < n; j++) {
      S[i][j] = (i==j)? 1.0:0.0;
    }
  }
  int dim = binomi(n + d, d);
  std::cout << "d = " << d <<", dim = " << dim <<std::endl;
  double** tensorhelp = myalloc2(1, dim);
  tensor_eval(tag, 1, n, d, n, x, tensorhelp, S);

  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) + (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;

#define PRINT_RESULTS
#ifdef PRINT_RESULTS
  int* multi = new int[d];
  for (i=0; i<d; ++i) { multi[i] = 0;}
  for (i=0; i<d; i++) {
    generate_multi(0, i, d, n, n, multi, tensorhelp);
  }
#endif

  printf("eval_func time elapsed = %.10f\n", func_time);
  printf("total time elapsed = %.10f\n", time_elapsed);
/*
  for (i=0; i<d; ++i) { multi[i] = 0;}
  int add;
  double ret;
  for (i=0; i<n; ++i) {
    multi[0] = i + 1;
    add = tensor_address(d, multi);
    ret = tensorhelp[0][add];
    if (fabs(ret) > def_tol) {
      std::cout << "A[" << i <<  "] = " << ret << std::endl;
    }
  } 
  for (i=0; i<n; ++i) {
    multi[0] = i + 1;
    for(j=0; j<=i; ++j) {
      multi[1] = j + 1;
      add = tensor_address(d, multi);
      ret = tensorhelp[0][add];
      if (fabs(ret) > def_tol) {
        std::cout << "H[" << i << "," << j << "] = " << ret << std::endl;
      }
    }
  } 
  for (i=0; i<n; ++i) {
    multi[0] = i + 1;
    for(j=0; j<=i; ++j) {
      multi[1] = j + 1;
      for(k=0; k<=j; ++k) {
        multi[2] = k + 1;
        add = tensor_address(d, multi);
        ret = tensorhelp[0][add];
//      tensor_value(d, 1, &ret, tensorhelp, multi);
        if (fabs(ret) > def_tol) {
          std::cout << "T[" << i << ","
                    << j << "," << k << "] = " << ret << std::endl;
        }
      }
    }
  } 
*/

  delete[] x;
  delete[] xad;
  return 0;
}


