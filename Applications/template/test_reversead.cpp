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
#include "sys/time.h"

#include "func_eval.hpp"

#include "reversead/reversead.hpp"
using namespace ReverseAD;

//#define DUMMY_SCALAR
//#define PRINT_RESULT

//#define ORDER 2
//#define METHOD 2

double k_getTime() {
  struct timeval v;
  struct timezone z;
  gettimeofday(&v, &z);
  return ((double)v.tv_sec)+((double)v.tv_usec)/1000000;
}

/*
   options[0] : the order to be evaluated
   options[1] : 0 : specific, Adjoint, Hessian, Third
                1 : generic
                2 : flat code generator
   options[2] : 0/1, F/T, dump sparsity pattern
*/
double evaluate_derivatives(int n, int m, std::shared_ptr<TrivialTrace<double>> trace, int* options) {
  size_t order = options[0];
  double t1 = k_getTime();
  std::shared_ptr<DerivativeTensor<size_t, double>> tensor;
  if (options[1] == 0) { // specific
    if (order == 1) {
      BaseReverseAdjoint<double> adjoint(trace);
      tensor = adjoint.compute(n, m);
    } else if (order == 2) {
      BaseReverseHessian<double> hessian(trace);
      tensor = hessian.compute(n, m);
    } else if (order == 3) {
      BaseReverseThird<double> third(trace);
      tensor = third.compute(n, m);
    } else {
      fprintf(stderr, "Options[0] = %d : 1, 2, 3 when Options[1] == %d\n", options[0], options[1]); 
    }
  } else if (options[1] == 1) { // generic
    BaseReverseGeneric<double> reverse_tensor(trace, order);
    tensor = reverse_tensor.compute(n, m);
  } else if (options[1] == 2) { // flat code generator
    BaseReverseTensor<double> reverse_tensor(trace, order);
    tensor = reverse_tensor.compute(n, m);
  } else {
    fprintf(stderr, "Options[1] = %d : 0 (Specific), 1 (Generic), 2 (Flat Code)\n", options[1]);
  }
  double time_elapsed = k_getTime() - t1;
  size_t size;
  size_t** tind;
  double* values;
  tensor->get_internal_coordinate_list(0, order, &size, &tind, &values);

#ifdef PRINT_RESULT
  for (int i = 0; i < size; i++) {
    std::cout << "T[ ";
    for (int j = 0; j < order; j++) {
      std::cout << tind[i][j] << " ";
    }
    std::cout << "] = " << values[i] << std::endl;
  }
#endif

  printf("ReverseAD nnz[%u] method[%d] order[%d] timing = %.6f\n", size, options[1], options[0], time_elapsed);
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

  trace_on<double>();
  for (int i = 0; i < n; i++) {
    xad[i] <<= ind[i];
  }
  func_eval(n, xad, m, yad);
#ifdef DUMMY_SCALAR
  adouble zad = 0;
  double z = 0;
  for (int i = 0; i < m; i++) {
    zad = zad+yad[i]; 
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
  
  std::shared_ptr<TrivialTrace<double>> trace = trace_off<double>(); 
  tear_down();
  // Evaluate derivatives;
  int options[2] = {ORDER, METHOD};
  double t;
  std::cout << "m = " << m << "  n = " << n << std::endl;
  t = evaluate_derivatives(n, m, trace, options);

  delete[] xad;
  delete[] yad;

}


