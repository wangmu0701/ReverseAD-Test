#include <memory>
#include <iostream>
#include <sys/time.h>

#include "reversead/reversead.hpp"

//#include "./../function.h"
//#include "./../physics.hpp"
#include "./../synthetic.hpp"


using ReverseAD::adouble;
using ReverseAD::TrivialTrace;
using ReverseAD::BaseReverseHessian;
using ReverseAD::BaseReverseThird;
using ReverseAD::BaseReverseGeneric;
using ReverseAD::BaseReverseGenerator;
using ReverseAD::BaseReverseTensor;
using ReverseAD::DerivativeTensor;
using ReverseAD::trace_on;
using ReverseAD::trace_off;

int main() {
  int n = NUM_IND;
  adouble *xad;
  adouble yad;
  double y;
  double *x;
  x=new double[n];
  xad=new adouble[n];
  struct timeval tv1, tv2;
  double time_elapsed = 0;

  get_initials(x, n);

  printf("evaluating the function...");
  gettimeofday(&tv1, NULL);

  trace_on<double>(); // begin tracing
  for (int i = 0; i < n; i++) {
    xad[i] <<= x[i];
  }
  yad = eval_func<adouble>(xad, n);
  yad >>= y;
  std::shared_ptr<TrivialTrace<double>> trace =
      trace_off<double>(); // end tracing
  //trace->dump_trace();
  printf("done!\n");
  std::cout << "y = " << y << std::endl;
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) + (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  double func_time = time_elapsed;

  BaseReverseGeneric<double> generic(trace, DERIVATIVE_ORDER);
//  BaseReverseThird<double> generic(trace);
  std::shared_ptr<DerivativeTensor<int, double>> tensor =
      generic.compute(NUM_IND, 1).get_tensor();

  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) + (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;

#ifdef PRINT_RESULTS
  // retrieve results
  int size;
  int** tind;
  double* values;
  // adjoints : dep[0].order[1]
  for (int order = 1; order <=DERIVATIVE_ORDER; order++) {
    tensor->get_internal_coordinate_list(0, order, &size, &tind, &values);
    for (int i = 0; i < size; i++) {
      std::cout << "T[ " << tind[i][0];
      for (int j = 1; j < order; j++) {
        std::cout << " " << tind[i][j]; 
      }
      std::cout << " ] = " << values[i] << std::endl;
    }
  }
#endif
/*
  std::cout << "size of adjoints = " << size << std::endl;
  for (int i = 0; i < size; i++) {
    std::cout << "A["<< tind[i][0] << "] = " << values[i] << std::endl;
  }
  // hessian : dep[0].order[2]
  tensor->get_internal_coordinate_list(0, 2, &size, &tind, &values);
  std::cout << "size of hessian = " << size << std::endl;
  for (int i = 0; i < size; i++) {
    std::cout << "H["<< tind[i][0] << ", " << tind[i][1]
              << "] = " << values[i] << std::endl;
  }
*/
  printf("eval_func time elapsed = %.10f\n", func_time);
  printf("total time elapsed = %.10f\n", time_elapsed);
}
