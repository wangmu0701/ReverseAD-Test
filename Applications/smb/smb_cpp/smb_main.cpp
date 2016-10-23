#include <cmath>
#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include "eval_func_chem.hpp"

using namespace std;

#ifdef USING_REVERSEAD
#include "reversead/reversead.hpp"
using namespace ReverseAD;
#endif
#ifdef USING_ADOLC
#include "adolc/adolc.h"
#include "adolc/adolc_sparse.h"
#endif

double ep = 1.0e-5;

double k_getTime() {
   struct timeval v;
   struct timezone z;
   gettimeofday(&v, &z);
   return ((double)v.tv_sec)+((double)v.tv_usec)/1000000;
}

void init_dim(int *n, int *m);
void init_startpoint(double *x, int n);

int main(int argc, char* argv[])
{
  int i, j, k, l;
  constexpr int n = get_n();
  constexpr int m = get_m();
  double f;
  double *x, *c;

  double t1, t2;

  //init_dim(&n,&m); // initialize n and m

  std::cout << "n = " << n << ", m = " << m << std::endl;
  x = new double[n];
  c = new double[m];

  init_startpoint(x,n);

  t1 = k_getTime();
  f = feval(x,n);
  t2 = k_getTime();
  std::cout << "Function value : " << f << std::endl;
  std::cout << "The time needed for scalar function evaluation : " << t2 - t1 << std::endl;
/*
  t1 = k_getTime();
  for(i=0;i<repnum;i++) {
    ceval(x,c,n);
  }
  t2 = k_getTime();
  std::cout << "The time needed for vector constraint evaluation : " << t2 - t1 << std::endl;
*/
// Active Section
#ifdef USING_ADTOOL
  adouble fad, *xad;
  xad = new adouble[n];
  t1 = k_getTime();
#ifdef USING_REVERSEAD
  trace_on<double>();
#endif
#ifdef USING_ADOLC
  trace_on(1);
#endif
  for (int i =0 ; i < n; i++) {
    xad[i] <<= x[i];
  }
  fad = feval<adouble>(xad, n);
  fad >>= f;
  t2 = k_getTime();
#endif  

#ifdef USING_ADTOOL
  std::cout << "Overloaded function value : " << f << std::endl;
  std::cout << "The time needed for overloaded function evaluation : " << t2 - t1 << std::endl;
#ifdef USING_REVERSEAD
  std::shared_ptr<TrivialTrace<double>> trace = trace_off<double>();
#ifdef COMPUTE_GRADIENT
  double gradient_ad[n];
  t1 = k_getTime();
  BaseReverseAdjoint<double> adjoint(trace);
  std::shared_ptr<DerivativeTensor<size_t, double>> tensor = adjoint.compute(n, 1);
  t2 = k_getTime();
  std::cout << "ReverseAD Gradient cost : " << t2 - t1 << std::endl;
  size_t size;
  size_t** tind;
  double* values;
  tensor->get_internal_coordinate_list(0, 1, &size, &tind, &values);
  for (int i = 0; i < n; i++) { gradient_ad[i] = 0.0;}
  for (int i = 0; i < size; i++) {
    //cout << "A["<<tind[i][0]<<"] = " << values[i] << endl;
    if (tind[i][0] < n) {
      gradient_ad[tind[i][0]] = values[i];
    }
  }
#endif // REVERSEAD_GRADIENT 
#ifdef COMPUTE_HESSIAN
  double hessian_ad[n][n];
  t1 = k_getTime();
  BaseReverseHessian<double> hessian(trace);
  std::shared_ptr<DerivativeTensor<size_t, double>> tensor = hessian.compute(n, 1);
  t2 = k_getTime();
  std::cout << "ReverseAD Hessian cost : " << t2 - t1 << std::endl;
  size_t size;
  size_t** tind;
  double* values;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      hessian_ad[i][j] = 0;
    }
  }
  tensor->get_internal_coordinate_list(0, 2, &size, &tind, &values);
  std::cout << "size of hessian = " << size << std::endl;
  for (int i = 0; i < size; i++) {
    //cout << "H["<<tind[i][0]<<","<<tind[i][1]<<"] = " << values[i]<<endl;
    hessian_ad[tind[i][0]][tind[i][1]] = values[i];
  }
#endif // COMPUTE_HESSIAN
#endif // USING_REVERSEAD
#ifdef USING_ADOLC
  trace_off();
#ifdef COMPUTE_GRADIENT
  double gradient_ad[n]; 
  t1 = k_getTime();
  gradient(1, n, x, gradient_ad);
  t2 = k_getTime();
  std::cout << "ADOLC Gradient cost : " << t2 - t1 << std::endl;
#endif // COMPUTE_GRADIENT
#ifdef COMPUTE_HESSIAN
  double hessian_ad[n][n];
  for (int i =0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      hessian_ad[i][j] = 0.0;
    }
  }
  int options[2] = {0, 1}; 
  if (argc > 1) {
    options[1] = atoi(argv[1]);
  }
  unsigned int * rind = NULL;
  unsigned int * cind = NULL;
  double * values = NULL;
  int nnz;
  t1 = k_getTime();
  sparse_hess(1, n, 0, x, &nnz, &rind, &cind, &values, options);
  t2 = k_getTime();  
  if (options[1] == 0) {
    std::cout << "(direct)";
  } else if (options[1] == 1) {
    std::cout << "(indirect)";
  }
  std::cout << "ADOLC Hessian cost : " << t2 - t1 << std::endl;
  std::cout << "size of hessian = " << nnz << std::endl;
  for (int i = 0; i < nnz; i++) {
    //cout << "H["<<rind[i]<<","<<cind[i]<<"] = "<<values[i]<< std::endl;
    hessian_ad[cind[i]][rind[i]] = values[i];
  } 
#endif
#endif // USING_ADOLC
#endif


#ifdef CHECK_GRADIENT
  double fci, fdg;
  std::cout<<setfill('-')<<setw(10)<<"index"<<setw(20)<<"FiniteDiff";
#ifdef COMPUTE_GRADIENT
  std::cout<<setw(20)<<"ReverseAD"<<setw(20)<<"(Abs/Rel)Error";
#endif
  std::cout<<std::endl<<setfill(' ')<<setprecision(12);
  for (int i = 0; i < n; i++) {
    init_startpoint(x, n);
    x[i] += ep;
    fci = feval<double>(x, n);
    fdg = (fci - f) / ep;
    std::cout<<setw(10)<<i<<setw(20)<<fdg;
#ifdef COMPUTE_GRADIENT
    std::cout<<setw(20)<<gradient_ad[i];
    double e = 0.0;
    if (fabs(fdg+gradient_ad[i]) > ep) {
      e = (fabs(fdg-gradient_ad[i])/fabs(fdg+gradient_ad[i]));
    } else {
      e = fabs(fdg - gradient_ad[i]);
    }
    std::cout<<setw(20)<<e;
#endif
    std::cout << std::endl;
  }
#endif
#ifdef CHECK_HESSIAN
  double fci, fcj, fcij;
  double fdh;
  cout<<setfill('-')<<setw(5)<<"i"<<setw(5)<<"j"<<setw(20)<<"FiniteDiff";
#ifdef COMPUTE_HESSIAN
  cout<<setw(20)<<"ReverseAD"<<setw(20)<<"(Abs/Rel)Error";
#endif
  cout<<endl<<setfill(' ')<<setprecision(12);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <=i; j++) {
#ifdef COMPUTE_HESSIAN
      if (fabs(hessian_ad[i][j]) != 0.0) {
#endif
        init_startpoint(x, n);
        x[i] += ep;
        fci = feval<double>(x,n);
        x[j] += ep;
        fcij = feval<double>(x,n);
        x[i] -= ep;
        fcj = feval<double>(x,n);
        fdh = (fcij-fci-fcj+f) / (ep * ep);
        cout<<setw(5)<<i<<setw(5)<<j<<setw(20)<<fdh;
#ifdef COMPUTE_HESSIAN
        cout<<setw(20)<<hessian_ad[i][j];
        double e = 0.0;
        if (fabs(fdh+hessian_ad[i][j]) > ep) {
          e = (fabs(fdh-hessian_ad[i][j])/fabs(fdh+hessian_ad[i][j]));
        } else {
          e = fabs(fdh - hessian_ad[i][j]);
        }
        cout<<setw(20)<<e;
#endif
        cout << endl;
#ifdef COMPUTE_HESSIAN
      }
#endif
    }
  }
#endif 
  return 0;

}

