#include <cmath>
#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include "mesh.h"
#include "mesh_eval.hpp"

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
double ep2 = 1.0e-10;

double k_getTime() {
   struct timeval v;
   struct timezone z;
   gettimeofday(&v, &z);
   return ((double)v.tv_sec)+((double)v.tv_usec)/1000000;
}
#ifdef USING_ADOLC
void generate_multi(int curr_level, int max_level, int d, int n,
                    int prev_value, int* multi, double** tensorhelp) {
  int i;
  if (curr_level > max_level) {
    int add;
    double ret;
    add = tensor_address(d, multi);
    ret = tensorhelp[0][add];
    if (fabs(ret) > ep2) {
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
#endif // USING_ADOLC

void init_dim(int *n, int *m);
void init_startpoint(double *x, int n);

int main(int argc, char* argv[])
{
  Mesh *m = nullptr;
  std::cout << "reading mesh ... " << meshFile << " ";
  if (readMesh(meshFile, &m)) {
    freeMesh(&m);
    return -1;
  }
  std::cout << "done " << std::endl;
  int i, j, k, l;
  int n = 3 * m->nv;

  double f;
  double *x, *c;

  double t1, t2;

  std::cout << "n = " << n << std::endl;
  x = new double[n];

  for (int i = 0; i < 3 * m->nv; i++) {
    x[i] = m->v[i];
  }  

  t1 = k_getTime();
  f = aFcn(x, m->ne, m->e);
  t2 = k_getTime();
  std::cout << "Function value : " << f << std::endl;
  std::cout << "The time needed for scalar function evaluation : " << t2 - t1 << std::endl;

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
  fad = aFcn<adouble>(xad, m->ne, m->e);
  fad >>= f;
  t2 = k_getTime();
#endif  

#ifdef USING_ADTOOL
  std::cout << "Overloaded function value : " << f << std::endl;
  std::cout << "The time needed for overloaded function evaluation : " << t2 - t1 << std::endl;
#ifdef USING_REVERSEAD
  std::shared_ptr<TrivialTrace<double>> trace = trace_off<double>();
  int order = 0;
#ifdef COMPUTE_GRADIENT
  order = 1;
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
  order = 2;
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
#ifdef COMPUTE_THIRD
  order = 3;
  t1 = k_getTime();
  BaseReverseThird<double> third(trace);
  std::shared_ptr<DerivativeTensor<size_t, double>> tensor = third.compute(n, 1);
  t2 = k_getTime();
  std::cout << "ReverseAD Third cost : " << t2 - t1 << std::endl;
  size_t size;
  size_t** tind;
  double* values;
  tensor->get_internal_coordinate_list(0, 3, &size, &tind, &values);
  std::cout << "size of third = " << size << std::endl;
#endif // COMPUTE_THIRD
#ifdef COMPUTE_HIGHER_ORDER
  order = atoi(argv[1]);
  t1 = k_getTime();
  BaseReverseTensor<double> reverse_tensor(trace, order);
  std::shared_ptr<DerivativeTensor<size_t, double>> tensor = reverse_tensor.compute(n, 1);
  t2 = k_getTime();
  std::cout << "ReverseAD Order["<<order<<"] cost : " << t2 - t1 << std::endl;
  size_t size;
  size_t** tind;
  double* values;
  tensor->get_internal_coordinate_list(0, order, &size, &tind, &values);
  std::cout << "size of Order["<<order<<"] = " << size << std::endl;
#ifdef CHECK_HIGHER_ORDER
  for (int i = 0; i < size; i++) {
    std::cout << "T[ ";
    for (int j = 0; j < order; j++) {
      std::cout << tind[i][j] << " "; 
    }
    std::cout << "] = " << values[i] << std::endl;
  }
#endif // CHECK_HIGHER_ORDER
#endif // COMPUTE_HIGHER_ORDER
if (argc > 2) { // the second argument
  std::ofstream of("sparsity.pat");
  of << order << " " << size << std::endl;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < order; j++) {
     of << tind[i][j] << " ";
    }
    of << std::endl;
  }
  of.close();
}
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
  if (argc <= 2) {
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
  } else { // argc > 2
    int num_rows = atoi(argv[2]);
/*
    double** v = new double*[n];
    double** hv = new double*[n];
    for (int i = 0; i < n; i++) {
      v[i] = new double[num_rows];
      hv[i] = new double[num_rows];
      for (int j = 0; j < num_rows; j++) {
        hv[i][j] = 0.0;
        v[i][j] = 0.0;
      }
    }
    for (int i = 0; i < num_rows; i++) {
      v[i][i] = 1.0;
    }
    t1 = k_getTime();
    hess_mat(1, n, num_rows, x, v, hv);
    t2 = k_getTime();
    std::cout << "ADOLC Hessian-Mat cost : " << t2 - t1 << std::endl;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < num_rows; j++) {
        hessian_ad[i][j] = hv[i][j];
      }
    }
*/
    double* v = new double[n];
    double* hv = new double[n];
    for (int i = 0; i < n; i++) {
      v[i] = 0.0; hv[i] = 0.0;
    }
    double dummy_sum = 0.0;
    t1 = k_getTime();
    for (int i = 0; i < num_rows; i++) {
      v[i] = 1.0;
      hess_vec(1, n, x, v, hv);
      v[i] = 0.0;
      for (int j = 0; j < n; j++) {
        hessian_ad[j][i] = hv[j];
        dummy_sum += hv[j];
      }
    }
    t2 = k_getTime();
    std::cout << dummy_sum << " ADOLC Hessian-Vec cost : " << t2 - t1 << std::endl;
  }
#endif // COMPUTE_HESSIAN
#ifdef COMPUTE_HIGHER_ORDER
  int order = atoi(argv[1]);
  double** seed = new double*[n];
  for (int i = 0; i < n; i++) {
    seed[i] = new double[n];
    for (int j = 0; j < n; j++) {
      seed[i][j] = ((i==j)?1.0:0.0);
    }
  }
  int dim = binomi(n+order, order);
  double** tensorhelp = myalloc2(1, dim);
  std::cout << "order = " << order << ", dim = " << dim << std::endl;
  t1 = k_getTime();
  tensor_eval(1, 1, n, order, n, x, tensorhelp, seed);
  t2 = k_getTime();
  std::cout << "ADOLC Order["<<order<<"] cost : " << t2 - t1 << std::endl;
#ifdef CHECK_HIGHER_ORDER
  int multi[order];
  for (int i = 0; i < order; i++) {multi[i] = 0;}
  for (int i = 0; i < order; i++) {
    generate_multi(0, i, order, n, n, multi, tensorhelp);
  }
#endif // CHECK_HIGHER_ORDER
 
#endif // COMPUTE_HIGHER_ORDER
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
    fci = aFcn<double>(x, m->ne, m->e);
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
        fci = aFcn<double>(x,m->ne, m->e);
        x[j] += ep;
        fcij = aFcn<double>(x,m->ne, m->e);
        x[i] -= ep;
        fcj = feval<double>(x,m->ne, m->e);
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

