#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/time.h>

#include "mesh.h"
#include "mesh_eval.hpp"

using namespace std;

#include "ctaylor.hpp"

double ep = 1.0e-5;
double ep2 = 1.0e-10;

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

#ifdef COMPUTE_GRADIENT
  t1 = k_getTime();
  typedef ctaylor<double, 1> ttype;
  ttype* xad = new ttype[n];
  ttype yad = 0;
  double gradient_ad[n];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      xad[j] = x[j];
    } 
    xad[i].set(VAR0, 1);
    yad = aFcn<ttype>(xad, m->ne, m->e);
    gradient_ad[i] = yad.get(VAR0);
  }
  t2 = k_getTime();
  std::cout << "LibTaylor Graident cost : " << t2 - t1 << std::endl;
#endif // COMPUTE_GRADIENT

#ifdef COMPUTE_HESSIAN
  typedef ctaylor<double, 2> ttype;
  ttype* xad = new ttype[n];
  ttype yad = 0;
  double hessian_ad[n][n];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      hessian_ad[i][j] = 0.0;
    }
  }
  for (int i = 0; i < n; i++) {
    xad[i] = x[i];
  }
  t1 = k_getTime();
if (argc >= 2) {
  std::ifstream inf("sparsity.pat", std::ifstream::in);
  size_t order;
  size_t size;
  size_t tind[2];
  inf >> order >> size;  
  std::cout << "size of Hessian = " << size << std::endl;
  assert(order == 2);
  for (int i = 0; i < size; i++) {
    inf >> tind[0] >> tind[1];
    //std::cout << tind[0] << ", " << tind[1] << std::endl;
    xad[tind[0]].set(VAR0, 1);
    xad[tind[1]].set(VAR1, 1);
    yad = aFcn<ttype>(xad, m->ne, m->e);
    hessian_ad[tind[0]][tind[1]] = yad.get(VAR0|VAR1);
    xad[tind[0]] = x[tind[0]];
    xad[tind[1]] = x[tind[1]];
  }
} else {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      xad[i].set(VAR0, 1);
      xad[j].set(VAR1, 1);
      yad = aFcn<ttype>(xad, m->ne, m->e);
      hessian_ad[i][j] = yad.get(VAR0|VAR1);
      xad[i] = x[i];
      xad[j] = x[j];
    }
  }
}
  t2 = k_getTime();
  std::cout << "LibTaylor Hessian cost : " << t2 - t1 << std::endl;
#endif // COMPUTE_HESSIAN

#ifdef COMPUTE_THIRD
  typedef ctaylor<double, 3> ttype;
  ttype* xad = new ttype[n];
  ttype yad = 0;
  double third_ad[n][n][n];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
        third_ad[i][j][k] = 0.0;
      }
    }
  }
  for (int i = 0; i < n; i++) {
    xad[i] = x[i];
  }
  t1 = k_getTime();
if (argc >= 2) {
  std::ifstream inf("sparsity.pat", std::ifstream::in);
  size_t order;
  size_t size;
  size_t tind[3];
  inf >> order >> size;  
  std::cout << "size of Third = " << size << std::endl;
  assert(order == 3);
  for (int i = 0; i < size; i++) {
    inf >> tind[0] >> tind[1] >> tind[2];
    //std::cout << tind[0] << ", " << tind[1] << std::endl;
    xad[tind[0]].set(VAR0, 1);
    xad[tind[1]].set(VAR1, 1);
    xad[tind[2]].set(VAR2, 1);
    yad = aFcn<ttype>(xad, m->ne, m->e);
    third_ad[tind[0]][tind[1]][tind[2]] = yad.get(VAR0|VAR1|VAR2);
    xad[tind[0]] = x[tind[0]];
    xad[tind[1]] = x[tind[1]];
    xad[tind[2]] = x[tind[2]];
  }
} else {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      for (int k = 0; k <= j; k++) {
        xad[i].set(VAR0, 1);
        xad[j].set(VAR1, 1);
        xad[k].set(VAR2, 1);
        yad = aFcn<ttype>(xad, m->ne, m->e);
        third_ad[i][j][k] = yad.get(VAR0|VAR1|VAR2);
        xad[i] = x[i];
        xad[j] = x[j];
        xad[k] = x[k];
      }
    }
  }
}
  t2 = k_getTime();
  std::cout << "LibTaylor Third cost : " << t2 - t1 << std::endl;
#endif // COMPUTE_THIRD


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
        fcj = aFcn<double>(x,m->ne, m->e);
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

