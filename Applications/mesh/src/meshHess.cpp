#include "fcn.h"
#include <iostream>
#include <sys/time.h>

#ifdef USING_REVERSEAD
#include "reversead/reversead.hpp"
using namespace ReverseAD;
#endif

#ifdef USING_ADOLC
#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#endif

static double fTime=0.0;
static double hTime=0.0;
static double afTime=0.0;
static double ahTime=0.0;

#define rcbrt(x) pow(x,-3.333333333333333333333333333e-01)
#define sqrt3   5.77350269189625797959429519858e-01        /*  1.0/sqrt(3.0) */
#define sqrt6   4.08248290463863052509822647505e-01        /*  1.0/sqrt(6.0) */
#define a       3.33333333333333333333333333333e-01        /*  1.0/3.0       */
#define b      -6.66666666666666666666666666667e-01        /* -2.0/3.0       */

#if defined(USING_REVERSEAD) || defined(USING_ADOLC)
int a_fcn(adouble *obj, const adouble x[12])
{
  static adouble matr[9], f;
  static adouble g;
  f       = x[1] + x[0];
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - f)*sqrt3;
  matr[2] = (3.0*x[3] - x[2] - f)*sqrt6;

  f       = x[5] + x[4];
  matr[3] = x[5] - x[4];
  matr[4] = (2.0*x[6] - f)*sqrt3;
  matr[5] = (3.0*x[7] - x[6] - f)*sqrt6;

  f       = x[9] + x[8];
  matr[6] = x[9] - x[8];
  matr[7] = (2.0*x[10] - f)*sqrt3;
  matr[8] = (3.0*x[11] - x[10] - f)*sqrt6;

  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
      matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);
  if (g <= epsilonf) { *obj = g; return 1; }

  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  (*obj) = a * f * rcbrt(g*g);
  return 0;
}

int aFcn(adouble *obj, adouble *v,const Mesh *m)
{
  const int e = m->ne;
  adouble *w;
  int    *t = m->e;

  adouble  x[12];
  adouble  o, f;
  int     v1, v2, v3, v4;
  int     i;
  
  *obj = 0.0;

  o = 0.0;
  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];
    if (a_fcn(&f, x)) return 1;

    o += f;
  }
  *obj = o;
  return 0;
}
#endif

int main(int argc, char **argv){
    Mesh *m = NULL;
    char defaultFile[20]="gear.mesh";
    char* meshFile = defaultFile;
    if (argc >= 2) {
      meshFile = argv[1];
    }
    std::cout << "reading mesh .... ";
    if (readMesh(meshFile, &m)) {
        freeMesh(&m);
        return -1;
    }
    std::cout << "done" << std::endl;
    hMesh(m);

    struct timeval tv1,tv2;
    int nv=m->nv;
    int n=3*nv;
    double f;
gettimeofday(&tv1,NULL);
    oFcn(&f,m);
gettimeofday(&tv2,NULL);
    fTime+=(tv2.tv_sec-tv1.tv_sec)+(tv2.tv_usec-tv1.tv_usec)/1000000.0;

gettimeofday(&tv1,NULL);
    hOnly(m);
gettimeofday(&tv2,NULL);
    hTime+=(tv2.tv_sec-tv1.tv_sec)+(tv2.tv_usec-tv1.tv_usec)/1000000.0;

#if defined(USING_REVERSEAD) || defined(USING_ADOLC)
    adouble *v=new adouble[3*nv];
    adouble fObj;
    std::cout << "Overloaded function evaluation .... ";
gettimeofday(&tv1,NULL);
//Active Section
#ifdef USING_REVERSEAD
    trace_on<double>();
#else // ADOLC
    trace_on(1);
#endif
    for(int i = 0;i < 3 * nv; i++){
        v[i]<<=m->v[i];
    }
    aFcn(&fObj,v,m);
    fObj>>=f;
#ifdef USING_REVERSEAD
    std::shared_ptr<TrivialTrace<double>> trace = trace_off<double>();
#else // ADOLC
    trace_off(1);
#endif
gettimeofday(&tv2,NULL);
    std::cout << "done" << std::endl;
    afTime+=(tv2.tv_sec-tv1.tv_sec)+(tv2.tv_usec-tv1.tv_usec)/1000000.0;

gettimeofday(&tv1, NULL);
#ifdef USING_REVERSEAD
    BaseReverseHessian<double> hessian(trace);
    std::shared_ptr<DerivativeTensor<size_t, double>> tensor = hessian.compute(n, 1);
    size_t nnz;
    size_t** tind;
    double* values;
    tensor->get_internal_coordinate_list(0, 2, &nnz, &tind, &values);
#else // ADOLC
    int options[2] = {0, 1};
    unsigned int    *rind  = NULL;
    unsigned int    *cind  = NULL;
    double *values = NULL;
    int nnz;
    sparse_hess(1, n, 0, m->v, &nnz, &rind, &cind, &values, options);
#endif
gettimeofday(&tv2,NULL);
    ahTime+=(tv2.tv_sec-tv1.tv_sec)+(tv2.tv_usec-tv1.tv_usec)/1000000.0;
    std::cout << "n = " << n << std::endl;
    std::cout << "nnz = " << nnz << std::endl;
#endif

    printf("Total Plain      Function Time Estimate=<%10.6f>\n",fTime);
    printf("Total Analytic    Hessian Time Estimate=<%10.6f>\n",hTime);
    printf("Total Overloaded Function Time Estimate=<%10.6f>\n",afTime);
    printf("Total Overloaded  Hessian Time Estimate=<%10.6f>\n",ahTime);
    return 0;
}
