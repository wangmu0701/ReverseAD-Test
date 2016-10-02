#include <cstdio>
#include <memory>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sys/time.h>

using namespace std;

#ifdef USING_REVERSEAD
#include "reversead/reversead.hpp"
using namespace ReverseAD;
typedef adouble scalar;
#else
typedef double scalar;
#endif

#include "stream_vel.hpp"

#ifdef USING_REVERSEAD
void init_adouble(adouble* h, adouble* u) {
  for (int i = 0; i < N; i++) {
    h[i] <<= h_left + (h_right - h_left) / Lx * ((i+1) - 0.5) * dx;

  }
  for (int i = 0; i <= N; i++) {
    u[i] = 0.0;
  }
}
#endif

void init_double(double* h, double* u) {
  for (int i = 0; i < N; i++) {
    h[i] = h_left + (h_right - h_left) / Lx * ((i+1) - 0.5) * dx;
    u[i] = 0.0;
  }
  u[N] = 0.0;
}
int main() {
  scalar h[N];
  scalar u[N+1];
  scalar fc;

  double fc_0;
  double ep = 1.0e-4;
  double dh[N];
  double du[N+1];
  struct timeval tv1, tv2;
  double time_elapsed;
  gettimeofday(&tv1, NULL);
#ifdef USING_REVERSEAD
#if defined(REVERSEAD_GRADIENT) || defined(REVERSEAD_HESSIAN)
  trace_on<double>();
#endif
  init_adouble(h, u); 
#else
  init_double(h, u);
#endif
  // call strem_vel_timedep(h, u, bb, fc);
  stream_vel_timedep<scalar>(h, u, fc);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                        (tv2.tv_usec - tv1.tv_usec) / 1000000.0;
  std::cout << "Function evaluation cost : " << time_elapsed << endl;

#ifdef USING_REVERSEAD
  fc >>= fc_0;
#if defined(REVERSEAD_GRADIENT) || defined(REVERSEAD_HESSIAN)
  shared_ptr<TrivialTrace<double>> trace = trace_off<double>();
#ifdef REVERSEAD_GRADIENT
  double gradient_reversead[N];
  gettimeofday(&tv1, NULL);
  BaseReverseAdjoint<double> adjoint(trace);
  shared_ptr<DerivativeTensor<size_t, double>> tensor = adjoint.compute(N, 1);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                        (tv2.tv_usec - tv1.tv_usec) / 1000000.0;
  std::cout << "ReverseAD Gradient cost : " << time_elapsed << endl;
  size_t size;
  size_t** tind;
  double* values;
  tensor->get_internal_coordinate_list(0, 1, &size, &tind, &values);
  for (int i = 0; i < N; i++) { gradient_reversead[i] = 0.0;}
  for (int i = 0; i < size; i++) {
    //cout << "A["<<tind[i][0]<<"] = " << values[i] << endl;
    if (tind[i][0] < N) {
      gradient_reversead[tind[i][0]] = values[i];
    }
  }
#endif
#ifdef REVERSEAD_HESSIAN
  double hessian_reversead[N][N];
  gettimeofday(&tv1, NULL);
  BaseReverseHessian<double> hessian(trace);
  shared_ptr<DerivativeTensor<size_t, double>> tensor = hessian.compute(N, 1);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                        (tv2.tv_usec - tv1.tv_usec) / 1000000.0;
  std::cout << "ReverseAD Hessian cost : " << time_elapsed << endl;
  size_t size;
  size_t** tind;
  double* values;
  tensor->get_internal_coordinate_list(0, 2, &size, &tind, &values);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      hessian_reversead[i][j] = 0;
    }
  }
  for (int i = 0; i < size; i++) {
    //cout << "H["<<tind[i][0]<<","<<tind[i][1]<<"] = " << values[i]<<endl;
    if (tind[i][0] < N && tind[i][1] < N) {
      hessian_reversead[tind[i][0]][tind[i][1]] = values[i];
    }
  }
#endif

#endif
#else
  fc_0 = fc;
#endif
/*
  cout<<setfill('-')<<setw(10)<<"position"<<setw(20)<<"velocity"<<setw(20)<<"thickness"<<endl;
  cout<<setfill(' ')<<setprecision(12);
  for (int i = 0; i < N; i++) {
    cout<<setw(10)<<i<<setw(20)<<h[i]<<setw(20)<<u[i+1]<<endl; 
  }

#ifdef CHECK_GRADIENT
  double fci, fdg;
  cout<<setfill('-')<<setw(10)<<"position"<<setw(20)<<"FiniteDiff";
#ifdef REVERSEAD_GRADIENT
  cout<<setw(20)<<"ReverseAD"<<setw(20)<<"RelError";
#endif
  cout<<endl<<setfill(' ')<<setprecision(12);
  for (int i = 0; i < N; i++) {
    init_double(dh, du);
    dh[i] += ep;
    //call stream_vel_timedep(h, u, bb, fc);
    stream_vel_timedep<double>(dh, du, fci);
    fdg = (fci - fc_0) / ep;
    cout<<setw(10)<<i<<setw(20)<<fdg;
#ifdef REVERSEAD_GRADIENT
    cout<<setw(20)<<gradient_reversead[i];
    scalar e = (fabs(fdg-gradient_reversead[i])/fabs(fdg+gradient_reversead[i]));
    cout<<setw(20)<<e;
#endif
    cout << endl;
  }
#endif


#ifdef CHECK_HESSIAN
  double fci, fcj, fcij;
  double fdh;
  cout<<setfill('-')<<setw(5)<<"i"<<setw(5)<<"j"<<setw(20)<<"FiniteDiff";
#ifdef REVERSEAD_HESSIAN
  cout<<setw(20)<<"ReverseAD"<<setw(20)<<"RelError";
#endif
  cout<<endl<<setfill(' ')<<setprecision(12);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j <=i; j++) {
      init_double(dh, du);
      dh[i] += ep;
      //call stream_vel_timedep(h, u, bb, fc);
      stream_vel_timedep<double>(dh, du, fci);
      dh[j] += ep;
      stream_vel_timedep<double>(dh, du, fcij);
      dh[i] -= ep;
      stream_vel_timedep<double>(dh, du, fcj);
      fdh = (fcij-fci-fcj+fc_0) / (ep * ep);
if (fabs(i-j) < 2) {
      cout<<setw(5)<<i<<setw(5)<<i<<setw(20)<<fdh;
#ifdef REVERSEAD_HESSIAN
      cout<<setw(20)<<hessian_reversead[i][j];
      scalar e = (fabs(fdh-hessian_reversead[i][j])/fabs(fdh+hessian_reversead[i][j]));
      cout<<setw(20)<<e;
#endif
      cout << endl;
}
    }
  }
#endif
*/
}
