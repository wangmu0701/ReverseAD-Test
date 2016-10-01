#include <cstdio>
#include <memory>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sys/time.h>

using namespace std;

#ifdef USING_ADOLC
#include "adolc/adolc.h"
typedef adouble scalar;
#else
typedef double scalar;
#endif

#include "stream_vel.hpp"

#ifdef USING_ADOLC
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
#ifdef USING_ADOLC
#if defined(CHECK_GRADIENT) || defined(CHECK_HESSIAN)
  trace_on(1);
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

#ifdef USING_ADOLC
  fc >>= fc_0;
#if defined(CHECK_GRADIENT) || defined(CHECK_HESSIAN)
  trace_off();
  init_double(dh, du);
#ifdef CHECK_GRADIENT
  double gradient_adolc[N];
  gettimeofday(&tv1, NULL);
  gradient(1, N, dh, gradient_adolc);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                        (tv2.tv_usec - tv1.tv_usec) / 1000000.0;
  std::cout << "ADOLC Gradient cost : " << time_elapsed << endl;

#endif
#ifdef CHECK_HESSIAN
  double** hessian_adolc = new double*[N];
  for (int i = 0; i < N; i++) {
    hessian_adolc[i] = new double[N];
  }
  gettimeofday(&tv1, NULL);
  hessian(1, N, dh, hessian_adolc); 
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) +
                        (tv2.tv_usec - tv1.tv_usec) / 1000000.0;
  std::cout << "ADOLC Hessian cost : " << time_elapsed << endl;

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
  cout<<setw(20)<<"ADOLC"<<setw(20)<<"RelError";
  cout<<endl<<setfill(' ')<<setprecision(12);
  for (int i = 0; i < N; i++) {
    init_double(dh, du);
    dh[i] += ep;
    //call stream_vel_timedep(h, u, bb, fc);
    stream_vel_timedep<double>(dh, du, fci);
    fdg = (fci - fc_0) / ep;
    cout<<setw(10)<<i<<setw(20)<<fdg;
    cout<<setw(20)<<gradient_adolc[i];
    double e = (fabs(fdg-gradient_adolc[i])/fabs(fdg+gradient_adolc[i]));
    cout<<setw(20)<<e;
    cout << endl;
  }
#endif


#ifdef CHECK_HESSIAN
  double fci, fcj, fcij;
  double fdh;
  cout<<setfill('-')<<setw(5)<<"i"<<setw(5)<<"j"<<setw(20)<<"FiniteDiff";
  cout<<setw(20)<<"ADOLC"<<setw(20)<<"RelError";
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
      //cout<<setw(20)<<fci<<setw(20)<<fcj<<setw(20)<<fcij<<setw(20)<<fc_0<<endl;
      cout<<setw(5)<<i<<setw(5)<<i<<setw(20)<<fdh;
      cout<<setw(20)<<hessian_adolc[i][j];
      double e = (fabs(fdh-hessian_adolc[i][j])/fabs(fdh+hessian_adolc[i][j]));
      cout<<setw(20)<<e;
      cout << endl;
}
    }
  }
#endif
*/
}
