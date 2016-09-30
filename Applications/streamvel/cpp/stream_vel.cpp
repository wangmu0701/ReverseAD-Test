#include <cstdio>
#include <memory>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "stream_vel.hpp"

using namespace std;

#ifdef USING_REVERSEAD
using namespace ReverseAD;

void init_adouble(adouble* h, adouble* u) {
  for (int i = 0; i < N; i++) {
    h[i] <<= h_left + (h_right - h_left) / Lx * ((i+1) - 0.5) * dx;

  }
  for (int i = 0; i <= N; i++) {
    u[i] <<= 0.0;
  }
}
#else
void init_double(double* h, double* u) {
  for (int i = 0; i < N; i++) {
    h[i] = h_left + (h_right - h_left) / Lx * ((i+1) - 0.5) * dx;
    u[i] = 0.0;
  }
  u[N] = 0.0;
}
#endif
int main() {
  scalar h[N];
  scalar u[N+1];
  scalar fc;
  scalar fdfc;
  double fc_0;
  double gradient_reversead[N];

#ifdef USING_REVERSEAD
#ifdef CHECK_GRADIENT
  trace_on<double>();
#endif
  init_adouble(h, u); 
#else
  init_double(h, u);
#endif
  // call strem_vel_timedep(h, u, bb, fc);
  stream_vel_timedep<scalar>(h, u, fc);
#ifdef USING_REVERSEAD
  fc >>= fc_0;

#ifdef CHECK_GRADIENT
  shared_ptr<TrivialTrace<double>> trace = trace_off<double>();
  BaseReverseAdjoint<double> adjoint(trace);
  shared_ptr<DerivativeTensor<size_t, double>> tensor = adjoint.compute(N*2+1, 1);
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
#else
  fc_0 = fc;
#endif
  cout<<setfill('-')<<setw(10)<<"position"<<setw(20)<<"velocity"<<setw(20)<<"thickness"<<endl;
  cout<<setfill(' ')<<setprecision(12);
  for (int i = 0; i < N; i++) {
    cout<<setw(10)<<i<<setw(20)<<h[i]<<setw(20)<<u[i+1]<<endl; 
  }

#ifdef CHECK_GRADIENT
  scalar ep = 1.0e-4;
  cout<<setfill('-')<<setw(10)<<"position"<<setw(20)<<"FiniteDiff"<<setw(20)<<"ReverseAD"<<setw(20)<<"RelError"<<endl;
  cout<<setfill(' ')<<setprecision(12);
  for (int i = 0; i < N; i++) {
#ifdef USING_REVERSEAD
    init_adouble(h, u);
#else
    init_double(h, u);
#endif
    h[i] += ep;
    //call stream_vel_timedep(h, u, bb, fc);
    stream_vel_timedep(h, u, fc);
    fdfc = (fc - fc_0) / ep;
    cout<<setw(10)<<i<<setw(20)<<fdfc;
#ifdef CHECK_GRADIENT
    cout<<setw(20)<<gradient_reversead[i];
    scalar e = (fabs(fdfc-gradient_reversead[i])/fabs(fdfc+gradient_reversead[i]));
    cout<<setw(20)<<e;
#endif
    cout << endl;
  }
#endif
}
