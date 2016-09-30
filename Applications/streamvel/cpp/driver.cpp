#include <cstdio>
#include <iostream>
#include "stream_vel.hpp"

void stream_vel_timedep(double*, double*, double*, double&);

int main() {
  double bb[N];
  double h[N];
  double u[N+1];
  double fc;
  double fc_0, accuracyAD, fdfc;
  double ep = 1.0e-4;
  for (int i = 0; i < N; i++) {
    bb[i] = 0.0;
  }
  u[N] = 0.0;
  fc = 0.0;
  // call strem_vel_timedep(h, u, bb, fc);
  stream_vel_timedep(h, u, bb, fc);

  fc_0 = fc;
  printf("  position  velocity   thickness \n");
  for (int i = 0; i < N; i++) {
    printf("%d  %.10f  %.10f\n", i, u[i+1], h[i]);
  }

  printf(" position  dfc/DM [ADM], dfc/DM[f.d.] relative accuracy");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      bb[j] = 0.0;
    }
    bb[i] = ep;
    //call stream_vel_timedep(h, u, bb, fc);
    stream_vel_timedep(h, u, bb, fc);
    fdfc = (fc - fc_0) / ep;
    printf("%d  %.10f\n", i, fdfc);
  }

}
