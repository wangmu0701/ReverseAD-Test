#include <cmath>
#include <iostream>

#include "stream_vel.hpp"

int n_nl = 60;
int n_timesteps = 10;
double dt = 0.04;
double Lx = 79.e3;
double tol = 1.e-7;
double adjtol = 1.e-7;
double ep_glen = 1.e-7;
double eps = 1.e-5;
double Aglen = 5.0002e-17;
double nglen = 3.0;
double g = 9.81;
double rhoi = 910.0;
double rhow = 1035.0;
double R_bed = -900;
double beta_const = 5.0;
double h_left = 1050.0;
double h_right = 1050.0;
constexpr double PI = 3.14159265358979323844;
double r[N];
double r_old[N];
double p[N];
double p_old[N];
double ax[N];
double x_old[N];
double tridiag_0[N][3];
double dx;
int isinloop0;
int isinloop1;
int isinloop2;

void stream_vel_init(double*, double*);
void stream_vel_timedep_stage(double*, double*, double*);
void stream_vel(double*, double*, double*);
void stream_vel_taud(double*, double*, double&);
void stream_vel_visc(double*, double*, double*);
void stream_vel_assemble(double*, double*, double[N][3]);
void forward_step(double*, double*, double*);
void phistage(double*, double*, double*, double*, double*);
void phi(double*, double*, double*, double*, double*);
void solve(double*, double*, double[N][3]);

// here reference on fc mimics the intent(inout) for non-pointer types
void stream_vel_timedep(double* h, double* u, double* bb, double& fc) {
  double h0[N];
  double beta_0[N];
  double beta_fric[N];
  // call stream_vel_init(h0, beta_0);
  stream_vel_init(h0, beta_0);

  // beta_fric = beta_0
  for (int i = 0; i < N; i++) {beta_fric[i] = beta_0[i];}
  // h = h0+bb
  for (int i = 0; i < N; i++) {h[i] = h0[i] + bb[i];}
  // call stream_vel(u, h, beta_fric);

  stream_vel(u, h, beta_fric);
  // call stream_vel_timedep_stage(h, u, beta_fric);

  stream_vel_timedep_stage(h, u, beta_fric);

  fc = 0; 
  for (int i = 1; i <= N; i++) {
    fc += u[i]*u[i];
  } 

}

void forward_step(double* h, double* u, double* beta_fric) {
  h[0] = h[0] - dt/dx*u[1]*h[1]; 
  for (int i = 1; i < N; i++) {
    h[i] = h[i] - dt/dx*(u[i+1]*h[i] - u[i]*h[i-1]);
  }
  stream_vel(u, h, beta_fric);
}

void stream_vel_timedep_stage(double* h, double* u, double* beta_fric) {
  for (int i =0 ; i < n_timesteps; i++ ){
    forward_step(h, u, beta_fric);
  }
}


void stream_vel_init(double* h, double* beta) {
  dx = Lx / double(N);
  for (int i = 0; i < N; i++) {
    beta[i] = beta_const;
    h[i] = h_left + (h_right - h_left) / Lx * ((i+1) - 0.5) * dx;
  }
}

void stream_vel(double* u, double* h, double* beta_fric) {
  double f[N];
  double utmp[N];
  double b[N];
  double unew[N+1];
  double fend;
  
  isinloop0 = 0;
  isinloop1 = 1;
  isinloop2 = 2;
  
  // call stream_vel_taud(h, f, fend);
  stream_vel_taud(h, f, fend);

  for (int i = 0; i <= N; i++) {u[i] = 0;}
// ----- driving stress -----

  for (int i = 0; i < N; i++) {
    b[i] = -dx * f[i];
    if (i < N-1) {
      b[i] = b[i] - f[i+1] * dx;
    }
  }
  b[N-1] = b[N-1] + fend;

// --------------------------
  // call phistage(u, unew, b, h, beta_fric, ininloop0)
  phistage(u, unew, b, h, beta_fric); 
  // u = unew;
  for (int i = 0; i <=N; i++) {u[i] = unew[i];}

  // I dont' underwhat is isinloop does in original fortran code
  // But it seems that only one call of the phistage is enough?
}

void phistage(double* u, double* u_ip1, double* b, double* h,
              double* beta_fric) {
  phi(u, u_ip1, b, h, beta_fric);
}

void phi(double* u_i, double* u_ip1, double* b, double* h, double* beta_fric) {
  double nu[N];
  double utmp[N];
  double A[N][3];
  // stream_vel_visc(h, u_i, nu)
  stream_vel_visc(h, u_i, nu);
  
  stream_vel_assemble(nu, beta_fric, A);
  for (int i = 0; i < N; i++) {utmp[i] = 0.0;}
  
  // solve(utmp, b, A);
  solve(utmp, b, A);
  for (int i = 0; i < N; i++) {
    u_ip1[i+1] = utmp[i];
  }
}

void stream_vel_visc(double* h, double* u, double* nu) {
  double ux, tmp;
  for (int i = 0; i < N; i++) {
    ux = (u[i+1] - u[i]) / dx;
    tmp = ux * ux + ep_glen * ep_glen;
    nu[i] = 0.5*h[i]*pow(Aglen, -1.0/nglen)*pow(tmp, ((1-nglen)/2./nglen));
  }

}

void stream_vel_assemble(double* nu, double* beta_fric, double A[N][3]) {
  for (int i = 0; i < N; i++) {
    A[i][0] = A[i][2] = 0.0;
    A[i][1] = 4*nu[i]/dx + dx/3.0 * beta_fric[i] * beta_fric[i];
    if (i > 0) {
      A[i][0] = -4*nu[i]/dx + dx/6.0 * beta_fric[i] * beta_fric[i];
    }
    if (i < N-1) {
      A[i][1] = A[i][1] + 4*nu[i+1]/dx+dx/3.0*beta_fric[i+1]*beta_fric[i+1];
      A[i][2] = -4*nu[i+1]/dx + dx/6.0*beta_fric[i+1]*beta_fric[i+1];
    }
  }
}

void stream_vel_taud(double* h, double* f, double& fend) {
  for (int i = 0; i < N; i++) {
    if (i > 0 && i < N-1) {
      f[i] = rhoi * g * h[i] * (h[i+1] - h[i-1]) / 2.0 / dx;
    } else if (i == 0) {
      f[i] = rhoi * g * h[i] * (h[i+1] - h[i]) / dx;
    } else if (i == N-1) {
      f[i] = rhoi * g * h[i] * (h[i] - h[i-1]) / dx;
    }
  }
  fend = 0.5 * (rhoi * g * h[N-1] * h[N-1] - rhow * g * R_bed*R_bed);
}


void solve(double* x, double* b, double A[N][3]) {
  double alpha, beta, dp1, dp2, res_init, resid;
  int k_iter = 0;
  for (int i = 0; i < N; i++) {
    r[i] = r_old[i] = p[i] = p_old[i] = x[i] = 0.0;
  }
  for (int i = 0; i < N; i++) {
    r[i] = b[i] - A[i][1] * x[i];
    if (i > 0) {
      r[i] = r[i] - A[i][0] * x[i-1];
    }
    if (i < N-1) {
      r[i] = r[i] - A[i][2] * x[i+1];
    }
  }  
  dp1 = 0;
  for (int i = 0; i < N; i++) {
    dp1 += r[i]*r[i];
  }
  res_init = sqrt(dp1);
  resid = res_init;
  for (int i = 0; i < N; i++) {
    p[i] = r[i];
  }

  while (k_iter < 200 * N && resid > 1.0e-10*res_init) {
    k_iter = k_iter +1;
    for (int i = 0; i < N; i++) {
      ax[i] = A[i][1] * p[i];
      if (i > 0) {
        ax[i] = ax[i] + A[i][0] * p[i-1];
      }
      if (i < N-1) {
        ax[i] = ax[i] + A[i][2] * p[i+1];
      }
    }
    dp1 = 0;
    dp2 = 0;
    for (int i = 0; i < N; i++) {
      dp1 += r[i] * r[i];
      dp2 += p[i] * ax[i];
    }
    alpha = dp1 / dp2;

    for (int i = 0; i < N; i++) {
      r_old[i] = r[i];
      x[i] = x[i] + alpha*p[i];
      r[i] = r[i] - alpha*ax[i];
    }
/*
    dp1 = 0;
    for (int i = 0; i < N; i++) {
      dp1 += r[i] * r[i];
    } 
    resid = sqrt(dp1);
*/
    dp1 = 0;
    dp2 = 0;
    for (int i = 0; i < N; i++) {
      dp1 += r[i] * r[i];
      dp2 += r_old[i] * r_old[i];
    }
    beta = dp1 / dp2;
    resid = sqrt(dp1);    

    for (int i = 0; i < N; i++) {
      p[i] = r[i] + beta * p[i];
    }
    
  }
}
