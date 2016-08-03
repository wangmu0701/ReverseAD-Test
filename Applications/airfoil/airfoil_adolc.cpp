#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<memory>
    
#include"const.h"

#include "adolc/adolc.h"
#include "adolc/adolc_sparse.h"

#include "routines_r.hpp"

double gam, gm1, cfl, eps, mach, alpha;    
    
int main()
{
  int i;
  const int maxnode = 721801;
  const int maxcell = 720001;
  const int maxedge = 1441400;
  int **ecell = new int*[maxedge];
  for (i = 0; i < maxedge; i++)
    ecell[i] = new int[2];

  int *boun = new int[maxedge];

  int **edge = new int*[maxedge];
  for (i = 0; i < maxedge; i++)
    edge[i] =  new int[2];

  int **cell = new int*[maxcell];
  for (i = 0; i < maxcell; i++)
    cell[i] = new int[4];

  double **x = new double*[maxnode];
  for (i = 0; i < maxnode; i++)
    x[i] = new double[2];

  double **q = new double*[maxcell];
  for (i = 0; i < maxcell; i++)
    q[i] = new double[4];

  double **qold = new double*[maxcell];
  for (i = 0; i < maxcell; i++)
    qold[i] = new double[4];

  double *adt = new double[maxcell];

  double **res = new double*[maxcell];
  for (i = 0; i < maxcell; i++)
    res[i] = new double[4];

  int in1, in2, in3, in4, ic, ic1, ic2, ie, ipde, iter, niter, k, nnode, ncell, nedge;
 
//------read in grid and flow data-------------
  input(maxnode, maxcell, maxedge, nnode, ncell, nedge, x, q, cell, edge, ecell, boun);
 
  std::cout << "nnode = " << nnode << std::endl;
  std::cout << "ncell = " << ncell << std::endl;
  std::cout << "nedge = " << nedge << std::endl;

  std::cout << "using adolc ..." << std::endl;
  typedef adouble active_type;
  active_type** x_ad = new active_type*[nnode];
  for (int i = 0; i < nnode; i++) 
    x_ad[i] = new active_type[2];

  active_type** q_ad = new active_type*[ncell];
  for (int i = 0; i < ncell; i++)
    q_ad[i] = new active_type[4];

  active_type** qold_ad = new active_type*[ncell];
  for (int i = 0; i < ncell; i++)
    qold_ad[i] = new active_type[4];

  active_type* adt_ad = new active_type[ncell];

  active_type** res_ad = new active_type*[ncell];
  for (int i = 0; i < ncell; i++)
    res_ad[i] = new active_type[4];

  active_type lift_ad, rms_ad;

// Independent : x and q

#ifdef COMPUTE_DERIVATIVE
  std::cout << "Will evaluate the derivative ... " << std::endl;
  trace_on(1);
#endif
  int num_ind = 0;
  for (int i = 0; i < nnode; i++) {
    x_ad[i][0] <<= x[i][0];
    x_ad[i][1] <<= x[i][1];
    num_ind += 2;
  }
  for (int i = 0; i < ncell; i++) {
    for (int j = 0; j < 4; j++) {
      q_ad[i][j] <<= q[i][j];
      num_ind++;
    }
  }
#ifdef COMPUTE_DERIVATIVE
  double* dummy_x = new double[num_ind];
  num_ind = 0;
  for (int i = 0; i < nnode; i++) {
    dummy_x[num_ind++] = x[i][0];
    dummy_x[num_ind++] = x[i][1];
  }
  for (int i = 0; i < ncell; i++) {
    for (int j = 0; j < 4; j++) {
      dummy_x[num_ind++] = q[i][j];
    }
  }
#endif

//----main time-marching loop------------------
  niter = MAX_ITER;
  for(iter = 0;iter < niter;++iter) {
 
//----save old flow solution--------------------
    for(ic=0; ic< ncell-1; ++ic) {
      for( ipde = 0;ipde < 4; ++ipde) {
        qold_ad[ic][ipde] = q_ad[ic][ipde];
      }
    }
//-----predictor/corrector update loop----------

    for(k = 0; k < 2 ; ++k) { 
      for(ic = 0; ic < ncell; ++ic) {
        for(ipde = 0;ipde < 4; ++ipde) {
          res_ad[ic][ipde] = 0.0;
        }
      }
//-----calculate area/timstep------------------
      for (ic = 0; ic < ncell - 1; ++ic) {
        in1 = cell[ic][0] -1;
        in2 = cell[ic][1] -1;
        in3 = cell[ic][2] -1;
        in4 = cell[ic][3] -1;
        time_cell<active_type>(&x_ad[in1][0], &x_ad[in2][0], &x_ad[in3][0],&x_ad[in4][0],&q_ad[ic][0],adt_ad[ic]);
      }

      adt[ncell] = 0.0;
//-----flux evaluation loop--------------------
      for (ie = 0;ie < nedge; ++ie) {
        in1 = edge[ie][0] -1;
        in2 = edge[ie][1] -1;
        ic1 = ecell[ie][0]-1;
        ic2 = ecell[ie][1]-1;
        if (boun[ie] == 0)  {
          flux_face<active_type>(&x_ad[in1][0], &x_ad[in2][0], &q_ad[ic1][0], &q_ad[ic2][0],adt_ad[ic1], adt_ad[ic2], &res_ad[ic1][0], &res_ad[ic2][0]);
        } else if (boun[ie] == 1) {
          flux_wall<active_type>(&x_ad[in1][0], &x_ad[in2][0], &q_ad[ic2][0], &res_ad[ic2][0]);
        } else if (boun[ie] == 2) {
          std::cout << "exit on boun[ie] == 2, something wrong?" << std::endl;
          exit(0);
        }
      }
    }
//-----flow field update-----------------------
    rms_ad = 0.0;
    for (ic = 0;ic < ncell - 1; ++ic) {
      for(ipde = 0;ipde < 4; ++ipde) {
        q_ad[ic][ipde] = qold_ad[ic][ipde] - res_ad[ic][ipde]/adt_ad[ic];
        rms_ad = rms_ad + (pow((q_ad[ic][ipde] - qold_ad[ic][ipde]),2));
      }
    }
    rms_ad = sqrt(rms_ad/ncell);

//-----print iteration history, including lift calculation
    if( (iter % PRINT_HISTORY_GAP) == 0) {
      lift_ad = 0.0;
      for(ie = 0;ie < nedge; ++ie) {
        if((boun[ie]) == 1) {
          in1 = edge[ie][0] -1;
          in2 = edge[ie][1] -1;
          ic2 = ecell[ie][1]-1;
          lift_wall<active_type>(&x_ad[in1][0],&x_ad[in2][0],&q_ad[ic2][0], lift_ad);
        }
      }
      std::cout << "iter = " << iter
                << ", rms = " << rms_ad.getValue()
                << ", lift = " << lift_ad.getValue() << std::endl;
    }
  }


#ifdef COMPUTE_DERIVATIVE
  // make lift to be the dependent variable
  double lift;
  lift_ad >>= lift;
  trace_off();
 
 
  unsigned int *rind = NULL;
  unsigned int *cind = NULL;
  double * values = NULL;
  int nnz;
  int options[2] = {0, 1}; // direct recovery
  sparse_hess(1, num_ind, 0, dummy_x, &nnz, &cind, &rind, &values, options);
  std::cout << "number of independes = " << num_ind << std::endl;
  std::cout << "size of Hessian = " << nnz << std::endl;
/*
// let's first compute gradient
  double* g = new double[num_ind];
  std::cout << "number of independes = " << num_ind << std::endl;
  gradient(1, num_ind, dummy_x, g);
  int g_nnz = 0;
  for (int i = 0; i < num_ind; i++) {
    if (g[i] != 0.0) {
      std::cout << "A["<<i<<"] = " << g[i] << std::endl;
      g_nnz++;
    }
  }
  std::cout << "size of gradient = " << g_nnz << std::endl;
*/
#endif

  for (ie = 0; ie < maxedge; ie++)
    delete[] ecell[ie];
  delete[] ecell;

  delete[] boun;

  for (i = 0; i < maxedge; i++)
    delete[] edge[i];
  delete[] edge;

  for (i = 0; i < maxcell; i++)
    delete[] cell[i];
  delete[] cell;

  for (i = 0; i < maxnode; i++)
    delete[] x[i];
  delete[] x;

  for (i = 0; i < maxcell; i++)
    delete[] q[i];
  delete[] q;

  for (i = 0; i < maxcell; i++)
    delete[] qold[i];
   delete[] qold;

  delete[] adt;

  for (i = 0; i < maxcell; i++)
    delete[] res[i];
  delete[] res;

  return 0;
}

