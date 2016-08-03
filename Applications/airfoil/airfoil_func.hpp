//----------------------------------------- program airfoil------------------------------------------------
    
    
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<memory>
    
#include"const.h"

#ifdef USING_REVERSEAD
#include "reversead/reversead.hpp"
#endif

#include "routines_r.hpp"

template <typename T>
T airfoil_func(T* dummy_x, int num_ind, int nnode, int ncell, int nedge,
                  int** ecell, int* boun, int** edge, int** cell) {
  int in1, in2, in3, in4, ic, ic1, ic2, ie, ipde, iter, i, k;
  T rms_ad, lift_ad;

  std::cout << "nnode = " << nnode << std::endl;
  std::cout << "ncell = " << ncell << std::endl;
  std::cout << "nedge = " << nedge << std::endl;

  T** x_ad = new T*[nnode];
  for (int i = 0; i < nnode; i++) 
    x_ad[i] = new T[2];

  T** q_ad = new T*[ncell];
  for (int i = 0; i < ncell; i++)
    q_ad[i] = new T[4];

  T** qold_ad = new T*[ncell];
  for (int i = 0; i < ncell; i++)
    qold_ad[i] = new T[4];

  T* adt_ad = new T[ncell+1];

  T** res_ad = new T*[ncell];
  for (int i = 0; i < ncell; i++)
    res_ad[i] = new T[4];

  int l = 0;
  for (int i = 0; i < nnode; i++) {
    x_ad[i][0] = dummy_x[l++];
    x_ad[i][1] = dummy_x[l++];
  }
  for (int i = 0; i < ncell; i++) {
    for (int j = 0; j < 4; j++) {
      q_ad[i][j] = dummy_x[l++];
    }
  }

//----main time-marching loop------------------
  for(iter = 0;iter < MAX_ITER; ++iter) {
 
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
        time_cell<T>(&x_ad[in1][0], &x_ad[in2][0], &x_ad[in3][0],&x_ad[in4][0],&q_ad[ic][0],adt_ad[ic]);
      }

      adt_ad[ncell] = 0.0;
//-----flux evaluation loop--------------------
      for (ie = 0;ie < nedge; ++ie) {
        in1 = edge[ie][0] -1;
        in2 = edge[ie][1] -1;
        ic1 = ecell[ie][0]-1;
        ic2 = ecell[ie][1]-1;
        if (boun[ie] == 0)  {
          flux_face<T>(&x_ad[in1][0], &x_ad[in2][0], &q_ad[ic1][0], &q_ad[ic2][0],adt_ad[ic1], adt_ad[ic2], &res_ad[ic1][0], &res_ad[ic2][0]);
        } else if (boun[ie] == 1) {
          flux_wall<T>(&x_ad[in1][0], &x_ad[in2][0], &q_ad[ic2][0], &res_ad[ic2][0]);
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
          lift_wall<T>(&x_ad[in1][0],&x_ad[in2][0],&q_ad[ic2][0], lift_ad);
        }
      }

      std::cout << "iter = " << iter
                << ", rms = " << rms_ad
                << ", lift = " << lift_ad << std::endl;

    }
  }

  for (i = 0; i < nnode; i++)
    delete[] x_ad[i];
  delete[] x_ad;

  for (i = 0; i < ncell; i++)
    delete[] q_ad[i];
  delete[] q_ad;

  for (i = 0; i < ncell; i++)
    delete[] qold_ad[i];
   delete[] qold_ad;

  delete[] adt_ad;

  for (i = 0; i < ncell; i++)
    delete[] res_ad[i];
  delete[] res_ad;

  return lift_ad;
}
