//----------------------------------------- program airfoil------------------------------------------------
    
    
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<memory>
    
#include "reversead/reversead.hpp"
using namespace ReverseAD;

#include "const.h"
#include "routines_r.hpp"

double gam, gm1, cfl, eps, mach, alpha;    
int nnode;
int ncell;
int nedge;
int** ecell;
int* boun;
int** edge;
int** cell;
int iter = 0;

adouble** x_ad;
adouble** q_ad;
adouble** qold_ad;
adouble* adt_ad;
adouble** res_ad;

void dummy_func(){
}

void airfoil_setup(adouble* dummy_x, size_t num_ind, adouble* dummy_t, size_t num_t) {

  std::cout << "nnode = " << nnode << std::endl;
  std::cout << "ncell = " << ncell << std::endl;
  std::cout << "nedge = " << nedge << std::endl;
  for (int i = 0; i < num_ind; i++) {
    dummy_t[i] = dummy_x[i];
  }
  x_ad = new adouble*[nnode];
  for (int i = 0; i < nnode; i++)
    x_ad[i] = new adouble[2];

  q_ad = new adouble*[ncell];
  for (int i = 0; i < ncell; i++)
    q_ad[i] = new adouble[4];

  qold_ad = new adouble*[ncell];
  for (int i = 0; i < ncell; i++)
    qold_ad[i] = new adouble[4];

  adt_ad = new adouble[ncell+1];

  res_ad = new adouble*[ncell];
  for (int i = 0; i < ncell; i++)
    res_ad[i] = new adouble[4];

}

//----main time-marching loop------------------
void airfoil_run(adouble* dummy_t, size_t num_t) {
  int l = 0;
  for (int i = 0; i < nnode; i++) {
    x_ad[i][0] = dummy_t[l++];
    x_ad[i][1] = dummy_t[l++];
  }
  for (int i = 0; i < ncell; i++) {
    for (int j = 0; j < 4; j++) {
      q_ad[i][j] = dummy_t[l++];
    }
  }

  int i, k, ic, ipde, in1, in2, in3, in4, ie, ic1, ic2;
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
        time_cell<adouble>(&x_ad[in1][0], &x_ad[in2][0], &x_ad[in3][0],&x_ad[in4][0],&q_ad[ic][0],adt_ad[ic]);
      }

      adt_ad[ncell] = 0.0;
//-----flux evaluation loop--------------------
      for (ie = 0;ie < nedge; ++ie) {
        in1 = edge[ie][0] -1;
        in2 = edge[ie][1] -1;
        ic1 = ecell[ie][0]-1;
        ic2 = ecell[ie][1]-1;
        if (boun[ie] == 0)  {
          flux_face<adouble>(&x_ad[in1][0], &x_ad[in2][0], &q_ad[ic1][0], &q_ad[ic2][0],adt_ad[ic1], adt_ad[ic2], &res_ad[ic1][0], &res_ad[ic2][0]);
        } else if (boun[ie] == 1) {
          flux_wall<adouble>(&x_ad[in1][0], &x_ad[in2][0], &q_ad[ic2][0], &res_ad[ic2][0]);
        } else if (boun[ie] == 2) {
          std::cout << "exit on boun[ie] == 2, something wrong?" << std::endl;
          exit(0);
        }
      }
    }
//-----flow field update-----------------------
    adouble rms_ad = 0.0;
    for (ic = 0;ic < ncell - 1; ++ic) {
      for(ipde = 0;ipde < 4; ++ipde) {
        q_ad[ic][ipde] = qold_ad[ic][ipde] - res_ad[ic][ipde]/adt_ad[ic];
        rms_ad = rms_ad + (pow((q_ad[ic][ipde] - qold_ad[ic][ipde]),2));
      }
    }
    rms_ad = sqrt(rms_ad/ncell);


//-----print iteration history, including lift calculation
    if( (iter % PRINT_HISTORY_GAP) == 0) {
      adouble lift_ad = 0.0;
      for(ie = 0;ie < nedge; ++ie) {
        if((boun[ie]) == 1) {
          in1 = edge[ie][0] -1;
          in2 = edge[ie][1] -1;
          ic2 = ecell[ie][1]-1;
          lift_wall<adouble>(&x_ad[in1][0],&x_ad[in2][0],&q_ad[ic2][0], lift_ad);
        }
      }

      std::cout << "iter = " << iter++
                << ", rms = " << rms_ad
                << ", lift = " << lift_ad << std::endl;

    }
    l = 0;
    for (int i = 0; i < nnode; i++) {
      dummy_t[l++] = x_ad[i][0];
      dummy_t[l++] = x_ad[i][1];
    }
    for (int i = 0; i < ncell; i++) {
      for (int j = 0; j < 4; j++) {
        dummy_t[l++] = q_ad[i][j];
      }
    }
}
void airfoil_final(adouble* dummy_t, size_t num_t, adouble* dummy_y, size_t num_y) {
  int l = 0;
  for (int i = 0; i < nnode; i++) {
    x_ad[i][0] = dummy_t[l++];
    x_ad[i][1] = dummy_t[l++];
  }
  for (int i = 0; i < ncell; i++) {
    for (int j = 0; j < 4; j++) {
      q_ad[i][j] = dummy_t[l++];
    }
  }
  int i, ie, in1, in2, ic2;
  adouble lift_ad = 0.0;
  for(ie = 0;ie < nedge; ++ie) {
    if((boun[ie]) == 1) {
      in1 = edge[ie][0] -1;
      in2 = edge[ie][1] -1;
      ic2 = ecell[ie][1]-1;
      lift_wall<adouble>(&x_ad[in1][0],&x_ad[in2][0],&q_ad[ic2][0], lift_ad);
    }
  }

  dummy_y[0] = lift_ad;

  std::cout << "final lift = " << dummy_y[0] << std::endl;
}

void tear_down() {
  int i;
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
} 

int main(int argc, char* argv[]) {
  int i;
  const int maxnode = 721801;
  const int maxcell = 720001;
  const int maxedge = 1441400;
  ecell = new int*[maxedge];
  for (i = 0; i < maxedge; i++)
    ecell[i] = new int[2];

  boun = new int[maxedge];

  edge = new int*[maxedge];
  for (i = 0; i < maxedge; i++)
    edge[i] =  new int[2];

  cell = new int*[maxcell];
  for (i = 0; i < maxcell; i++)
    cell[i] = new int[4];

  double **x = new double*[maxnode];
  for (i = 0; i < maxnode; i++)
    x[i] = new double[2];

  double **q = new double*[maxcell];
  for (i = 0; i < maxcell; i++)
    q[i] = new double[4];

//------read in grid and flow data-------------
  input(maxnode, maxcell, maxedge, nnode, ncell, nedge, x, q, cell, edge, ecell, boun);

  std::cout << "using IterativeFuncFixed ..." << std::endl;
  typedef ReverseAD::adouble active_type;
  size_t num_ind = nnode*2+ncell*4;
  double* dummy_x = new double[num_ind];
  int l = 0;
  for (int i = 0; i < nnode; i++) {
    dummy_x[l++] = x[i][0];
    dummy_x[l++] = x[i][1];
  }
  for (int i = 0; i < ncell; i++) {
    for (int j = 0; j < 4; j++) {
      dummy_x[l++] = q[i][j];
    }
  }

  IterativeFuncFixed iter_func_fixed(
      num_ind, num_ind, 1,
      &dummy_func, &tear_down,
      &(airfoil_setup), &(airfoil_run), &(airfoil_final), MAX_ITER);
  double lift = 0;
  //iter_func_fixed.run(dummy_x, num_ind, &lift, 1);
  //exit(-1);
  std::shared_ptr<ReverseAD::DerivativeTensor<size_t, double>> tensor =
      iter_func_fixed.compute(dummy_x, num_ind, 2);
  size_t size;
  size_t** tind;
  double* values;
  std::cout << "number of independents = " << num_ind << std::endl;
  tensor->get_internal_coordinate_list(0, 1, &size, &tind, &values);
  std::cout << "size of adjoints = " << size << std::endl;

/*
  for (int i = 0; i < size; i++) {
    std::cout << "A["<<tind[i][0] << "] = " << values[i] << std::endl;
  }
*/
  tensor->get_internal_coordinate_list(0, 2, &size, &tind, &values);
  std::cout << "size of hessian = " << size << std::endl;
/*
  double* hv = new double[num_ind];
  for (int i = 0; i < num_ind; i++) {
    hv[i] = 0;
  }
  for (int i = 0; i < size; i++) {
    if (tind[i][0] != tind[i][1]) {
      hv[tind[i][0]] += values[i] * (tind[i][1]+1.0);
      hv[tind[i][1]] += values[i] * (tind[i][0]+1.0);
    } else {
      hv[tind[i][0]] += values[i] * (tind[i][0]+1.0);
    }
  }
  std::ofstream fp;
  fp.open("hessian-v.mm");
  for (int i = 0; i < size; i++) {
    if (hv[i] != 0.0) {
      fp << "hv["<<i<<"] = " << hv[i] << std::endl;
    }
  }
  fp.close();
*/

  std::ofstream fp;
  fp.open("hessian-i.mm");
  // banner
  fp << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
  fp << num_ind << " " << num_ind << " " << size << std::endl;
  for (int i = 0; i < size; i++) {
    fp << tind[i][0]+1 << " " << tind[i][1]+1 << " " << values[i] << std::endl; 
  }
  fp.close();


  for (int ie = 0; ie < maxedge; ie++)
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

  return 0;
}

