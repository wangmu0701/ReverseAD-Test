//----------------------------------------- program airfoil------------------------------------------------
    
    
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<memory>
    
#ifdef USING_REVERSEAD
#include "reversead/reversead.hpp"
#endif

#include "airfoil_func.hpp"

double gam, gm1, cfl, eps, mach, alpha;    
    
int main(int argc, char* argv[]) {
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

  int nnode, ncell, nedge; 
//------read in grid and flow data-------------
  input(maxnode, maxcell, maxedge, nnode, ncell, nedge, x, q, cell, edge, ecell, boun);

#ifdef USING_REVERSEAD
  std::cout << "using ReverseAD ..." << std::endl;
  typedef ReverseAD::adouble active_type;
  int num_ind = nnode*2+ncell*4;
#ifdef COMPUTE_DERIVATIVE
  std::cout << "Will evaluate the derivative ... " << std::endl;
  ReverseAD::trace_on<double>();
#endif
  active_type* dummy_x = new active_type[num_ind];
  int l = 0;
  for (int i = 0; i < nnode; i++) {
    dummy_x[l++] <<= x[i][0];
    dummy_x[l++] <<= x[i][1];
  }
  for (int i = 0; i < ncell; i++) {
    for (int j = 0; j < 4; j++) {
      dummy_x[l++] <<= q[i][j];
    }
  }
#else
  int num_ind = nnode*2+ncell*4;
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
 
  typedef double active_type;
#endif

  active_type lift_ad = airfoil_func<active_type>(dummy_x, num_ind, nnode, ncell, nedge, ecell, boun, edge, cell);


#ifdef USING_REVERSEAD

#ifdef COMPUTE_DERIVATIVE
  // make lift to be the dependent variable
  double lift;
  lift_ad >>= lift;

  std::shared_ptr<ReverseAD::TrivialTrace<double>> trace = ReverseAD::trace_off<double>();
  std::cout << "Function tracing done" << std::endl;

  ReverseAD::BaseReverseHessian<double> hessian(trace);
  std::shared_ptr<ReverseAD::DerivativeTensor<size_t, double>> tensor =
      hessian.compute(num_ind, 1);
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
  fp.open("hessian.mm");
  // banner
  fp << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
  fp << num_ind << " " << num_ind << " " << size << std::endl;
  for (int i = 0; i < size; i++) {
    fp << tind[i][0]+1 << " " << tind[i][1]+1 << " " << values[i] << std::endl; 
  }
  fp.close();
#endif
#else 

#ifdef HESSIAN_CHECK
// Now we can run some checks
  if (argc < 3) {
    std::cout << "Usage: ./check_hess rind cind [eps (optional)]" << std::endl;
    exit(-1);
  }
  int r = atoi(argv[1]);
  int c = atoi(argv[2]);
  double eps = 1.0e-8;
  if (argc > 3) {
    eps = atof(argv[3]);
  } 
  double fxy = lift_ad;
  dummy_x[r]+=eps;
  double fdxy = airfoil_func<double>(dummy_x, num_ind, nnode, ncell, nedge, ecell, boun, edge, cell);
  dummy_x[c]+=eps;
  double fdxdy = airfoil_func<double>(dummy_x, num_ind, nnode, ncell, nedge, ecell, boun, edge, cell);
  dummy_x[r]-=eps;
  double fxdy = airfoil_func<double>(dummy_x, num_ind, nnode, ncell, nedge, ecell, boun, edge, cell);
  double he = (fxy+fdxdy-fdxy-fxdy)/(eps*eps);
  printf("fxy = %.10f\n", fxy);
  printf("fdxy = %.10f\n", fdxy);
  printf("fxdy = %.10f\n", fxdy);
  printf("fdxdy = %.10f\n", fdxdy);
  std::cout << "Finite difference H["<<r<<"]["<<c<<"] = " << he << std::endl;
#endif
#endif

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

