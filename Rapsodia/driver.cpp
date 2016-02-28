//*********************************************************
// This file is part of Rapsodia released under the LGPL. *
// The full COPYRIGHT notice can be found in the top      *
// level directory of the Rapsodia distribution           *
//*********************************************************
#include <vector>
#include <stdio.h>
#include <cmath>
#include <cassert>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include "RAinclude.ipp"
#include "HigherOrderTensor.hpp"

#include "../function.h"

int main() { 
  int rc=0;
  double myEps=1.0E-10;
  unsigned short  n,o;
  n=NUM_IND;
  o=DERIVATIVE_ORDER;	
  HigherOrderTensor T(n,o); 
  int dirs=T.getDirectionCount();
  std::cout << "Number of directions: " << dirs << std::endl;

  int i,j,k;
  // argument values
  RAfloatD* x = new RAfloatD[n];
  double* xv = new double[n];

  get_initials(xv, n); 

  for(i=1;i<=n;i++) {
    x[i-1] = xv[i-1];
  }

  // get the seed matrix
  Matrix<unsigned int> SeedMatrix=T.getSeedMatrix();
  for(i=1;i<=dirs;i++) { 
    for(j=1;j<=n;j++) { 
      x[j-1].set(i,1,SeedMatrix[j-1][i-1]);
    }
  }
   // compute the target function
  struct timeval tv1, tv2;
  double time_elapsed = 0;
  gettimeofday(&tv1, NULL);
  RAfloatD y = eval_func<RAfloatD>(x, n);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) + (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  printf("func_ret = %.10f, time elapsed = %.10f\n", y.v, time_elapsed);
  // transfer the taylor coefficients
  Matrix<double> TaylorCoefficients(o,dirs);
  for(i=1;i<=o;i++) { 
    for(j=1;j<=dirs;j++) { 
      TaylorCoefficients[i-1][j-1]=y.get(j,i);
    }
  }
  T.setTaylorCoefficients(TaylorCoefficients);
  // harvest the compressedTensor
  int nnz = 0;
  for (k=1; k<=o; k++){ 
    std::cout << "order: " << k << std::endl;
    double entry=0.0;
    nnz = 0;
    std::vector<double> compressedTensor=T.getCompressedTensor(k);
    HigherOrderTensor Helper(n,k); 
    dirs=Helper.getDirectionCount();
    // get the helper seed matrix
    Matrix<unsigned int> HelperSeedMatrix=Helper.getSeedMatrix();
    for(i=1;i<=dirs;i++) {
      std::vector<unsigned short> index(n);
      if (fabs(compressedTensor[i-1]) < myEps) {
        continue;
      }
      for(j=1;j<=n;j++) { 
	index[j-1]=HelperSeedMatrix[j-1][i-1];
      }
      std::cout << "T";
      for(j=1;j<=n;j++) { 
        std::cout << "[" << std::setw(2) << HelperSeedMatrix[j-1][i-1] << "]";
      }
      if (compressedTensor[i-1] != 0) {
        nnz++;
      }
      std::cout << " = " << compressedTensor[i-1];
      std::cout << std::endl;
    }
    std::cout << "nnz = " << nnz << std::endl;
  }

  printf("func_ret = %.10f, time elapsed = %.10f\n", y.v, time_elapsed);
  gettimeofday(&tv2, NULL);
  time_elapsed = (tv2.tv_sec - tv1.tv_sec) + (double)(tv2.tv_usec - tv1.tv_usec) / 1000000;
  printf("retrieve results, time elapsed = %.10f\n", time_elapsed);
  return rc;
}

