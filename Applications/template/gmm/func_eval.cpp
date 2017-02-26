#include <vector>
#include <cassert>
#include <iostream>
#include <string>
#include <fstream>
#include <string.h>
#include "defs.h"

#include "func_eval.hpp"


using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::getline;
using std::vector;

int d = 0;
int k = 0;
int nn;
Wishart wishart;
vector<double> x;

void read_gmm_instance(const string& fn,
                       int *d, int *k, int *n,
                       vector<double>& alphas,
                       vector<double>& means,
                       vector<double>& icf,
                       vector<double>& x,
                       Wishart& wishart,
                       bool replicate_point)
{
    FILE *fid = fopen(fn.c_str(), "r");
    
    fscanf(fid, "%i %i %i", d, k, n);
    
    int d_ = *d, k_ = *k, n_ = *n;
    
    int icf_sz = d_*(d_ + 1) / 2;
    alphas.resize(k_);
    means.resize(d_*k_);
    icf.resize(icf_sz*k_);
    x.resize(d_*n_);
    
    for (int i = 0; i < k_; i++)
    {
        fscanf(fid, "%lf", &alphas[i]);
    }
    
    for (int i = 0; i < k_; i++)
    {
        for (int j = 0; j < d_; j++)
        {
            fscanf(fid, "%lf", &means[i*d_ + j]);
        }
    }
    
    for (int i = 0; i < k_; i++)
    {
        for (int j = 0; j < icf_sz; j++)
        {
            fscanf(fid, "%lf", &icf[i*icf_sz + j]);
        }
    }
    
    if (replicate_point)
    {
        for (int j = 0; j < d_; j++)
        {
            fscanf(fid, "%lf", &x[j]);
        }
        for (int i = 0; i < n_; i++)
        {
            memcpy(&x[i*d_], &x[0], d_ * sizeof(double));
        }
    }
    else
    {
        for (int i = 0; i < n_; i++)
        {
            for (int j = 0; j < d_; j++)
            {
                fscanf(fid, "%lf", &x[i*d_ + j]);
            }
        }
    }
    
    fscanf(fid, "%lf %i", &(wishart.gamma), &(wishart.m));
    
    fclose(fid);
}

void set_up(int argc, char** argv, vector<double>& ind, int& n, int& m) {
  string fn_in(argv[1]);
  vector<double> alphas, means, icf;

  read_gmm_instance(fn_in, &d, &k, &nn,
      alphas, means, icf, x, wishart, false);
  n = (k*(d+1)*(d+2)) / 2;
  m = 1;

  ind.clear();
  ind.reserve(n);
  int icf_sz = d * (d+1) / 2;
  for (int i = 0; i < k; i++) {
    ind.push_back(alphas[i]);
  }
  for (int i = 0; i < d*k; i++) {
    ind.push_back(means[i]);
  }
  for (int i = 0; i < icf_sz*k; i++) {
    ind.push_back(icf[i]);
  }
  assert(ind.size() == n);
}




void tear_down() { // do nothing?
}

