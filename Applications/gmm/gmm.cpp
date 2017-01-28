#include <iostream>
#include <string>
#include <chrono>
#include <set>

#include "utils.h"
#include "defs.h"

#include "gmm.h"



using std::cout;
using std::endl;
using std::string;
using namespace std::chrono;

#ifdef USING_REVERSEAD
#include "reversead/reversead.hpp"
using namespace ReverseAD;
#endif

#ifdef USING_ADOLC
#include "adolc/adolc.h"
#endif

/*
double compute_gmm_J(int nruns,
	int d, int k, int n, double *alphas, 
	double *means, double *icf, double *x, 
	Wishart wishart, double& err, double *J)
{
	int tapeTag = 1;
	int icf_sz = d*(d + 1) / 2;
	int Jrows = 1;
	int Jcols = (k*(d + 1)*(d + 2)) / 2;
	adouble *aalphas, *ameans, *aicf, aerr;

	// Record on a tape
	trace_on(tapeTag);

	aalphas = new adouble[k];
	for (int i = 0; i < k; i++)
	{
		aalphas[i] <<= alphas[i];
	}
	ameans = new adouble[d*k];
	for (int i = 0; i < d*k; i++)
	{
		ameans[i] <<= means[i];
	}
	aicf = new adouble[icf_sz*k];
	for (int i = 0; i < icf_sz*k; i++)
	{
		aicf[i] <<= icf[i];
	}

	gmm_objective(d, k, n, aalphas, ameans, 
		aicf, x, wishart, &aerr);

	aerr >>= err;

	trace_off();

	delete[] aalphas;
	delete[] ameans;
	delete[] aicf;


	high_resolution_clock::time_point start, end;
	start = high_resolution_clock::now();

	// Compute J
	double *in = new double[Jcols];
	memcpy(in, alphas, k*sizeof(double));
	int off = k;
	memcpy(in + off, means, d*k*sizeof(double));
	off += d*k;
	memcpy(in + off, icf, icf_sz*k*sizeof(double));

	for (int i = 0; i < nruns; i++)
	{
		gradient(tapeTag, Jcols, in, J);

		//int keepValues = 1;
		//double errd = 1;
		//zos_forward(tapeTag, Jrows, Jcols, keepValues, in, &err);
		//fos_reverse(tapeTag, Jrows, Jcols, &errd, J[0]);
	}

	end = high_resolution_clock::now();

	delete[] in;

  return duration_cast<duration<double>>(end - start).count() / nruns;
}
*/

double compute_gmm_H(int nruns,
	int d, int k, int n, double *alphas, 
	double *means, double *icf, double *x, 
	Wishart wishart, double& err)
{
	int icf_sz = d*(d + 1) / 2;
	int num_ind = (k*(d + 1)*(d + 2)) / 2;
	adouble *aalphas, *ameans, *aicf, aerr;

	high_resolution_clock::time_point start, end;
	start = high_resolution_clock::now();

	// Record on a tape
#ifdef USING_REVERSEAD
	trace_on<double>();
#endif
#ifdef USING_ADOLC
        trace_on(1);
#endif

	aalphas = new adouble[k];
	for (int i = 0; i < k; i++)
	{
		aalphas[i] <<= alphas[i];
	}
	ameans = new adouble[d*k];
	for (int i = 0; i < d*k; i++)
	{
		ameans[i] <<= means[i];
	}
	aicf = new adouble[icf_sz*k];
	for (int i = 0; i < icf_sz*k; i++)
	{
		aicf[i] <<= icf[i];
	}

	gmm_objective(d, k, n, aalphas, ameans, 
		aicf, x, wishart, &aerr);

	aerr >>= err;

#ifdef USING_REVERSEAD
	std::shared_ptr<TrivialTrace<double>> trace = trace_off<double>();
#endif

#ifdef USING_ADOLC
        trace_off();
#endif

	end = high_resolution_clock::now();
        printf("overload function tracing time = %.6f\n", 
                 duration_cast<duration<double>>(end-start).count());
	delete[] aalphas;
	delete[] ameans;
	delete[] aicf;




#ifdef USING_ADOLC
#ifdef ADOLC_FULLHESS
        double** H = new double*[num_ind];
        for (int i = 0; i < num_ind; i++) {
          H[i] = new double[num_ind];
        }
#else
        double* v = new double[num_ind];
        double* z = new double[num_ind];
        for (int i = 0; i < num_ind; i++) {v[i] = 0;}
#endif
        double *in = new double[num_ind];
	memcpy(in, alphas, k*sizeof(double));
	int off = k;
	memcpy(in + off, means, d*k*sizeof(double));
	off += d*k;
	memcpy(in + off, icf, icf_sz*k*sizeof(double));
#endif
	// Compute H
        std::cout << "num_ind = " << num_ind << std::endl;
	start = high_resolution_clock::now();
	for (int i = 0; i < nruns; i++)
	{
#ifdef USING_REVERSEAD
          BaseReverseHessian<double> hessian(trace);
          std::shared_ptr<DerivativeTensor<size_t, double>> tensor = hessian.compute(num_ind, 1);
          size_t size;
          size_t** tind;
          double* values;
          tensor->get_internal_coordinate_list(0, 2, &size, &tind, &values);
          std::cout << i << " : " << " hessian size = " << size;
#endif
#ifdef USING_ADOLC
#ifdef ADOLC_FULLHESS
          hessian(1, num_ind, in, H);
#else
          v[i] = 1;
          hess_vec(1, num_ind, in, v, z);
#endif

#endif
	  end = high_resolution_clock::now();
          std::cout << "accumute time = "
            << duration_cast<duration<double>>(end-start).count() << std::endl;


	}

	end = high_resolution_clock::now();
  return duration_cast<duration<double>>(end - start).count() / nruns;
}

void test_gmm(const string& fn_in, int nruns_f, int nruns_H, bool replicate_point)
{
  int d, k, n;
  vector<double> alphas, means, icf, x;
  double err;
  Wishart wishart;

  // Read instance
  read_gmm_instance(fn_in, &d, &k, &n,
    alphas, means, icf, x, wishart, replicate_point);

  int num_ind = (k*(d + 1)*(d + 2)) / 2;


  // Test
  high_resolution_clock::time_point start, end;
  double tf, tH = 0.;

  start = high_resolution_clock::now();
  for (int i = 0; i < nruns_f; i++)
  {
    gmm_objective(d, k, n, alphas.data(), means.data(),
      icf.data(), x.data(), wishart, &err);
  }
  end = high_resolution_clock::now();
  tf = duration_cast<duration<double>>(end - start).count() / nruns_f;
  cout << "err: " << err << endl;

  string name = "ADOLC";
  tH = compute_gmm_H(nruns_H, d, k, n, alphas.data(), means.data(), icf.data(), 
    x.data(), wishart, err);

  cout << "err: " << err << endl;

  printf("Plain objective function time = %.6f\n", tf);
#ifdef USING_REVERSEAD
  printf("ReverseAD Hessian time = %.6f\n", tH);
#endif

#ifdef USING_ADOLC
#ifdef ADOLC_FULLHESS
  printf("ADOLC fullhess time = %.6f\n", tH);
#else
  printf("ADOLC hess_vec time = %.6f\n", tH);
#endif
#endif

}


int main(int argc, char *argv[])
{
  string fn(argv[1]);
  int nruns_f = std::stoi(string(argv[2]));
  int nruns_H = std::stoi(string(argv[3]));

  // read only 1 point and replicate it?
  bool replicate_point = (argc >= 5 && string(argv[4]).compare("-rep") == 0);
  
  test_gmm(fn, nruns_f, nruns_H, replicate_point);
}
