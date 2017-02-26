#include <cmath>
#include <cassert>
#include <memory>
#include <iostream>
#include <vector>

#include "mesh.hpp"

////////////////////////////////////////////////////////////
//////////////////// Declarations //////////////////////////
////////////////////////////////////////////////////////////

using std::vector;

extern Mesh *grid;
extern double sqrt3;
extern double sqrt6;
extern double a;
extern double b;

void set_up(int argc, char* argv[], vector<double>& ind_value, int& n, int& m);
void tear_down();

#define rcbrt(x) pow(x,-3.333333333333333333333333333e-01)
#define epsilonf 1.0e-14

template <typename T>
bool a_fcn(T& obj, const T x[12])
{
  static T matr[9], f;
  static T g;
  f       = x[1] + x[0];
  matr[0] = x[1] - x[0];
  matr[1] = (2.0*x[2] - f)*sqrt3;
  matr[2] = (3.0*x[3] - x[2] - f)*sqrt6;
 f       = x[5] + x[4];
  matr[3] = x[5] - x[4];
  matr[4] = (2.0*x[6] - f)*sqrt3;
  matr[5] = (3.0*x[7] - x[6] - f)*sqrt6;

  f       = x[9] + x[8];
  matr[6] = x[9] - x[8];
  matr[7] = (2.0*x[10] - f)*sqrt3;
  matr[8] = (3.0*x[11] - x[10] - f)*sqrt6;

  g = matr[0]*(matr[4]*matr[8] - matr[5]*matr[7]) +
      matr[1]*(matr[5]*matr[6] - matr[3]*matr[8]) +
      matr[2]*(matr[3]*matr[7] - matr[4]*matr[6]);

  // This line is a problem for libtaylor, remove it because in valid inputs
  // this line never get invoked.
  //if (g <= epsilonf) { obj = g; return true; }

  f = matr[0]*matr[0] + matr[1]*matr[1] + matr[2]*matr[2] +
      matr[3]*matr[3] + matr[4]*matr[4] + matr[5]*matr[5] +
      matr[6]*matr[6] + matr[7]*matr[7] + matr[8]*matr[8];

  obj = a * f * rcbrt(g*g);
  return false;
}

template <typename T>
T aFcn(T *v, const int e, const int *t)
{
  //const int e = m->ne;
  T *w;
  //int    *t = m->e;

  T  x[12];
  T  o, f;
  int     v1, v2, v3, v4;
  int     i;

  T obj = 0.0;

  o = 0.0;
  for (i = 0; i < e; ++i) {
    v1 = t[0];
    v2 = t[1];
    v3 = t[2];
    v4 = t[3];
    t += 4;

    w = v + 3*v1;
    x[0] = w[0];
    x[4] = w[1];
    x[8] = w[2];

    w = v + 3*v2;
    x[1] = w[0];
    x[5] = w[1];
    x[9] = w[2];

    w = v + 3*v3;
    x[2] = w[0];
    x[6] = w[1];
    x[10]= w[2];

    w = v + 3*v4;
    x[3] = w[0];
    x[7] = w[1];
    x[11]= w[2];
    if (a_fcn(f, x)) {
      std::cerr << "something is wrong, check the input mesh!" << std::endl;
      exit(-1);
    };

    o += f;
  }
  obj = o;
  return obj;
}

template<typename T>
void func_eval(int n, T* ind_val, int m, T* dep_val) {
  assert(m == 1);
  T y = aFcn<T>(ind_val, grid->ne, grid->e);
  dep_val[0] = y;
}

