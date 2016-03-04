#include <cmath>
#include <cstdlib>

//#define PRINT_RESULTS

#define DERIVATIVE_ORDER 3


// The number of SAC of the code is given by :
// NUM_IND + LIVE_SIZE * (NUM_DOUBLE_LAYERS + 1) * 2
// We decide to make LIVE_SIZE * (NUM_DOUBLE_LAYERS + 1) = 4900
// So the unbalanced SAC number is only 1%.
// [2,244], [3, 162], [4, 122], [5, 97], [6, 81], [7, 69], [8, 61], [9,53]
#define NUM_IND 20
#define LIVE_SIZE NUM_IND
#define NUM_DOUBLE_LAYERS 80

#define RANDOM_SEED 12345

void get_initials(double* x, int n) {
  for (int i = 0; i < n; i++) {
    x[i] = 1.0 - 1.0 / (i + 2.0);
  }
  srand(RANDOM_SEED);
}


template <typename T>
T eval_random_identity(T& rhs) {
  int M = 6;
  int l = rand() % M;
  switch (l) {
    case 0:
      return sqrt(rhs*rhs);
    case 1:
      return 2.0 + rhs - 2.0;
    case 2:
      return rhs * 2.0 * 0.5;
    case 3:
      return log(exp(rhs));
    case 4:
      return 1.0 / (1.0 / rhs);
    case 5:
      return sin(asin(rhs));
    default:
      return rhs;
  }
  return rhs;
}

template <typename T>
T eval_func(T* xad, int n) {
  T fad;
  T z[LIVE_SIZE];
  T t = 1.0;
  for (int i = 0; i < NUM_IND; i++) {
    t = t * xad[i];
  }
  for (int i = 0; i < LIVE_SIZE; i++) {
    z[i] = t;
  }

  for (int i = 0; i < NUM_DOUBLE_LAYERS; i++) {
    for (int j = 0; j < LIVE_SIZE; j++) {
      z[j] = eval_random_identity(z[j]);
    }
  }
  fad = 1.0;
  for (int i = 0; i < LIVE_SIZE; i++) {
    fad = fad * z[i];
  }
  return sqrt(fad);
}

