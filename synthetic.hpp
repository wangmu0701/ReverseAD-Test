#include <cmath>

#define PRINT_RESULTS

#define DERIVATIVE_ORDER 5

#define NUM_IND 30
#define LIVE_SIZE 15
#define NUM_DOUBLE_LAYERS 5

void get_initials(double* x, int n) {
  for (int i = 0; i < n; i++) {
    x[i] = 1.0 / (i + 2.0);
  }
}

template <typename T>
T eval_func(T* xad, int n) {
  T fad;
  T z[LIVE_SIZE];
  for (int i = 0; i < LIVE_SIZE; i++) {
    z[i] = i;
    for (int j = 0; j < NUM_IND; j++) {
      z[i] += xad[j]; 
    }
  }
  for (int i = 0; i < NUM_DOUBLE_LAYERS; i++) {
    for (int j = 0; j < LIVE_SIZE; j++) {
      z[j] = z[j] * z[j];
    }
    for (int j = 0; j < LIVE_SIZE; j++) {
      z[j] = sqrt(z[j]);
    }
  }
  fad = 1.0;
  for (int i = 0; i < LIVE_SIZE; i++) {
    fad = fad * z[i];
  }
  return fad;
}
