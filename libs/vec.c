#include <stdio.h>
#include <stdlib.h>
#include "vec.h"

int vec_print( const char *name, const int n, const double *x) {
  printf( "Printing %s (len %i):\n", name, n);
  for( int i = 0; i < n; i++) {
    printf("%f\n", x[i]);
  }
  printf("\n");

  return 0;
}

double rand_double(void) {
  return ( (double) rand() / (double) RAND_MAX);
}

//p \gets random
void vec_random( const int n, const int seed, double *p) {
  srand( seed);
  for( int i = 0; i < n; i++) {
    p[i] = rand_double();
  }
}

//p \gets r
void vec_copy( const int n, const double *r, double *p) {
  for( int i = 0; i < n; i++) {
    p[i] = r[i];
  }
}

//p \gets \beta p
void vec_scale( const int n, const double beta, double *p) {
  for( int i = 0; i < n; i++) {
    p[i] *= beta;
  }
}

//r \gets r + \alpha w
void vec_axpy( const int n, const double alpha, const double *w, double *r) {
  for( int i = 0; i < n; i++) {
    r[i] += alpha * w[i];
  }
}
