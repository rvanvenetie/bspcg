#include <stdio.h>
#include <math.h>
#include <Accelerate/Accelerate.h>

float cg( const int kmax, const float eps, const int n, const float **A, const float *b, const float *x0, float *x) {
  printf("Solving system Ax=b, dim=%i, kmax=%i, eps=%f, with\n", n, kmax, eps);
  printf("A = \n");
  for( int i = 0; i < n; i++) {
    for( int j = 0; i < n; i++) {
      printf("%.4f", A[i][j]);
    }
  }
  int k = 0;
  float r[n];
  float rho;
  float rho_old;
  float normbsq;

  cblas_scopy( n, b, 1, r, 1); // r \to b
  rho = cblas_sdot( n, r, 1, r, 1); //rho \to <r, r>
  normbsq = cblas_sdot( n, b, 1, b, 1); //find <b, b>
  cblas_scopy( n, x0, 1, x, 1); // x \to x_0


  while( rho > eps*eps*normbsq && k < kmax) {
    float p[n];
    float w[n];
    float alpha;
    float beta;
    float gamma;

    if( k == 0) {
      cblas_scopy( n, r, 1, p, 1); //p \to r
    } else {
      beta = rho/rho_old;
      cblas_sscal( n, beta, p, 1); // p \to \beta p
      cblas_saxpy( n, 1.0, r, 1, p, 1); // p \to r + p
    }

    cblas_ssymv( CblasRowMajor, 122, n, 1.0, (float *)A, n, p, 1, 0.0, w, 1); // w \to Ap
    gamma = cblas_sdot( n, p, 1, w, 1); // \gamma \to <p, w>
    alpha = rho/gamma;
    cblas_saxpy( n, alpha, p, 1, x, 1); // x \to x + \alpha p
    cblas_saxpy( n, - alpha, w, 1, r, 1); // r \to r - \alpha w;
    rho_old = rho;
    rho = cblas_sdot( n, r, 1, r, 1); // \rho \to <r, r>

    k = k+1;
  }

  return rho;
}

int main() {
  float x0[2] = {1.0, 2.0};
  float a[2][2] = {{4.0, 1.0}, {1.0, 3.0}};
  float b[2] = {1.0, 2.0};
  float x[2] = {0.0, 0.0};

  float end_roundoff = cg( 2, 0.001, 2, a, b, x0, x);

  printf("%f\n", x[0]);
  printf("%f\n", x[1]);
  printf("%f\n", end_roundoff);
  return 0;
}
