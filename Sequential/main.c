#include <stdio.h>
#include <math.h>
#include <Accelerate/Accelerate.h>

float cg( const int kmax, const float eps, //settings
          const int n, const float A[n][n], const float b[n], const float x0[n], //input data
          float x[n]) { //output data
  printf("Solving system Ax=b, dim=%i, kmax=%i, eps=%f, with\n", n, kmax, eps);

  printf("A = \n");
  for( int i = 0; i < n; i++) {
    for( int j = 0; j < n; j++) {
      printf("%.4f ", A[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("b = \n");
  for( int i = 0; i < n; i++) {
    printf("%.4f\n", b[i]);
  }
  printf("\n");

  printf("x0 = \n");
  for( int i = 0; i < n; i++) {
    printf("%.4f\n", x0[i]);
  }
  printf("\n");
  int k = 0;
  float r[n];
  float rho;
  float rho_old;
  float normbsq;

  cblas_scopy( n, x0, 1, x, 1); // x \gets x_0
  cblas_scopy( n, b, 1, r, 1); // r \gets b
  cblas_ssymv( CblasRowMajor, 122, n, -1.0, (float *)A, n, x, 1, 1.0, r, 1); // r \gets r - Ax
  rho = cblas_sdot( n, r, 1, r, 1); //rho \gets <r, r>
  normbsq = cblas_sdot( n, b, 1, b, 1); //find <b, b>



  while( rho > eps*eps*normbsq && k < kmax) {
    float p[n];
    float w[n];
    float alpha;
    float beta;
    float gamma;

    if( k == 0) {
      cblas_scopy( n, r, 1, p, 1); //p \gets r
    } else {
      beta = rho/rho_old;
      cblas_sscal( n, beta, p, 1); // p \gets \beta p
      cblas_saxpy( n, 1.0, r, 1, p, 1); // p \gets r + p
    }

    cblas_ssymv( CblasRowMajor, 122, n, 1.0, (float *)A, n, p, 1, 0.0, w, 1); // w \gets Ap
    gamma = cblas_sdot( n, p, 1, w, 1); // \gamma \gets <p, w>
    alpha = rho/gamma;
    cblas_saxpy( n, alpha, p, 1, x, 1); // x \gets x + \alpha p
    cblas_saxpy( n, - alpha, w, 1, r, 1); // r \gets r - \alpha w;
    rho_old = rho;
    rho = cblas_sdot( n, r, 1, r, 1); // \rho \gets <r, r>

    k = k+1;
  }

  printf("Solution found! rho = %f\n", rho);
  printf("x = \n");
  for( int i = 0; i < n; i++) {
    printf("%.4f\n", x[i]);
  }
  printf("\n");

  float b2[n];
  cblas_ssymv( CblasRowMajor, 122, n, 1.0, (float *)A, n, x, 1, 0.0, b2, 1); //b2 \gets Ax

  printf("Ax = \n");
  for( int i = 0; i < n; i++) {
    printf("%.4f\n", b2[i]);
  }
  printf("\n");
  return rho;
}

int main() {
  float x0[2] = {2.0, 1.0};
  float a[2][2] = {{2.0, 1.0}, {1.0, 3.0}};
  float b[2] = {1.0, 2.0};
  float x[2] = {0.0, 0.0};

  float end_roundoff = cg( 2, 0.001, 2, a, b, x0, x);
  return 0;
}
