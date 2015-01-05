#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef __APPLE__
  #include <Accelerate/Accelerate.h>
#else
  #include <cblas.h>
#endif
#include "mmio.h"

#define UPPER 121
#define LOWER 122
#define PRINTMATS 0

void zero( int n, double *w) {
  for( int i = 0; i < n; i++) {
    w[i] = 0;
  }
}

int print_vec( const char *name, const int n, const double *x) {
  printf( "Printing %s (len %i):\n", name, n);
  for( int i = 0; i < n; i++) {
    printf("%f\n", x[i]);
  }
  printf("\n");

  return 0;
}

//find sparse y \gets \alpha A x + \beta y
int sparse_matvec( const double alpha, const double beta, const int n, const int nz, //input numbers
          const int *I, const int *J, const double *val, const double *x, //input vectors
          double *y) { //output data
  for( int i = 0; i < n; i++) {
    y[i] = beta * y[i];
  }

  for( int i = 0; i < nz; i++) {
    y[I[i]] = alpha*val[i]*x[J[i]] + y[I[i]];
    if( I[i] != J[i])
      y[J[i]] = alpha*val[i]*x[I[i]] + y[J[i]];
  }

  return 0;
}

double sparse_cg( const int kmax, const double eps, //settings
          const int n, const int nz, const int *I, const int *J, const double *val, const double *b, const double *x0, //input data
          double *x) { //output data
  printf("Solving sparse system Ax=b, dim=%i, nz=%i, kmax=%i, eps=%f, with\n", n, nz, kmax, eps);

  if( PRINTMATS) {
    printf("A = \n");
    for( int i = 0; i < nz; i++) {
      printf("A[%i][%i] = %.4f\n", I[i], J[i], val[i]);
      if( I[i] != J[i])
        printf("A[%i][%i] = %.4f\n", J[i], I[i], val[i]);
    }
    printf("\n");
  }

  int k = 0;
  double r[n], rho, rho_old, normbsq;

  cblas_dcopy( n, x0, 1, x, 1); // x \gets x_0
  cblas_dcopy( n, b, 1, r, 1); // r \gets b
  sparse_matvec( -1.0, 1.0, n, nz, I, J, val, x, r); // r \gets r - Ax
  rho = cblas_ddot( n, r, 1, r, 1); //rho \gets <r, r>
  printf("rho = %g\n", rho);
  normbsq = cblas_ddot( n, b, 1, b, 1); //find <b, b>

  while( rho > eps*eps*normbsq && k < kmax) {
    double p[n], w[n], alpha, beta, gamma;

    if( k == 0) {
      cblas_dcopy( n, r, 1, p, 1); //p \gets r
    } else {
      beta = rho/rho_old;
      cblas_dscal( n, beta, p, 1); // p \gets \beta p
      cblas_daxpy( n, 1.0, r, 1, p, 1); // p \gets r + p
    }

    sparse_matvec( 1.0, 0.0, n, nz, I, J, val, p, w);
    gamma = cblas_ddot( n, p, 1, w, 1); // \gamma \gets <p, w>
    alpha = rho/gamma;
    cblas_daxpy( n, alpha, p, 1, x, 1); // x \gets x + \alpha p
    cblas_daxpy( n, - alpha, w, 1, r, 1); // r \gets r - \alpha w;
    rho_old = rho;
    rho = cblas_ddot( n, r, 1, r, 1); // \rho \gets <r, r>
    printf("rho = %g\n", rho);

    k = k+1;
  }

  printf("Solution found after %i iterations! rho = %f\n", k, rho);
  print_vec("x", n, x);

  double b2[n];
  sparse_matvec( 1.0, 0.0, n, nz, I, J, val, x, b2); // b2 \gets Ax
  print_vec("Ax", n, b2);
  return rho;
} 

double dense_cg( const int kmax, const double eps, //settings
          const int n, const double *A, const double *b, const double *x0, //input data
          double *x) { //output data
  printf("Solving system Ax=b, dim=%i, kmax=%i, eps=%f, with\n", n, kmax, eps);

  if( PRINTMATS) {
    printf("A = \n");
    for( int i = 0; i < n; i++) {
      for( int j = 0; j < n; j++) {
        printf("%.4f ", A[i*n + j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  int k = 0;
  double r[n], rho, rho_old, normbsq;

  cblas_dcopy( n, x0, 1, x, 1); // x \gets x_0
  cblas_dcopy( n, b, 1, r, 1); // r \gets b
  cblas_dsymv( CblasRowMajor, LOWER, n, -1.0, A, n, x, 1, 1.0, r, 1); // r \gets r - Ax
  rho = cblas_ddot( n, r, 1, r, 1); //rho \gets <r, r>
  printf("rho = %g\n", rho);
  normbsq = cblas_ddot( n, b, 1, b, 1); //find <b, b>

  while( rho > eps*eps*normbsq && k < kmax) {
    double p[n], w[n], alpha, beta, gamma;

    if( k == 0) {
      cblas_dcopy( n, r, 1, p, 1); //p \gets r
    } else {
      beta = rho/rho_old;
      cblas_dscal( n, beta, p, 1); // p \gets \beta p
      cblas_daxpy( n, 1.0, r, 1, p, 1); // p \gets r + p
    }

    cblas_dsymv( CblasRowMajor, LOWER, n, 1.0, A, n, p, 1, 0.0, w, 1); // w \gets Ap
    gamma = cblas_ddot( n, p, 1, w, 1); // \gamma \gets <p, w>
    alpha = rho/gamma;
    cblas_daxpy( n, alpha, p, 1, x, 1); // x \gets x + \alpha p
    cblas_daxpy( n, - alpha, w, 1, r, 1); // r \gets r - \alpha w;
    rho_old = rho;
    rho = cblas_ddot( n, r, 1, r, 1); // \rho \gets <r, r>
    printf("rho = %g\n", rho);

    k = k+1;
  }

  printf("Solution found after %i iterations! rho = %f\n", k, rho);
  print_vec("x", n, x);

  double b2[n];
  cblas_dsymv( CblasRowMajor, LOWER, n, 1.0, A, n, x, 1, 0.0, b2, 1); //b2 \gets Ax
  print_vec("Ax", n, b2);
  return rho;
}

int main( int argc, char *argv[]) {
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;
  int i, *I, *J;
  double *val;
  double *A, *b, *x0, *x;

  f = fopen("../mtxMatrices/LFAT5.mtx", "r");
  if( mm_read_banner( f, &matcode) != 0) {
    printf("Could not process Matrix Market banner.\n");
    exit(1);
  }

  if( mm_is_complex( matcode) && mm_is_matrix( matcode) && mm_is_sparse( matcode)) {
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    exit(1);
  }

  if( (ret_code = mm_read_mtx_crd_size( f, &M, &N, &nz)) != 0) {
    exit(1);
  }

  if( M != N) {
    printf("Square matrix required!\n");
    exit(1);
  }

  I = (int *) malloc( nz * sizeof( int));
  J = (int *) malloc( nz * sizeof( int));
  val = (double *) malloc( nz * sizeof( double));

  for( i = 0; i < nz; i++) {
    fscanf( f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
    I[i]--;
    J[i]--;
  }

  fclose( f);

  A = (double *) malloc( M*N*sizeof( double));
  b = (double *) malloc( M * sizeof( double));
  x0 = (double *) malloc( M * sizeof( double));
  x = (double *) malloc( M * sizeof( double));

  for( i = 0; i < M; i++) {
    b[i] = 1;
    x0[i] = 0;
    x[i] = 0;
  }

  for( i = 0; i < M*N; i++) {
    A[i] = 0;
  }

  for( i = 0; i < nz; i++) {
    A[I[i]*M + J[i]] = val[i];
    A[J[i]*M + I[i]] = val[i];
  }

  /*
  M = 2;
  N = 2;
  nz = 3;
  I = (int *) malloc( nz * sizeof( int));
  J = (int *) malloc( nz * sizeof( int));
  val = (double *) malloc( nz * sizeof( double));
  A = (double *) malloc( M*N*sizeof( double));
  b = (double *) malloc( M * sizeof( double));
  x0 = (double *) malloc( M * sizeof( double));
  x = (double *) malloc( M * sizeof( double));

  A[0] = 4.0;
  A[1] = 1.0;
  A[2] = 1.0;
  A[3] = 3.0;
  I[0] = 0;
  J[0] = 0;
  val[0] = 4.0;
  I[1] = 0;
  J[1] = 1;
  val[1] = 1.0;
  I[2] = 1;
  J[2] = 1;
  val[2] = 3.0;
  b[0] = 1;
  b[1] = 2;
  x0[0] = 2;
  x0[1] = 1;
  x[0] = 0;
  x[1] = 0;
  */

  dense_cg( 14, 0.00001, M, A, b, x0, x);
  sparse_cg( 14, 0.00001, M, nz, I, J, val, b, x0, x);

  return 0;
}
