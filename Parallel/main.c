#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cblas.h>
#include "io.h"
#include "bspmv.h"
#include "bspedupack.h"


int print_vec( const char *name, const int n, const double *x) {
  printf( "Printing %s (len %i):\n", name, n);
  for( int i = 0; i < n; i++) {
    printf("%f\n", x[i]);
  }
  printf("\n");

  return 0;
}
//p \gets r
void copy_vec( const int n, const double *r, double *p) {
  for( int i = 0; i < n; i++) {
    p[i] = r[i];
  }
}

//p \gets \beta p
void scale_vec( const int n, const double beta, double *p) {
  for( int i = 0; i < n; i++) {
    p[i] *= beta;
  }
}

//r \gets r + \alpha w
void axpy_vec( const int n, const double alpha, const double *w, double *r) {
  for( int i = 0; i < n; i++) {
    r[i] += alpha * w[i];
  }
}
/*
 * Edited version of bspip that uses vectors that are stored
 * according to a certain distribution (both must have the same).

 * Input:
 * p - number of processors.
 * s - current processor number.
 * x - locally stored vector elements of vector x
 * y - locally stored vector elements of vector y
 * nv - amount of vector elements locally stored.
 */
double bspip_dist(int p, int s,  //Processor information  
    double * x, double * y,      //Vector data
    int nv)                      //Distribution information
{
  /*
   * Every processor calculates its own local inner product.
   * It then puts this value in all the other processors.
   * After which every processor calculates the total inner
   * product
   */
  //Initialize array that holds all partial sums
  double inprod_local;
  double * inprod;
  inprod = vecallocd(p);
  bsp_push_reg(inprod, p * SZDBL);
  bsp_sync();

  //Calculate local inproduct
  inprod_local = 0;
  for (int i = 0; i < nv; i++)
    inprod_local += x[i] * y[i];

  //Put the local value in all the processors
  for (int t = 0; t < p; t++) 
    bsp_put(t, &inprod_local, inprod, s * SZDBL, SZDBL);

  bsp_sync();
  //Calculate the actual inner product
  inprod_local = 0;
  for (int t = 0; t < p; t++)
    inprod_local += inprod[t];

  //Free data stuff
  bsp_pop_reg(inprod);
  vecfreed(inprod);

  return inprod_local;
} 
int P;
#define kmax 200
#define eps 0.01
#define mat_file  "../mtxMatrices/bodyy5.mtx-P%d"
#define dist_file "../mtxMatrices/bodyy5.mtx-u%d"
#define b_file    "../mtxMatrices/bodyy5.mtx-b"

#define DEBUG(...) { if (s == 0) printf(__VA_ARGS__); }

#define PRINTMATS 1
void bspcg() {
  int p,s; //Defaults for bsp
  //Vectors
  double *b; //Holds rhs  
  double *r; //Holds the residual
  double *x; //Holds the approximate solution
  double *u; //Used in the algorithm to update x
  double *w; //Used in the algorithm to calculate alpha
  
  bsp_begin(P);
  p= bsp_nprocs(); /* p = number of processors obtained */
  s= bsp_pid();    /* s = processor number */
  
  /* Read matrix */
  char matbuffer[1024], vecdisbuffer[1024], vecvalbuffer[1024];
  snprintf( matbuffer, 1024, mat_file, p);
  snprintf( vecdisbuffer, 1024, dist_file, p);
  snprintf( vecvalbuffer, 1024, b_file);

  distributed_matrix mat = load_symm_distributed_matrix_from_file( matbuffer, p, s);
  vector_distribution dis = load_vector_distribution_from_file( vecdisbuffer, vecvalbuffer, p, s, &b);
  bsp_sync();
  if (s == 0) {
    printf("Solving sparse system Ax=b, dim=%i, nz=%i, kmax=%i, eps=%f, with\n", mat.n, mat.nz, kmax, eps);
    printf("Using %d processors\n",p);
  }
  /* Matrix, vector dist and b is loaded */
  /* Loading phase done, bsp_sync, add timer?*/

  /* Initalize bsp_mv */
  //This are arrays needed for bspmv
  int *srcprocv, *srcindv, *destprocu, *destindu;

  srcprocv  = vecalloci(mat.ncols);
  srcindv   = vecalloci(mat.ncols);
  destprocu = vecalloci(mat.nrows);
  destindu  = vecalloci(mat.nrows);

  bspmv_init(p, s, mat.n, mat.nrows, mat.ncols,
             dis.nv, dis.nv, mat.rowindex, mat.colindex,
	     dis.vindex, dis.vindex, 
	     srcprocv, srcindv, destprocu, destindu);

  int k = 0; //Iterations
  //Initalize vectors used in algorithm, b is already initialized
  r = vecallocd(dis.nv); 
  u = vecallocd(dis.nv);
  x = vecallocd(dis.nv);
  w = vecallocd(dis.nv);
  //u and w are initialized in the algorithm.

  //Initial value of x = 0, r = b - Ax = b.
  for (int i = 0; i < dis.nv; i++) {
    x[i] = 0;
    r[i] = b[i];
  }

  //Rho used in algo, initial rho = <r,r> = <b,b>
  double rho, rho_old;
  rho = bspip_dist(p,s,r,r, dis.nv);

  //We store |b|^2 to calculate the relative norm of r
  //|b|^2 = <b,b> = rho
  double normbsq = rho;
  //bsp_sync?
  while( rho > eps*eps*normbsq && k < kmax) {
    double beta, gamma, alpha;
    DEBUG("Iteration %d, rho = %g\n", k + 1, rho);
    
    if( k == 0) {
      copy_vec(dis.nv, r, u); //u \gets r
    } else {
      beta = rho/rho_old;
      //u \gets r + beta u
      scale_vec(dis.nv, beta, u); // u \gets \beta p
      axpy_vec(dis.nv, 1.0, r, u); // u \gets r + p
    }
    //print_vec("u", dis.nv, u);
    //print_vec("A", mat.nz, mat.val);
    // w \gets Au
    bspmv(p,s, mat.n, mat.nz, mat.nrows, mat.ncols, 
	  mat.val, mat.inc,
	  srcprocv, srcindv, destprocu, destindu,
	  dis.nv, dis.nv, u, w);
    //print_vec("w", dis.nv, w);
    // \gamma \gets <u, w>
    gamma = bspip_dist(p,s, u,w, dis.nv);
    alpha = rho/gamma;
    //  x \gets x + alpha u
    axpy_vec(dis.nv, alpha, u, x); 
    //  r \gets r - alpha w
    axpy_vec(dis.nv, - alpha, w, r); 
    rho_old = rho;
    // \rho \gets <r ,r>
    rho = bspip_dist(p,s, r,r, dis.nv);
    k++;
    //print_vec("x", dis.nv, x);
    //print_vec("r", dis.nv, r);
  }

  if (k < kmax) {
    DEBUG("Solution found after %i iterations! rho = %f\n", k, rho);
  } else {
    DEBUG("CG stopped, maximum iterations (%i) reached. rho = %g\n", k, rho);
  }
  //bsp_sync?
  //Calculate w = Ax
  //Print timing and resulting vector w/b 

  /*
   * We now give processor zero the solution vector.
   * Note that we probably do not want to do this in the real
   * application, as this vector might be to big to store
   * on one processor
   */
  
  double * x_glob;
  x_glob = vecallocd(dis.n);

  bsp_push_reg(x_glob, dis.n*SZDBL);
  bsp_sync();
  for (int i = 0; i < dis.nv; i++)
  {
    int index_glob = dis.vindex[i];
    bsp_put(0, &x[i], x_glob, index_glob * SZDBL, SZDBL);
  }

  bsp_sync();
  bsp_pop_reg(x_glob);
  /* 
  if (s == 0) {
    printf("Solution vector:\n");
    for (int i = 0; i <dis.n; i++)
      printf("%g\n", x_glob[i]);
  }
  */
  vecfreed(x_glob);

  vecfreed(u);
  vecfreed(x);
  vecfreed(r);
  vecfreed(w);
  vecfreed(b);
  vecfreei(srcprocv);
  vecfreei(srcindv);
  vecfreei(destprocu);
  vecfreei(destindu);
  bsp_end();
}


int main(int argc, char *argv[]) {
  bsp_init(bspcg,argc, argv); //Execute after we start with initalization
  //Handle arguments

  printf("How many processors do you want to use?\n"); 
  fflush(stdout);
  scanf("%d",&P);

  if (P > bsp_nprocs()){
    printf("Sorry, not enough processors available.\n");
    exit(1);
  }
  bspcg();
  exit(0);
}
