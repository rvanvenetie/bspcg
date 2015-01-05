#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cblas.h>
#include "bspedupack.h"


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

#define DEBUG(...) ( if (s == 0) printf(__VA_ARGS__); )

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
  
  if (s == 0) {
    bsp_begin(P);

    int p = bsp_nprocs(); /* p=P */
    int s = bsp_pid();

    char matbuffer[1024], vecdisbuffer[1024], *vecvalbuffer = NULL;
    snprintf( matbuffer, 1024, "../../mtxMatrices/bodyy5.mtx-P%d", p);
    snprintf( vecdisbuffer, 1024, "../../mtxMatrices/bodyy5.mtx-u%d", p);
    distributed_matrix mat = load_symm_distributed_matrix_from_file( matbuffer, p, s);
    vector_distribution dis = load_vector_distribution_from_file( vecdisbuffer, p, s);
    double *vals = vecallocd( dis.nv);
    load_vector_values_from_file( vecvalbuffer, dis, p, s, vals);
    bsp_sync();
    //Read the matrix here
    //bspinput2triple?
    //triple2icrs

    printf("Solving sparse system Ax=b, dim=%i, nz=%i, kmax=%i, eps=%f, with\n", n, nz, kmax, eps);
    printf("Using %d processors\n",p);
  }

  if( PRINTMATS) {
    //Print matrix here? Might be hard, because it is
    //distributed over the processors. 
  }
  /* Matrix is now loaded. Load the vector distribution */
  bspinputvec(p,s,vfilename,&n,&nv,&vindex);
  /* Input vector has local index. Make array v[i] with global index? */
  v= vecallocd(nv);
  for(i=0; i<nv; i++){
    iglob= vindex[i];
    v[i]= iglob+1;
  }

  /* Loading phase done, bsp_sync, add timer?*/
  bsp_sync();

  /* Initalize bsp_mv */
  bspmv_init

  int k = 0; //Iterations
  //Initalize vectors used in algorithm, b is already initialized
  *r = vecallocd(nv); 
  *u = vecallocd(nv);
  *x = vecallocd(nv);
  *w = vecallocd(nv);
  //u and w are initialized in the algorithm.

  //Initial value of x = 0, r = b - Ax = b.
  for (int i = 0; i < nv; i++) {
    x[i] = 0;
    r[i] = b[i];
  }

  //Rho used in algo, initial rho = <r,r> = <b,b>
  double rho, rho_old;
  rho = bppip_dist(p,s,r,r, nv);

  //We store |b|^2 to calculate the relative norm of r
  //|b|^2 = <b,b> = rho
  double normbsq = rho;

  //bsp_sync?
  while( rho > eps*eps*normbsq && k < kmax) {
    double beta, gamma, alpha;
    DEBUG("Iteration %d, rho = %g\n", k + 1, rho);
    
    if( k == 0) {
      cblas_dcopy(nv, r, 1, p, 1); //u \gets r
    } else {
      beta = rho/rho_old;
      //u \gets r + beta u
      cblas_dscal(nv, beta, p, 1); // u \gets \beta p
      cblas_daxpy(nv, 1.0, r, 1, p, 1); // u \gets r + p
    }
    // w \gets Au
    bspmv();
    // \gamma \gets <u, w>
    gamma = bspip(p,s, u,w, nv);
    alpha = rho/gamma;
    //  x \gets x + alpha u
    cblas_daxpy( n, alpha, u, 1, x, 1); 
    //  r \gets r - alpha w
    cblas_daxpy( n, - alpha, w, 1, r, 1); 
    rho_old = rho;
    // \rho \gets <r ,r>
    rho = bspip(p,s, r,r, nv);
    k++;
  }

  DEBUG("Solution found after %i iterations! rho = %f\n", k, rho);
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
  if (s == 0)
    x_glob = vecallocd(x_glob);

  bsp_push_reg(x_glob, n*SZDBL);
  bsp_sync();
  for (int i = 0; i < nv; i++)
  {
    index_glob = bla;
    bsp_put(0, &x[i], x_glob, index_glob * SZDBL, SZDBL);
  }
  bsp_pop_reg(x_glob);
  bsp_sync();
  
  if (s == 0) {
    printf("Solution vector:\n");
    for (int i = 0; i < n; i++)
      printf("%g\n", x_glob[i]);
  }

  vecfreed(u);
  vecfreed(x);
  vecfreed(r);
  vecfreed(w);
  vecfreed(b);
  if (s == 0)
    vecfreed(x_glob);
  bsp_end();
}


int main(int argc, char *argv[]) {
  bsp_init(bspprimes,argc, argv); //Execute bspprimes after we start with initalization
  if (argc == 3) {
    P = atoi(argv[1]);
    N = atoi(argv[2]);
  } else if (BENCHMARK) {
    P = bsp_nprocs();
    N = 10000;
  } else {
    printf("How many processors do you want to use?\n"); 
    fflush(stdout);
    scanf("%d",&P);
    printf("Please enter n:\n"); fflush(stdout);
    scanf("%d",&N);
  }
  if(N<0)
    bsp_abort("Error in input: n is negative");
  if (P > bsp_nprocs()){
      printf("Sorry, not enough processors available.\n");
      exit(1);
  }
  bspprimes();
  exit(0);
}
