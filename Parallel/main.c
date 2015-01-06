#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <libgen.h> 
#include "io.h"
#include "vec.h"
#include "bspmv.h"
#include "bspip.h"
#include "bspedupack.h"

int P;
int use_time, use_debug, kmax, load_b;
double eps;
char *mtx_file;

#define DEBUG(...) { if ( use_debug && s == 0) printf(__VA_ARGS__); }

void bspcg() {
  double time_init, time_done, time_total;
  double time_before_mv, time_mv;
  double time_before_ip, time_ip;
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
  snprintf( matbuffer, 1024, "%s-P%d", mtx_file, p);
  snprintf( vecdisbuffer, 1024, "%s-u%d", mtx_file, p);

  distributed_matrix mat = load_symm_distributed_matrix_from_file( matbuffer, p, s);
  vector_distribution dis;
  if( load_b) {
    dis = load_vector_distribution_from_file( vecdisbuffer, vecvalbuffer, p, s, &b);
  } else {
    dis = load_vector_distribution_from_file( vecdisbuffer, NULL, p, s, &b);
    b = vecallocd(dis.nv);
    for( int i = 0; i < dis.nv; i++) {
      b[i] = 1;
    }
  }
  bsp_sync();
  if (s == 0) {
    DEBUG("Solving sparse system Ax=b, dim=%i, nz=%i, kmax=%i, eps=%f, with\n", mat.n, mat.nz, kmax, eps);
    DEBUG("Using %d processors\n",p);
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
    u[i] = 0;
    r[i] = b[i];
  }

  //Rho used in algo, initial rho = <r,r> = <b,b>
  double rho, rho_old;
  rho = bspip_dist(p,s,r,r, dis.nv);

  //We store |b|^2 to calculate the relative norm of r
  //|b|^2 = <b,b> = rho
  double normbsq = rho;
  //bsp_sync?
  time_init = bsp_time();
  while( rho > eps*eps*normbsq && k < kmax) {
    double beta, gamma, alpha;
    DEBUG("Iteration %d, rho = %g\n", k + 1, rho);
    
    if( k > 0) {
      beta = rho/rho_old;
      //u \gets r + beta u
      vec_scale(dis.nv, beta, u); // u \gets \beta u
    }
    vec_axpy(dis.nv, 1.0, r, u); // u \gets r + u
    // w \gets Au
    time_before_mv = bsp_time();
    bspmv(p,s, mat.n, mat.nz, mat.nrows, mat.ncols, 
	  mat.val, mat.inc,
	  srcprocv, srcindv, destprocu, destindu,
	  dis.nv, dis.nv, u, w);
    time_mv += bsp_time() - time_before_mv;
    //vec_print("w", dis.nv, w);
    // \gamma \gets <u, w>
    time_before_ip = bsp_time();
    gamma = bspip_dist(p,s, u,w, dis.nv);
    time_ip += bsp_time() - time_before_ip;
    alpha = rho/gamma;
    //  x \gets x + alpha u
    vec_axpy(dis.nv, alpha, u, x); 
    //  r \gets r - alpha w
    vec_axpy(dis.nv, - alpha, w, r); 
    rho_old = rho;
    // \rho \gets <r ,r>

    time_before_ip = bsp_time();
    rho = bspip_dist(p,s, r,r, dis.nv);
    time_ip += bsp_time() - time_before_ip;
    k++;
    //vec_print("x", dis.nv, x);
    //vec_print("r", dis.nv, r);
  }
  time_done = bsp_time();

  if (k < kmax) {
    DEBUG("Solution found after %i iterations! rho = %f\n", k, rho);
  } else {
    DEBUG("CG stopped, maximum iterations (%i) reached. rho = %g\n", k, rho);
  }
  //bsp_sync?
  //Calculate w = Ax
  //Print timing and resulting vector w/b 

  if( 0) {
    /*
     * We now give processor zero the solution vector.
     * Note that we probably do not want to do this in the real
     * application, as this vector might be too big to store
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
  }

  vecfreed(u);
  vecfreed(x);
  vecfreed(r);
  vecfreed(w);
  vecfreed(b);
  vecfreei(srcprocv);
  vecfreei(srcindv);
  vecfreei(destprocu);
  vecfreei(destindu);
  time_total = bsp_time();
  if( use_time) {
    //p s n k init done total mv ip
    if( use_debug) {
      if( s== 0) {
        printf("mat: %s\np: %d\nn: %d\nk: %d\n", basename(mtx_file), p, mat.n, k);
      }
      printf(
          "s: %d\n"
          "\ttime_init: %g\n"
          "\ttime_mv: %g\n"
          "\ttime_ip: %g\n"
          "\ttime_done: %g\n"
          "\ttime_total: %g\n"
          , s, time_init, time_mv, time_ip, time_done, time_total);
    } else {
      printf("%s\t%d\t%d\t%d\t%d\t%6f\t%6f\t%6f\t%6f\t%6f\n", basename(mtx_file), p, s, mat.n, k, time_init, time_mv, time_ip, time_done, time_total);
    }
  }
  bsp_end();

}


int main(int argc, char *argv[]) {
  bsp_init(bspcg,argc, argv); //Execute after we start with initalization

  //Handle arguments
  if( argc == 6) {
    mtx_file = argv[1];
    load_b = (int) strtol( argv[2], (char **)NULL, 10);
    P = (int) strtol( argv[3], (char **)NULL, 10);
    kmax = (int) strtol( argv[4], (char **)NULL, 10);
    eps = atof( argv[5]);
    use_time = 1;
    use_debug = 1;
  } else if( argc == 8) {
    mtx_file = argv[1];
    load_b = (int) strtol( argv[2], (char **)NULL, 10);
    P = (int) strtol( argv[3], (char **)NULL, 10);
    kmax = (int) strtol( argv[4], (char **)NULL, 10);
    eps = atof( argv[5]);
    use_time = (int) strtol( argv[6], (char **)NULL, 10);
    use_debug = (int) strtol( argv[7], (char **)NULL, 10);
  } else {
    printf("Sorry, no workie. Usage: %s mtx_file load_b P kmax eps [use_time use_debug]\n", argv[0]);
    exit(1);
  }

  if (P > bsp_nprocs()){
    printf("Sorry, not enough processors available.\n");
    exit(1);
  }
  bspcg();
  exit(0);
}
