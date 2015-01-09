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


/*
 *
 * We have two options: Domain distribution or Matrix distrubtion.
 * - Domain: We split up the domain over the given processors. This leads to a natural
 *   subdivision of the matrix A, where we only need communication for the boundary vertices.
 *   We communicate these values once, which allows us to construct a matrix A_p.
 *   With the property that A = sum_{i=1}^p A_p. 
 *
 *   For CG we have to calculate Au = sum_{i=1}^p (A_p u). So we only have to
 *   calculate (A_p u) locally on the matrices, after which we have to communicate the values
 *   of u lying on the boundary.
 *
 *   Problem with this approach is that we have to find a distrbibution that balances the communication
 *   and computation.
 *
 *   We could solve this by taking a balanced initial distribution, for say p=2 processors. We now
 *   have A = A_1 + A_2, with A_1 and A_2 sparser than the initial matrices. We could then scale this process
 *   by using the Mondriaan package to calculate (A_1 u & A_2 u). We could then use the Mondriaan
 *   package to distribute the load of these vectors accordingly. 
 *
 *
 * - Matrix distribution: We calculate the FEM-matrix. We now
 *
*/
typedef struct {
	double *x, *y;
	int * b;
	int * t[3]; 

	int n_vert;
	int n_tri;
} mesh_t;

generate_element_matrix( double *x, double *y, double *t, int n) { //ofzo
  static double A[9] = {0.5, -0.5, 0, 
                        -0.5, 0.5, 0, 
                        0, 0, 0};
  static double B[9] = {1, -0.5, -0.5,
                        -0.5, 0, 0.5,
                        -0.5, 0.5, 0};
  static double C[9] = {0.5, 0, -0.5,
                        0, 0, 0,
                        -0.5, 0, 0.5};
  static double rhs[3] = {1.0/6.0};
  int v1 = t[n*3 + 0];
  int v2 = t[n*3 + 1];
  int v3 = t[n*3 + 2];
  double f11 = x[v2] - x[v1];
  double f12 = x[v3] - x[v1];
  double f21 = y[v2] - y[v1];
  double f22 = y[v3] - y[v1];
  double D = f11*f22 - f12*f21;
  double E1 = f12*f12 + f22*f22;
  double E2 = f11*f12 + f21*f22;
  double E3 = f11*f11 + f21*f21;

  for( int r = 0; r < 3; r++) {
    for( int s = 0; s <= r; s++) {
      double krs = 1/fabs( D) * (E1*A[r*3+s] - E2*B[r*3+s] + E3*C[r*3+s]); //TODO of dit moet een +E2 zijn
    }

    double gr = fabs(D) * rhs[r]; //local rhs component
  }
}



















#define B_ONES 0
#define B_LOAD 1
#define B_AX 2
#define B_RAND 3
/* Below we have set the default values for the parameters of this program */

/* These variables may be acces by any processor */
int kmax		 = 10000;
int b_mode	 = B_ONES;
double eps	 = 1E-8;
int use_debug= 1;
int use_time = 1;
/* These variables may only be acces by s = 0 */
int P        = 1;


char * matbuffer    = NULL;
char * vecdisbuffer = NULL;
char * vecvalbuffer = NULL;
#define BUFERSTRSIZE 1024

#define DEBUG(...) { if ( use_debug && s == 0) printf(__VA_ARGS__); }

void bspcg() {
  double time_init, time_done, time_total;
  double time_before_mv, time_mv;
  double time_before_ip, time_ip;
  int p,s; //Defaults for bsp
  //Vectors local
  double *b; //Holds rhs  
  double *r; //Holds the residual
  double *x; //Holds the approximate solution
  double *u; //Used in the algorithm to update x
  double *w; //Used in the algorithm to calculate alpha
  distributed_matrix mat;  //Stores the `local' matrix.
  vector_distribution dis; //Stores the vector distribution
	//Arrays needed for bspmv
	int *srcprocv, *srcindv, *destprocu, *destindu;

  bsp_begin(P);
  p= bsp_nprocs(); /* p = number of processors obtained */
  s= bsp_pid();    /* s = processor number */
	/* Communicate the parameters of this program */

	//Push variable names
	bsp_push_reg(&b_mode,			SZINT);
	bsp_push_reg(&kmax,				SZINT);
	bsp_push_reg(&eps,				SZDBL);
	bsp_push_reg(&use_debug,  SZINT);
	bsp_push_reg(&use_time,   SZINT);
	bsp_sync();

	//Read values from s = 0
	bsp_get(0, &b_mode, 0, &b_mode, SZINT);
	bsp_get(0, &kmax  , 0, &kmax  , SZINT);
	bsp_get(0, &eps   , 0, &eps   , SZDBL);

	bsp_get(0, &use_debug, 0, &use_debug, SZINT);
	bsp_get(0, &use_time , 0, &use_time , SZINT);
	bsp_sync();

	bsp_pop_reg(&b_mode);
	bsp_pop_reg(&kmax);
	bsp_pop_reg(&eps);
	bsp_pop_reg(&use_debug);
	bsp_pop_reg(&use_time);

  /* Read matrix/vector distribution&values from file */
  mat = load_symm_distributed_matrix_from_file( matbuffer, p, s);
	dis = load_vector_distribution_from_file( vecdisbuffer, vecvalbuffer, p, s, &b);

	if (use_debug) {
      printf( "s/p: %d/%d\n"
					    "\tkmax: %d\n"
							"\teps : %g\n"
							"\tmat.n: %d\n"
							"\tmat.nz: %d\n"
							"\tmat.nrows: %d\n"
							"\tmat.ncols: %d\n"
							"\tdis.n: %d\n"
							"\tdis.nv: %d\n"
							, s,p,kmax,eps, mat.n, mat.nz, mat.nrows, mat.ncols, dis.n, dis.nv);
	}
	bsp_sync();
	/* Matrix, vector dist and b is loaded */
	/* Loading phase done, bsp_sync, add timer?*/

	/* Initalize bsp_mv */
  srcprocv  = vecalloci(mat.ncols);
  srcindv   = vecalloci(mat.ncols);
  destprocu = vecalloci(mat.nrows);
  destindu  = vecalloci(mat.nrows);
  bspmv_init(p, s, mat.n, mat.nrows, mat.ncols,
             dis.nv, dis.nv, mat.rowindex, mat.colindex,
	     dis.vindex, dis.vindex, 
	     srcprocv, srcindv, destprocu, destindu);
	if (b_mode != B_LOAD) { //Not loaded from file, initialize as a vector of ones
		//Initialize b to be a vector of ones
    b = vecallocd(dis.nv);
    if( b_mode == B_AX) {
      double *x_temp = vecallocd(dis.nv);
      for( int i = 0; i < dis.nv; i++) 
        x_temp[i] = 1;
      bspmv(p,s, mat.n, mat.nz, mat.nrows, mat.ncols, 
            mat.val, mat.inc,
            srcprocv, srcindv, destprocu, destindu,
            dis.nv, dis.nv, x_temp, b);
      vecfreed( x_temp);
    } else if( b_mode == B_ONES) {
      for( int i = 0; i < dis.nv; i++) 
        b[i] = 1;
    } else if( b_mode == B_RAND) {
      for( int i = 0; i < dis.nv; i++) 
        b[i] = ((double)rand())/RAND_MAX;
    }
  }

  int k = 0; //Iterations
  //Allocate vectors used in algorithm, b is already allocated
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

  time_init = bsp_time();
	/* All initalization is done, start the main loop */
  while( rho > eps*eps*normbsq && k < kmax) {
    double beta, gamma, alpha;
    DEBUG("Iteration %d, rho = %g\n", k + 1, rho);
    if( k > 0) {
      beta = rho/rho_old;
      /* u \gets r + beta u */
			
			// u \gets \beta u
      vec_scale(dis.nv, beta, u);
    }

		// u \gets r + u
    vec_axpy(dis.nv, 1.0, r, u); 

    // w \gets Au
    time_before_mv = bsp_time();
    bspmv(p,s, mat.n, mat.nz, mat.nrows, mat.ncols, 
					mat.val, mat.inc,
					srcprocv, srcindv, destprocu, destindu,
					dis.nv, dis.nv, u, w);
    time_mv += bsp_time() - time_before_mv;

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
  }
  //bsp_sync?
  time_done = bsp_time();

  if (k < kmax) {
    DEBUG("Solution found after %i iterations! rho = %f\n", k, rho);
  } else {
    DEBUG("CG stopped, maximum iterations (%i) reached. rho = %g\n", k, rho);
  }

	/*
	 * We now give processor zero the solution vector.
	 * Note that we probably do not want to do this in the real
	 * application, as this vector might be too big to store
	 * on one processor.

	 * We have this step mainly for timing
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

  time_total = bsp_time();

	if (use_debug) {
		if (s == 0)
			printf("mat: %s\np: %d\nn: %d\nb_mode: %d\nk: %d\n", basename(matbuffer), p, mat.n, b_mode, k);
		printf( "s: %d\n"
						"\ttime_init: %g\n"
						"\ttime_mv: %g\n"
						"\ttime_ip: %g\n"
						"\ttime_done: %g\n"
						"\ttime_total: %g\n"
							, s, time_init, time_mv, time_ip, time_done, time_total);
	}
	if (s == 0 && use_time) {
		double time_iter = time_done - time_init;
		double time_glob = time_total - time_done;
		double time_it_local = time_iter - (time_mv + time_ip);
		double density = mat.nzA / ((double) mat.n * mat.n);
		if (use_debug)
			printf("mat_name, p,      mat_n,  mat_nzA, density, k, time_init, time_iter, time_glob,time_total,  time_it_local, time_mv, time_ip\n");
		//mat_name, p,      mat_n,  mat_nzA, density, k,  time_init, time_iter, time_glob, time_total time_it_local, time_mv, time_ip
		printf("%s" "\t%d" "\t%d" "\t%d"    "\t%6f"	"\t%d"	"\t%6f"		 "\t%6f" 		"\t%6f"	 "\t%6f"		"\t%6f"					"\t%6f"   "\t%6f\n",
				basename(matbuffer), p, mat.n, mat.nzA,density,k,time_init,time_iter,time_glob,time_total, time_it_local,time_mv,time_ip);
  }

  bsp_end();
}


int main(int argc, char *argv[]) {
  bsp_init(bspcg,argc, argv); //Execute after we start with initalization
	P = bsp_nprocs();
	if(!(argc == 2 || argc == 3 || argc == 6 || argc == 8))
	{
		fprintf(stderr, "Invalid arguments given\n"
				            "\tUsage: %s mtx_file [P [b_mode kmax eps [use_time use_debug]]]\n"
										"\tDefault parameters:\n"
										"\t\tP=%d\n"
										"\t\tb_mode=%d\n"
										"\t\tkmax=%d\n"
										"\t\teps=%g\n"
										"\t\tuse_time=%d\n"
										"\t\tise_debug=%d\n",
										argv[0],P,b_mode, kmax, eps, use_time, use_debug);
		exit(1);
	}
	if (argc > 2) {
		P =  atoi(argv[2]);
		if (P > bsp_nprocs()){
			fprintf( stderr, "Sorry, not enough processors available.\n");
			exit(1);
		}
	}
	if (argc > 3) {
		b_mode = atoi(argv[3]);
		kmax =   atoi(argv[4]);
		eps  =   atof(argv[5]);
	}

	if (argc > 6) {
		use_time =  atoi(argv[6]);
		use_debug = atoi(argv[7]);
	}


	matbuffer   =  malloc(BUFERSTRSIZE);
	vecdisbuffer=  malloc(BUFERSTRSIZE);
	snprintf( matbuffer,    1024, "%s-P%d", argv[1], P);
	snprintf( vecdisbuffer, 1024, "%s-u%d", argv[1], P);
	if (b_mode == B_LOAD) {
		vecvalbuffer=  malloc(BUFERSTRSIZE);
		snprintf( vecvalbuffer, 1024, "%s-b"  , argv[1]);
	}
  bspcg();
	free(matbuffer);
	free(vecdisbuffer);
	free(vecvalbuffer);
  exit(0);
}
