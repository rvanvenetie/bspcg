#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <libgen.h> 
#include <cblas.h>
#include "mesh.h"
#include "bspedupack.h"
#include "vec.h"
#include "fem.h"

/* These variables may be acces by any processor */
int kmax		 = 10000;
double eps	 = 1E-8;
int use_debug= 1;
int use_time = 1;
/* These variables may only be acces by s = 0 */
int P        = 1;
char * meshbuffer   = NULL;
#define BUFFERSTRSIZE 1024

void print_mesh_dist(mesh_dist mesh) {
	printf("(x,y) : boundary\n");
	for (int i = 0; i < mesh.n_vert; i++) {
		printf("(%3f, %3f): %d\n", mesh.x[i], mesh.y[i], mesh.b[i]);
	}
	printf("(v1,v2,v3) : proc\n");
	for (int i = 0; i < mesh.n_tri; i++) {
		printf("(%d, %d, %d) : %d\n", mesh.t[i][0], mesh.t[i][1], mesh.t[i][2], mesh.p[i]);
	}
	printf("n_vert: %d\n", mesh.n_vert);
	printf("n_tri:  %d\n", mesh.n_tri);
	printf("p: %d\n", mesh.P);

	printf("N = np.array([");
	for (int i = 0; i < mesh.n_vert; i++)
		printf("[%.6f, %.6f]\n", mesh.x[i], mesh.y[i]);
	printf("])\n");

	printf("T = np.array([");
	for (int i = 0; i < mesh.n_tri; i++)
		printf("[%d, %d, %d]\n", mesh.t[i][0], mesh.t[i][1], mesh.t[i][2]);
	printf("])\n");

	printf("G = np.array([");
	for (int i = 0; i < mesh.n_vert; i++)
		printf("%d,",mesh.b[i]);
	printf("])\n");
}

void print_fem_vect(int s, char * name, bsp_fem_data * fem, mesh_dist * mesh_total, double * x) {
	double * x_glob = (double * ) 1337; //Hacky solution, otherwise BSP will struggle
	if (s == 0)
		x_glob = vecallocd(fem->n_vert_total);

	bsp_push_reg(x_glob, fem->n_vert_total*SZDBL);
	bsp_sync();

	for (int i = 0; i < fem->n_shared; i++) //Only push shared vertex if we own it
		if (s == proc_owner(fem->p_shared[i]))  
			bsp_put(0, &x[i], x_glob, fem->i_glob[i] * SZDBL, SZDBL);
		
	for (int i = fem->n_shared; i < fem->dof; i++) 
		bsp_put(0, &x[i], x_glob, fem->i_glob[i] * SZDBL, SZDBL);
	
	bsp_sync();
	bsp_pop_reg(x_glob);
	if (s == 0) { 
		int dof_cntr = 0;
		for (int i = 0; i < mesh_total->n_vert; i++)
			if (!mesh_total->b[i])
				x_glob[dof_cntr++] = x_glob[i];
		fprintf(stdout,"FEM Vector %s:\n[", name);
		for (int i = 0; i < dof_cntr; i++)
			fprintf(stderr,"%.5f\n", x_glob[i]);
		fprintf(stdout, "]\n");
		//printf("]\n");
		vecfreed(x_glob);
	}
}

void print_fem_mat(int s, char * name, bsp_fem_data * fem, mesh_dist * mesh_total) {
	int tagsz= SZINT;
	bsp_set_tagsize(&tagsz);
	double * mat_full = (double *) 1337; //Hack solution
	if (s == 0) 
		mat_full = calloc(sizeof(double),fem->n_vert_total * fem->n_vert_total);

	bsp_push_reg(mat_full, fem->n_vert_total*fem->n_vert_total * SZDBL);
	bsp_sync();

	//Every processor sends its matrix values to the main processor
	for (int k = 0; k < fem->mat.nz; k++) {
		int i = fem->i_glob[fem->mat.I[k]];
		int j = fem->i_glob[fem->mat.J[k]];
		int mat_idx = i*fem->n_vert_total + j;
		bsp_send(0,  &mat_idx, &fem->mat.val[k], SZDBL);
		if (i != j) { //Send it's transpose as well
			mat_idx = j*fem->n_vert_total + i;
			bsp_send(0,  &mat_idx, &fem->mat.val[k], SZDBL);
		}
	}
	bsp_sync();
	bsp_pop_reg(mat_full);
	if (s == 0) {
		int status, tag;
		double tmp_val;
		bsp_get_tag(&status, &tag);
		while (status != -1) { //Process all the messages
			bsp_move(&tmp_val, SZDBL);

			//Tag holds the global index, tmp_val the value
			mat_full[tag] += tmp_val;
			bsp_get_tag(&status, &tag);
		}
		printf("FEM %s\n", name);
		printf("[");
		for (int i = 0; i < mesh_total->n_vert; i++)
		  if (!mesh_total->b[i])  {
				for (int j = 0; j < mesh_total->n_vert; j++)
					if (!mesh_total->b[j]) 
						printf("%.3f ", mat_full[i * mesh_total->n_vert + j]);
				printf(";\n");
			}
		printf("]");
	}

}
void print_fem_data(int s, bsp_fem_data result, mesh_dist * mesh_total) {
	printf("%d\n"
			    "\tn_own:%d\n"
					"\tn_shared:%d\n"
					"\tn_bound:%d\n"
					"\tn_vert:%d\n"
					"\tdof:%d\n"
					"\tn_tri:%d\n"
					"\tn_vert_total:%d\n",
					s, result.n_own, result.n_shared, result.n_bound, result.n_vert, result.dof, result.n_tri, result.n_vert_total);
	if (SUPER_DEBUG) {
    //Print vertices
		for (int i = 0; i < result.n_vert; i++) 
			printf( "%d : (%3f, %3f) \t loc:%d glob:%d\n",
					s, result.x[i], result.y[i], i, result.i_glob[i]);
	  //Print triangles
		for (int i = 0; i < result.n_tri; i++)
			printf("%d : (%d, %d, %d)\n", s, result.t[i][0], result.t[i][1], result.t[i][2]);

		printf("\n\n");
		print_fem_mat(s, "FEM Matrix", &result, mesh_total);
		print_fem_vect(s, "FEM Rhs", &result, mesh_total, result.rhs);
		printf("\n\n");
	}
}


#define DEBUG(...) { if ( use_debug && s == 0) printf(__VA_ARGS__); }

void bspfem() {
  double time_init, time_done, time_total;
  double time_before_mv, time_mv;
  double time_before_ip, time_ip;
  int p,s; //Defaults for bsp
  //Vectors local
  double *r; //Holds the residual
  double *x; //Holds the approximate solution
  double *u; //Used in the algorithm to update x
  double *w; //Used in the algorithm to calculate alpha
	bsp_fem_data fem; //Stores the `local' FEM data
	mesh_dist mesh_total; //Processor 0 uses this

  bsp_begin(P);
  p= bsp_nprocs(); /* p = number of processors obtained */
  s= bsp_pid();    /* s = processor number */
	/* Communicate the parameters of this program */

	//Push variable names
	bsp_push_reg(&kmax,				SZINT);
	bsp_push_reg(&eps,				SZDBL);
	bsp_push_reg(&use_debug,  SZINT);
	bsp_push_reg(&use_time,   SZINT);
	bsp_sync();

	//Read values from s = 0
	bsp_get(0, &kmax  , 0, &kmax  , SZINT);
	bsp_get(0, &eps   , 0, &eps   , SZDBL);
	bsp_get(0, &use_debug, 0, &use_debug, SZINT);
	bsp_get(0, &use_time , 0, &use_time , SZINT);
	bsp_sync();

	bsp_pop_reg(&kmax);
	bsp_pop_reg(&eps);
	bsp_pop_reg(&use_debug);
	bsp_pop_reg(&use_time);

	//TODO: Load mesh_from file
	if (s == 0) {
		FILE * mesh_buf = fopen(meshbuffer, "r"); //TODO: FILE --> FILENAME
		mesh_total = readfrommeshfile(mesh_buf);
		fclose(mesh_buf);
		if (SUPER_DEBUG) 
			print_mesh_dist(mesh_total);
	}
	fem = bsp_fem_init(s,p, &mesh_total);

	if (use_debug)
		print_fem_data(s, fem, &mesh_total);
	/* FEM matrix + RHS are loaded */

  int k = 0; //Iterations
  //Allocate vectors used in algorithm
  r = vecallocd(fem.dof); 
  u = vecallocd(fem.dof);
  x = vecallocd(fem.dof);
  w = vecallocd(fem.dof);

  //u and w are initialized in the algorithm.
  //Initial value of x = 0, r = b - Ax = b.
  for (int i = 0; i < fem.dof; i++) {
    x[i] = 0;
    u[i] = 0;
    r[i] = fem.rhs[i];
  }

		
  //Rho used in algo, initial rho = <r,r> = <b,b>
  double rho, rho_old;
	rho = bsp_fem_ip(p,s,r,r, &fem);
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
      // u \gets r + beta u 
			
			// u \gets \beta u
      vec_scale(fem.dof, beta, u);
    }

		// u \gets r + u
    vec_axpy(fem.dof, 1.0, r, u); 

    // w \gets Au
    time_before_mv = bsp_time();
		bsp_fem_mv(s,&mesh_total, &fem,w, u);
    time_mv += bsp_time() - time_before_mv;

    // \gamma \gets <u, w>
    time_before_ip = bsp_time();
		gamma = bsp_fem_ip(p,s, u,w, &fem);
    time_ip += bsp_time() - time_before_ip;
    alpha = rho/gamma;

    //  x \gets x + alpha u

    vec_axpy(fem.dof, alpha, u, x); 

    //  r \gets r - alpha w
    vec_axpy(fem.dof, - alpha, w, r); 

    rho_old = rho;
    // \rho \gets <r ,r>

    time_before_ip = bsp_time();
		rho = bsp_fem_ip(p,s, r,r, &fem);
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
	if (use_debug)
		print_fem_vect(s,"x_solution", &fem, &mesh_total, x);
	/*
	 * We now give processor zero the solution vector.
	 * Note that we probably do not want to do this in the real
	 * application, as this vector might be too big to store
	 * on one processor.

	 * We have this step mainly for timing
	 */
	/*
	double * x_glob;
	x_glob = vecallocd(dis.n);
	bsp_push_reg(x_glob, dis.n*SZDBL);
	bsp_sync();

	for (int i = 0; i < fem.dof; i++)
	{
		int index_glob = dis.vindex[i];
		bsp_put(0, &x[i], x_glob, index_glob * SZDBL, SZDBL);
	}

	bsp_sync();
	bsp_pop_reg(x_glob);
	vecfreed(x_glob);
	*/
	if (s == 0)
		mesh_free(&mesh_total);
  vecfreed(u);
  vecfreed(x);
  vecfreed(r);
  vecfreed(w);
  time_total = bsp_time();
	if (use_debug && s==0) {
/*		if (s == 0)
			printf("mat: %s\np: %d\nn: %d\nb_mode: %d\nk: %d\n", basename(matbuffer), p, mat.n, b_mode, k);
			*/
		printf( "s: %d\n"
						"\ttime_init: %g\n"
						"\ttime_mv: %g\n"
						"\ttime_ip: %g\n"
						"\ttime_done: %g\n"
						"\ttime_total: %g\n"
							, s, time_init, time_mv, time_ip, time_done, time_total);
	}
	if (s == 0 && use_time) {
		/*
		double time_iter = time_done - time_init;
		double time_glob = time_total - time_done;
		double time_it_local = time_iter - (time_mv + time_ip);
		double density = mat.nzA / ((double) mat.n * mat.n);
		if (use_debug)
			printf("mat_name, p,      mat_n,  mat_nzA, density, k, time_init, time_iter, time_glob,time_total,  time_it_local, time_mv, time_ip\n");
		//mat_name, p,      mat_n,  mat_nzA, density, k,  time_init, time_iter, time_glob, time_total time_it_local, time_mv, time_ip
		printf("%s" "\t%d" "\t%d" "\t%d"    "\t%6f"	"\t%d"	"\t%6f"		 "\t%6f" 		"\t%6f"	 "\t%6f"		"\t%6f"					"\t%6f"   "\t%6f\n",
				basename(matbuffer), p, mat.n, mat.nzA,density,k,time_init,time_iter,time_glob,time_total, time_it_local,time_mv,time_ip);
				*/
  }
	fem_data_free(&fem);
  bsp_end();
}


int main(int argc, char *argv[]) {
  bsp_init(bspfem,argc, argv); //Execute after we start with initalization
	P = bsp_nprocs();
	if(!(argc == 2 || argc == 3 || argc == 5 || argc == 6))
	{
		fprintf(stderr, "Invalid arguments given\n"
				            "\tUsage: %s mesh file [P [kmax eps [use_time use_debug]]]\n"
										"\tDefault parameters:\n"
										"\t\tP=%d\n"
										"\t\tkmax=%d\n"
										"\t\teps=%g\n"
										"\t\tuse_time=%d\n"
										"\t\tise_debug=%d\n",
										argv[0],P,kmax, eps, use_time, use_debug);
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
		kmax =   atoi(argv[3]);
		eps  =   atof(argv[4]);
	}

	if (argc > 6) {
		use_time =  atoi(argv[5]);
		use_debug = atoi(argv[6]);
	}

	meshbuffer  = malloc(BUFFERSTRSIZE);
	snprintf( meshbuffer, BUFFERSTRSIZE, "%s-P%d", argv[1], P);
	bspfem();
	free(meshbuffer);
  exit(0);
}
