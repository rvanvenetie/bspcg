#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <libgen.h> 
#include <cblas.h>
#include "mesh.h"
#include "io.h"
#include "vec.h"
#include "bspmv.h"
#include "bspip.h"
#include "bspedupack.h"


#define B_ONES 0
#define B_LOAD 1
#define B_AX 2
#define B_RAND 3
/* These variables may be acces by any processor */
int kmax		 = 10000;
int b_mode	 = B_ONES;
double eps	 = 1E-8;
int use_debug= 1;
int use_time = 1;
/* These variables may only be acces by s = 0 */
int P        = 1;
char * meshbuffer   = NULL;
#define BUFFERSTRSIZE 1024


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

	/*
	 * If a vertex is shared among multiple processors we have to decide who
	 * is going to be the 'owner' of the vertex. We simply give the processor
	 * with the smallest number the owner. 
	 * 
	 * TODO: Even out the amount of vertices owned by a processor, this implies
	 * a better distribution of the vector among the processors. 
	 */

/* Functions to create a processor set */
typedef long long unsigned int proc_set;
#define SZPROCSET sizeof(proc_set)
#define SZTRI     sizeof(triangle)
 
/* Data structure for the data of each processor */
typedef struct {
	/*
	 * Vertices used in this processor. We have the following order:
	 * 1) Shared vertices 
	 * 2) Owned vertices
	 * 3) Boundary vertices (all the vertices lying on the boundary
	 */
  double *x, *y;			 //Vertices 
	int n_shared;        //Amount of shared vertices (not bound)
	int n_own;           //Amount of owned vertices (not bound)
	int n_bound;         //Amount of vertices on the boundary
	int n_vert;          //Total amount of vertices (n_shared + n_own + n_bound)

	int dof;					 	 //Degrees of freedom, in our case n_shared + n_own
	proc_set * p_shared; //Processor sets for the shared vertices

	int n_vert_total;    //Total amount of vertices in the system
  int * i_glob;				 //Global index that corresponds to the local vertex
	int * glob2local;    //Lazy way to convert global to local indices 

  //int ** i_shared;		 //Array of (p,i_loc), giving the local index on processor p for this vertex

  /* Triangles on this processor. Should be LOCAL indices */
  int n_tri; 
	triangle * t;

  /* FEM matrix for this processor */
  matrix_s mat;
} bsp_fem_data;


void gen_element_matrix(double * res, double *x, double *y, triangle t) {
  //http://en.wikipedia.org/wiki/Stiffness_matrix (Die matrix D onderin)
  double D_mat[] = {x[t[2]] - x[t[1]], x[t[0]] - x[t[2]], x[t[1]] - x[t[0]],
	  							  y[t[2]] - y[t[1]], y[t[0]] - y[t[2]], y[t[1]] - y[t[0]]};

  //Area triangle times 2, also det(coor_mat)
  double area = fabs((x[t[1]] - x[t[0]]) * (y[t[2]] - y[t[0]]) - (x[t[2]] - x[t[0]]) * (y[t[1]] - y[t[0]]));
	double alpha = 1 / (2 * area);

	//Element matrix is now given by: alpha * D.T * D.
	//Only calculate lower tri part
	for (int i = 0; i < 3; i++)
		for (int j = 0; j <= i; j++)
		{ 
			//Row i of D.T times colum j of D
			res[i*3 + j] = alpha * ( D_mat[i] * D_mat[j] + D_mat[i + 3] * D_mat[j+3]);
		}
	
	/*
  for (int i = 0; i < 3; i++)
   for (int j = 0; j <= i; j++)
     fprintf(stderr,"Res[%d,%d] = %g\n", i,j, res[i*3 + j]);

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
  int v1 = t[0];
  int v2 = t[1];
  int v3 = t[2];
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
      fprintf(stderr,"Bla[%d,%d] = %g\n", r,s, krs);
    }

    double gr = fabs(D) * rhs[r]; //local rhs component
		fprintf(stderr,"Bla_rhs[%d] = %g\n", r, gr);
  }
	*/
}


matrix_s gen_fem_mat(double *x, double *y, int dof,
		                 triangle * t, int n_tri) {
	matrix_s result = mat_init(dof, dof);
	for (int k = 0; k < n_tri; k++) {
		//Calculate element matrix for triangle t[k]
		double el_mat[3*3];
		gen_element_matrix(el_mat, x,y,t[k]);
		for (int li = 0; li < 3; li++)
			for (int lj = 0; lj <= li; lj++) //Local indices in element matrix
			{
				int i = t[k][li];
				int j = t[k][lj];
				//i and j hold the vertex index, only store if they account for a DOF
				if (i < dof && j < dof)
					mat_append(&result, i, j, el_mat[i*3 + j]); 
			}
	}
	/* TODO: Matrix is a list of (i,j,val) with duplicates in the list. Ideally 
		 we want to use CRS storage, or at least filter the duplicates */
	return result;
}



/* Count the amount of processors in a set */
int proc_count(proc_set set) {
	return	__builtin_popcountll(set);
}

/* Add a processor to a set */
void proc_add(proc_set * set, int proc) {
  (*set) |= 1ull << proc;
}

/* Remove a processor from a set */
void proc_remove(proc_set * set, int proc) {
	(*set) &= (1ull << proc);
}

/* Returns whether a given processor is in a set */
int proc_test(proc_set set, int proc) {
	return ((set & (1ull << proc)) > 0);
}

/* Returns the next processor (>= proc) in the set, returns -1 if no processor is found */
int proc_next(proc_set set, int proc)  {
	for (int i = proc; i < 64; i++)
		if (proc_test(set, i))
			return i;

	return -1;
}

/* Returns the owner of a given processor set (processor with smallest proc number) */
int proc_owner(proc_set set)  {
	return proc_next(set, 0);
}

/* Calculates the processor sets for each of the processors and
	 the amount of owned/shared vertices per processor */
void mesh_proc_sets(int p, mesh_dist * mesh,
									  int * n_tri_proc,
		                proc_set * vert_proc, int * n_own_proc, int * n_shared_proc, int * n_bound_proc) {
		 						
	//Loop over triangles to determine processor sets and amount of tri per proc
	for (int i = 0; i < mesh->n_tri; i++)  {
		n_tri_proc[mesh->p[i]]++;
		//Loop over vertices of triangle t[i]
		for (int v = 0; v < 3; v++) 
			//Processor that owns triangle t[i] uses vertex v
			proc_add(&vert_proc[mesh->t[i][v]], mesh->p[i]);
	}
  /*
	 * For each processor we now calculate:
	 *  1) Amount of shared vertices (non bound)
	 *  2) Amount of owned verices (non bound)
	 *  3) Amount of boundary vertices  (shared/own)
	 */

	for (int v = 0; v < mesh->n_vert; v++)
	{
		int v_shared = (proc_count(vert_proc[v]) > 1); 
		int v_bound  = (mesh->b[v]);

		//Loop over all the processors using this vertex
		int proc = -1;
		while ((proc = proc_next(vert_proc[v], proc+1)) != -1) 
		{
			if (v_bound)
				n_bound_proc[proc]++;
			else if (v_shared)
				n_shared_proc[proc]++;
			else
				n_own_proc[proc]++;
		}
	}
}

bsp_fem_data bsp_fem_init(int s, int p, mesh_dist * mesh) {
	/*
	 * We assume that s=0 has the entire initial distribution.
	 * This might not be the case in real applications.. TODO:
	 * storing the mesh in a split way, must include ghost cells. Bla bla
	 */

	/*
	 * We need to distribute the triangles and vertices to the designated processors.
	 * We do not want to send unneccesary data. For every vertex we create a \emph{set}
	 * of processors that is using this vertex.
	 *
	 * In c++ this can be implemented easily using bitvector/set, we use a 64-bit int.
	 * Hacky/lazy solution. This implies a restriction of p <= 64!
	 */
	if (p > 64)
	{
		fprintf(stderr, "Cannot handle more than 64 processors.");
		exit(0);
	}
	bsp_fem_data result;
	bsp_push_reg(&result.n_own,   SZINT);
	bsp_push_reg(&result.n_shared, SZINT);
	bsp_push_reg(&result.n_bound, SZINT);
	bsp_push_reg(&result.n_tri,    SZINT);
	bsp_push_reg(&result.n_vert_total, SZINT);
	bsp_sync();

	proc_set * vert_proc = NULL; 
	int *     n_own_proc = NULL;
	int *	 n_shared_proc = NULL;
	int *   n_bound_proc = NULL;
	int * 	  n_tri_proc = NULL;
	if (s == 0) {
		vert_proc					= calloc(sizeof(proc_set), mesh->n_vert); 
		n_own_proc				= calloc(sizeof(int), p);
		n_shared_proc     = calloc(sizeof(int), p);
		n_bound_proc			= calloc(sizeof(int), p);
		n_tri_proc        = calloc(sizeof(int), p);

		//Calculate processor sets, amount of shared/owned vertices and amount of triangles
		mesh_proc_sets(p, mesh, n_tri_proc, vert_proc, n_own_proc, n_shared_proc, n_bound_proc);
		for (int proc = 0; proc < p; proc++)
		{
			bsp_put(proc, &n_own_proc[proc],    &result.n_own   , 0, SZINT);
			bsp_put(proc, &n_shared_proc[proc], &result.n_shared, 0, SZINT);
			bsp_put(proc, &n_bound_proc[proc], &result.n_bound  , 0, SZINT);
			bsp_put(proc, &n_tri_proc[proc],    &result.n_tri,    0, SZINT);
			bsp_put(proc, &mesh->n_vert   , &result.n_vert_total, 0, SZINT);
		}
	}
	bsp_sync();
	//Every processor now knows the amount of vertices (shared/owned) and amount of triangles it holds.
	bsp_pop_reg(&result.n_own);
	bsp_pop_reg(&result.n_shared);
	bsp_pop_reg(&result.n_bound);
	bsp_pop_reg(&result.n_tri);
	bsp_pop_reg(&result.n_vert_total);

	result.n_vert = result.n_shared + result.n_own + result.n_bound;
	result.dof    = result.n_shared + result.n_own;

	result.x			  = malloc(sizeof(double)     * result.n_vert);
	result.y			  = malloc(sizeof(double)     * result.n_vert);

	result.i_glob   = malloc(sizeof(int)        * result.n_vert);
	result.p_shared = malloc(sizeof(proc_set)   * result.n_shared);
	result.t			  = malloc(sizeof(triangle)* result.n_tri);

  bsp_push_reg(result.x,			  SZDBL * result.n_vert);
  bsp_push_reg(result.y,			  SZDBL * result.n_vert);
	bsp_push_reg(result.i_glob,   SZINT * result.n_vert);
	bsp_push_reg(result.p_shared, sizeof(proc_set) * result.n_shared);
	bsp_push_reg(result.t,        sizeof(triangle) * result.n_tri);
	fprintf(stderr,"%d\n"
			    "\tn_own:%d\n"
					"\tn_shared:%d\n"
					"\tn_bound:%d\n"
					"\tn_tri:%d\n"
					"\tn_vert_total:%d\n",
					s, result.n_own, result.n_shared, result.n_bound, result.n_tri, result.n_vert_total);
	bsp_sync(); 

	//Lets push the data to the corresponding processors!

	if (s == 0) {
		/* Keep track of a counter per proccesor, so we know what index we should put items */
		int vert_cntr[p], tri_cntr[p];
		memset(vert_cntr, 0, sizeof(int) * p);
		memset(tri_cntr, 0, sizeof(int) * p);
		//Push triangles
		for (int i = 0; i < mesh->n_tri; i++) 
			bsp_put(mesh->p[i], &mesh->t[i], result.t, (tri_cntr[mesh->p[i]]++)*SZTRI, SZTRI);
		
		//First push shared vertices (not bound)
		for (int v = 0; v < mesh->n_vert; v++)
			if (!mesh->b[v] && proc_count(vert_proc[v]) > 1) { //Shared vertex
				int proc = -1;
				while ((proc = proc_next(vert_proc[v], proc+1)) != -1) 
				{
					/* All the processors sharing vertex v need:
					 * - The coordinates of v
					 * - Boundary vertex
					 * - The global index of v (given by v itself, confusing huh ;-))
					 */
					bsp_put(proc, &mesh->x[v],	 result.x, vert_cntr[proc] * SZDBL, SZDBL);
					bsp_put(proc, &mesh->y[v],	 result.y, vert_cntr[proc] * SZDBL, SZDBL);
					bsp_put(proc, &vert_proc[v], result.p_shared, vert_cntr[proc] * SZPROCSET, SZPROCSET);
					bsp_put(proc, &v,            result.i_glob, vert_cntr[proc] * SZINT, SZINT);
					vert_cntr[proc]++;
			}
		}

		//Push non-shared vertices
		for (int v = 0; v < mesh->n_vert; v++)
			if (proc_count(vert_proc[v]) == 0) { //Non-shared vertex
				int proc = proc_next(vert_proc[v], -1);
				bsp_put(proc, &mesh->x[v],	 result.x, vert_cntr[proc] * SZDBL, SZDBL);
				bsp_put(proc, &mesh->y[v],	 result.y, vert_cntr[proc] * SZDBL, SZDBL);
				bsp_put(proc, &v,            result.i_glob, vert_cntr[proc] * SZINT, SZINT);
				vert_cntr[proc]++;
			}

		//Push boundary vertices
		for (int v = 0; v < mesh->n_vert; v++)
			if (mesh->b[v]) { //Boundary vertex
				int proc = -1;
				while ((proc = proc_next(vert_proc[v], proc+1)) != -1) 
				{
					bsp_put(proc, &mesh->x[v],	 result.x, vert_cntr[proc] * SZDBL, SZDBL);
					bsp_put(proc, &mesh->y[v],	 result.y, vert_cntr[proc] * SZDBL, SZDBL);
					bsp_put(proc, &v,            result.i_glob, vert_cntr[proc] * SZINT, SZINT);
					vert_cntr[proc]++;
				}
			}
	}
	bsp_sync();
	for (int i = 0; i < result.n_tri; i++)
		fprintf(stderr, "%d : (%d, %d, %d)\n", s, result.t[i][0], result.t[i][1], result.t[i][2]);
	/* Every processor has all the elements. Now convert global indexing to local indexing */
	bsp_pop_reg(result.x);
	bsp_pop_reg(result.y);
	bsp_pop_reg(result.i_glob);
	bsp_pop_reg(result.p_shared);
	
	result.glob2local = malloc(sizeof(int) * result.n_vert_total);
	for (int i = 0; i < result.n_vert_total; i++)
		result.glob2local[i] = -1;
	for (int i = 0; i < result.n_vert; i++) {
		fprintf(stderr, "%d : (%3f, %3f) \t loc:%d glob:%d\n",
				             s, result.x[i], result.y[i], i, result.i_glob[i]);
		result.glob2local[result.i_glob[i]] = i;
	}
	/* Convert triangles to local indices */
	for (int i =0; i < result.n_tri; i++)
		for (int v=0; v < 3; v++)
			result.t[i][v] = result.glob2local[result.t[i][v]];
	
	/* Memory structure is now all set, lets produce the FEM matrix */
	result.mat = gen_fem_mat(result.x, result.y, result.dof,
			                     result.t, result.n_tri);

	/* We are DONE! */
	return result;
}
//Symmetric sparse matrix vector multiplication v = Au
void ssmv(double * v, matrix_s * mat,  double * u) {
	for (int k = 0; k < mat->n; k++)
		v[k] = 0;

	for (int k = 0; k < mat->nz; k++) {
		int i = mat->I[k];
		int j = mat->J[k];
		v[i] += mat->val[k] * u[j];
		if (i != j)
			v[j] += mat->val[k] * u[i];
	}
}

//Calculate v = Au using our FEM-stuff
void bsp_fem_mv(int s, bsp_fem_data * fem, double * v, double * u) {
	int tagsz= SZINT;
	bsp_set_tagsize(&tagsz);
	bsp_sync();

	ssmv(v, &fem->mat, u);
	/*
	 * For vertices shared by other processors we must communicate.
	 * We therefore send the local value to all the shared processors with
	 * the global index as tag
	 */
	for (int k = 0; k < fem->n_shared; k++) {
		int proc = -1;
		while ((proc = proc_next(fem->p_shared[k], proc+1)) != -1) 
			if (proc != s) 
				bsp_send(proc, &fem->i_glob[k], &v[k], SZDBL);
	}
	bsp_sync();
	int status, tag;
	double tmp_val;
  bsp_get_tag(&status, &tag);
	while (status != -1) { //Process all the messages
		bsp_move(&tmp_val, SZDBL);
		//Tag holds the global index, tmp_val the value
		v[fem->glob2local[tag]] += tmp_val;
		bsp_get_tag(&status, &tag);
	}
	//Done!
}

double bsp_fem_ip(int p, int s,  //Processor information  
    double * x, double * y,      //Vector data
    bsp_fem_data * fem)          //Distribution information
{
  /*
   * Every processor calculates its own local inner product.
   * It then puts this value in all the other processors.
   * After which every processor calculates the total inner
   * product.
	 *
	 * We explicitely need to exclude vertices we shared, but that
	 * we do not own.
   */
  double inprod_local;
  double * inprod;
  inprod = vecallocd(p);
  bsp_push_reg(inprod, p * SZDBL);
  bsp_sync();

  //Calculate local inproduct
  inprod_local = 0;
	for (int i = 0; i < fem->n_shared; i++)
		if (s == proc_owner(fem->p_shared[i])) //We own this shared vertex
			inprod_local += fem->x[i] * fem->y[i];

	for (int i = fem->n_shared; i < fem->n_shared + fem->n_own; i++) //Loop over owned vertices
		inprod_local += fem->x[i] * fem->y[i];

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
	bsp_sync();

  return inprod_local;
} 

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
}

void print_fem_data(int s, bsp_fem_data result) {
	printf("%d\n"
			    "\tn_own:%d\n"
					"\tn_shared:%d\n"
					"\tn_bound:%d\n"
					"\tdof:%d\n"
					"\tn_tri:%d\n"
					"\tn_vert_total:%d\n",
					s, result.n_own, result.n_shared, result.n_bound, result.dof, result.n_tri, result.n_vert_total);
  //Print vertices
	for (int i = 0; i < result.n_vert; i++) 
		fprintf(stderr, "%d : (%3f, %3f) \t loc:%d glob:%d\n",
				             s, result.x[i], result.y[i], i, result.i_glob[i]);
	//Print triangles
	for (int i = 0; i < result.n_tri; i++)
		fprintf(stderr, "%d : (%d, %d, %d)\n", s, result.t[i][0], result.t[i][1], result.t[i][2]);
}

#define DEBUG(...) { if ( use_debug && s == 0) printf(__VA_ARGS__); }

void bspfem() {
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
	bsp_fem_data fem; //Stores the `local' FEM data
	mesh_dist mesh_total; //Processor 0 uses this

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

	//TODO: Load mesh_from file
	if (s == 0) {
		FILE * mesh_buf = fopen(meshbuffer, "r"); //TODO: FILE --> FILENAME
		mesh_total = readfrommeshfile(mesh_buf);
		fclose(mesh_buf);
		if (use_debug) 
			print_mesh_dist(mesh_total);
	}
	fem = bsp_fem_init(s,p, &mesh_total);

	if (use_debug)
		print_fem_data(s, fem);
	bsp_sync();
	/* Matrix, vector dist and b is loaded */
	/* Loading phase done, bsp_sync, add timer?*/

	/* Initalize bsp_mv */
	/*
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
	*/
	/* All initalization is done, start the main loop */
	/*
  while( rho > eps*eps*normbsq && k < kmax) {
    double beta, gamma, alpha;
    DEBUG("Iteration %d, rho = %g\n", k + 1, rho);
    if( k > 0) {
      beta = rho/rho_old;
      // u \gets r + beta u 
			
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
	*/
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
	*/
  bsp_end();
}


int main(int argc, char *argv[]) {
  bsp_init(bspfem,argc, argv); //Execute after we start with initalization
	P = bsp_nprocs();
	if(!(argc == 2 || argc == 3 || argc == 6 || argc == 8))
	{
		fprintf(stderr, "Invalid arguments given\n"
				            "\tUsage: %s mesh file [P [b_mode kmax eps [use_time use_debug]]]\n"
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

	meshbuffer  = malloc(BUFFERSTRSIZE);
	snprintf( meshbuffer, BUFFERSTRSIZE, "%s-P%d", argv[1], P);
	bspfem();
	free(meshbuffer);
  exit(0);
}
