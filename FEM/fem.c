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
#include "fem.h"

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
 * If a vertex is shared among multiple processors we have to decide who
 * is going to be the 'owner' of the vertex. We simply give the processor
 * with the smallest number the owner. 
 * 
 */

 

int icrs_remove_duplicates( matrix_icrs mat) {
  int n = mat.n, nz = mat.nz, *inc = mat.inc;
  double *a = mat.val;
  int pnrows = mat.pnc, pncols = mat.pnc, *prowindex = mat.pri, *pcolindex = mat.pci;
  int nnz = nz;
  int j = 0;
  for( int i = 1; i < nz; i++) {
    if (inc[i] == 0) {
      a[j] += a[i];
      nnz--;
    } else {
      a[++j] = a[i];
      inc[j] = inc[i];
    }
  }
  inc[nnz] = inc[nz];
  a[nnz] = 0;
  return nnz;
}

matrix_icrs coo2icrs( matrix_s mat, int remove_duplicates) {
  matrix_icrs imat;
  int *imatinc = vecalloci( mat.nz +1);
  int *imatja = vecalloci( mat.nz + 1);
  double *imatval = vecallocd( mat.nz +1);
  imat.inc = imatinc;
  imat.val = imatval;

  for( int i = 0; i < mat.nz; i++) {
    imat.inc[i] = mat.I[i];
    imatja[i] = mat.J[i];
    imat.val[i] = mat.val[i];
  }

  imat.n = mat.n;
  imat.nz = mat.nz;

  triple2icrs( mat.n, mat.nz, imat.inc, imatja, imat.val, &imat.pnr, &imat.pnc, &imat.pri, &imat.pci);

  vecfreei( imatja);

  if( remove_duplicates)
    mat.nz = icrs_remove_duplicates( imat);

  return imat;
}

void gen_element_matrix(double * res, double *x, double *y, triangle t, int dof, double * rhs) {
  //http://en.wikipedia.org/wiki/Stiffness_matrix (Die matrix D onderin)
  double D_mat[] = {x[t[2]] - x[t[1]], x[t[0]] - x[t[2]], x[t[1]] - x[t[0]],
	  							  y[t[2]] - y[t[1]], y[t[0]] - y[t[2]], y[t[1]] - y[t[0]]};

  //determinant affinetransformation times two.
  double det = fabs((x[t[1]] - x[t[0]]) * (y[t[2]] - y[t[0]]) - (x[t[2]] - x[t[0]]) * (y[t[1]] - y[t[0]]));
	double alpha = 1 / (2 * det);

	//Element matrix is now given by: alpha * D.T * D.
	//Only calculate lower tri part
	for (int i = 0; i < 3; i++) {
		if (t[i] < dof)  //Update rhs for nodal function i over this triangle
			rhs[t[i]] += det * 1.0/6.0;
		
		for (int j = 0; j <= i; j++)
		{ 
			//Row i of D.T times colum j of D
			res[i*3 + j] = alpha * ( D_mat[i] * D_mat[j] + D_mat[i + 3] * D_mat[j+3]);
		}
	}

}


matrix_s gen_fem_mat(bsp_fem_data * fem, double *x, double *y, int dof,
		                 triangle * t, int n_tri, double * rhs) {
	for (int i = 0; i < dof; i++)
		rhs[i] = 0.0;
	matrix_s result = mat_init(dof, dof);
	for (int k = 0; k < n_tri; k++) {
		if (SUPER_DEBUG) 
			printf("Element matrix triangle: (%d, %d, %d)\n",  fem->i_glob[t[k][0]],fem->i_glob[t[k][1]],fem->i_glob[t[k][2]]);
		//Calculate element matrix for triangle t[k]
		double el_mat[3*3];
		gen_element_matrix(el_mat, x,y,t[k],dof,  rhs);
		for (int li = 0; li < 3; li++)
			for (int lj = 0; lj <= li; lj++) //Local indices in element matrix
			{
				int i = t[k][li];
				int j = t[k][lj];
			
				if (SUPER_DEBUG)
				  printf("el_mat[%d,%d] = %.2f\n", li,lj, el_mat[li*3+lj]);

				//i and j hold the vertex index, only store if they account for a DOF
				if (i < dof && j < dof) {
					if (j <= i)
						mat_append(&result, i, j, el_mat[li*3 +lj]); 
					else
						mat_append(&result, j, i ,el_mat[li*3 +lj]);
				}
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
void fem_data_free(bsp_fem_data * fem) {
	free(fem->x);
	free(fem->y);
	free(fem->p_shared);
	free(fem->i_glob);
	free(fem->glob2local);
	free(fem->t);
	free(fem->rhs);
	mat_free(&fem->mat);
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
	if (SUPER_DEBUG)
		printf("%d\n"
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
		for (int i = 0; i < mesh->n_tri; i++) {
			if (SUPER_DEBUG)
				fprintf(stderr, "Processor %d gets triangle (%d,%d,%d) at local idx %d\n", mesh->p[i], mesh->t[i][0], mesh->t[i][1], mesh->t[i][2], tri_cntr[mesh->p[i]]);
			bsp_put(mesh->p[i], &mesh->t[i], result.t, (tri_cntr[mesh->p[i]]++)*SZTRI, SZTRI);
		}
		//First push shared vertices (not bound)
		for (int v = 0; v < mesh->n_vert; v++)
			if (!mesh->b[v] && proc_count(vert_proc[v]) > 1) { //Shared vertex
				int proc = -1;
				while ((proc = proc_next(vert_proc[v], proc+1)) != -1) 
				{
					if (SUPER_DEBUG)
						printf("Processor %d gets shared_vertex %d local idx %d\n", proc, v, vert_cntr[proc]);
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
			if (!mesh->b[v] && proc_count(vert_proc[v]) == 1) { //Non-shared vertex
				int proc = proc_next(vert_proc[v], -1);
				if (SUPER_DEBUG)
					printf("Processor %d gets own_vertex %d local idx %d\n", proc, v, vert_cntr[proc]);
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
					if (SUPER_DEBUG)
						printf("Processor %d gets bound_vertex %d local idx %d\n", proc, v, vert_cntr[proc]);
					bsp_put(proc, &mesh->x[v],	 result.x, vert_cntr[proc] * SZDBL, SZDBL);
					bsp_put(proc, &mesh->y[v],	 result.y, vert_cntr[proc] * SZDBL, SZDBL);
					bsp_put(proc, &v,            result.i_glob, vert_cntr[proc] * SZINT, SZINT);
					vert_cntr[proc]++;
				}
			}
	}
	bsp_sync();

	if (SUPER_DEBUG) 
		for (int i = 0; i < result.n_tri; i++)
			printf("%d : (%d, %d, %d)\n", s, result.t[i][0], result.t[i][1], result.t[i][2]);
	/* Every processor has all the elements. Now convert global indexing to local indexing */
	bsp_pop_reg(result.x);
	bsp_pop_reg(result.y);
	bsp_pop_reg(result.i_glob);
	bsp_pop_reg(result.p_shared);
	
	result.glob2local = malloc(sizeof(int) * result.n_vert_total);
	for (int i = 0; i < result.n_vert_total; i++)
		result.glob2local[i] = -1;
	for (int i = 0; i < result.n_vert; i++) {
		if (SUPER_DEBUG)
			printf("%d : (%3f, %3f) \t loc:%d glob:%d\n",
				           s, result.x[i], result.y[i], i, result.i_glob[i]);
		result.glob2local[result.i_glob[i]] = i;
	}
	/* Convert triangles to local indices */
	for (int i =0; i < result.n_tri; i++)
		for (int v=0; v < 3; v++)
			result.t[i][v] = result.glob2local[result.t[i][v]];
	
	/* Memory structure is now all set, lets produce the FEM matrix */
	result.rhs = malloc(sizeof(double) * result.dof);

	result.mat = gen_fem_mat(&result, result.x, result.y, result.dof,
			                     result.t, result.n_tri, result.rhs);
  
  //coo to icrs
  result.imat = coo2icrs( result.mat, 1);

	bsp_fem_shared_dof_sum(s, &result, result.rhs);
	/* We are DONE! */
	return result;
}

//Symmetric sparse matrix vector multiplication v = Au
void ssmv_coo(double * v, matrix_s * mat,  double * u) {
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

void ssmv_icrs( double *v, matrix_icrs *mat, double *u) {
	for (int k = 0; k < mat->n; k++)
		v[k] = 0;
  int j = mat->inc[0];
  int k = 0;

  for( int i = 0; i < mat->n; i++) {
    while( j < mat->n) {
      v[i] += mat->val[k]*u[j];
      if( i != j)
        v[j] += mat->val[k]*u[i];
      k++;
      j += mat->inc[k];
    }

    j -= mat->n;
	}
}

/*
 * Some vertices of the procesor are shared with others.
 * This function sums the values of the shared indices of a vector.
 * 
 * We do this by sendin the local value to all shared processors with the global index
 * as tag.
 */
void bsp_fem_shared_dof_sum(int s, bsp_fem_data * fem, double * v) {
	int tagsz= SZINT;
	bsp_set_tagsize(&tagsz);
	bsp_sync();

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
}

/*
 * Assemble fem vector on processor s = 0
 */
double * bsp_fem_ass_vect(int s, bsp_fem_data * fem, mesh_dist * mesh_total, double * x) {
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
	} else
		x_glob = NULL;
	return x_glob;
}

//Calculate v = Au using our FEM-stuff
void bsp_fem_mv(int s,mesh_dist * mesh_total,  bsp_fem_data * fem, double * v, double * u) {
  ssmv_icrs(v, &fem->imat, u);
	//ssmv_coo(v, &fem->mat, u);
	bsp_fem_shared_dof_sum(s, fem, v); //Sum the shared vertices

	if (SUPER_DEBUG) {
	  //print_fem_vect(s,"u in v=Au", fem, mesh_total, u);
	  //print_fem_vect(s,"v in v=Au", fem, mesh_total, v);
	}
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
			inprod_local += x[i] * y[i];

	for (int i = fem->n_shared; i < fem->dof; i++) //Loop over owned vertices
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
	bsp_sync();

  return inprod_local;
} 

