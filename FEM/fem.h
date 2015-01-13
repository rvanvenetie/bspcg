#ifndef FEM_H
#define FEM_H

#define SUPER_DEBUG 0
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

	/* RHS of the equation */
	double * rhs;
} bsp_fem_data;



void gen_element_matrix(double * res, double *x, double *y, triangle t, int dof, double * rhs);
matrix_s gen_fem_mat(bsp_fem_data * fem, double *x, double *y, int dof,
		                 triangle * t, int n_tri, double * rhs);

int proc_count(proc_set set);
void proc_add(proc_set * set, int proc);
void proc_remove(proc_set * set, int proc);
int proc_test(proc_set set, int proc);
int proc_next(proc_set set, int proc);
int proc_owner(proc_set set) ;

void fem_data_free(bsp_fem_data * fem);
bsp_fem_data bsp_fem_init(int s, int p, mesh_dist * mesh);

void bsp_fem_mv(int s,mesh_dist * mesh_total,  bsp_fem_data * fem, double * v, double * u);
void bsp_fem_shared_dof_sum(int s, bsp_fem_data * fem, double * v);
void print_fem_vect(int s, char * name, bsp_fem_data * fem, mesh_dist * mesh_total, double * x);
double bsp_fem_ip(int p, int s, double * x, double * y,     bsp_fem_data * fem);
#endif
