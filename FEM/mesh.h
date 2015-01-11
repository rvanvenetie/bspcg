#ifndef MESH_H
#define MESH_H
typedef struct {
  int n, m, nz, lennz;
  int *I, *J;
  double *val;
} matrix_s;

typedef int triangle[3];
/* Data structure for the (distributed) triangulation */
typedef struct {
  /* Vertices */
  double *x, *y;//Coordinates
	int * b;      //Lies on boundary of domain?
  int n_vert;   //Amount of vertices
  
  /* Triangles */
  triangle * t;   //Triangles consist of three vertices
  int * p;      //Processor that is owner of this triangle
  int n_tri;    //Amount of triangles
  int P; //number of processors in distribution
} mesh_dist;

matrix_s mat_init(int n, int m);
int mat_append( matrix_s *mat, int i, int j, double val);

//write mesh to mesh file, either distributed or not
int write2meshfile( FILE *fp, mesh_dist *mesh, int distributed);
//read mesh file into struct
mesh_dist readfrommeshfile( FILE *fp);

//write mesh info mtx file (to go into Mondriaan)
int write2mtx( FILE *fp, mesh_dist *mesh, matrix_s *hypergraph);

//read Mondriaan output and combine with existing mesh struct
int readfromvfile( FILE *fp, mesh_dist *mesh);
matrix_s *create_hypergraph_from_mesh( mesh_dist *mesh);
#endif
