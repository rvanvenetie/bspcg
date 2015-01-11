#ifndef MESH_H
#define MESH_H
typedef struct {
  int n, m, nz, lennz;
  int *I, *J;
  double *val;
} matrix_s;

typedef struct {
  int lenverts, lentris;
  int nverts, ntris, P;
  double *x, *y; //x and y coords
  int *b; //boundary
  int *t; //tris indices, len = ntris*3
  int *d; //distribution
  matrix_s *hypergraph;
} mesh_s;

//write mesh to mesh file, either distributed or not
int write2meshfile( FILE *fp, mesh_s *mesh, int distributed);

//read mesh file into struct
mesh_s *readfrommeshfile( FILE *fp);

//write mesh info mtx file (to go into Mondriaan)
int write2mtx( FILE *fp, mesh_s *mesh);

//read Mondriaan output and combine with existing mesh struct
int readfromvfile( FILE *fp, mesh_s *mesh);
#endif
