#include <stdlib.h>
#include <stdio.h>
#include "mesh.h"

matrix_s * mat_create( int n, int m) {
  matrix_s *mat = malloc( sizeof( matrix_s));
  mat->nz = 0;
  mat->n = n;
  mat->m = m;
  mat->lennz = 2;
  mat->I = malloc( 2* sizeof( int));
  mat->J = malloc( 2* sizeof( int));
  mat->val = malloc( 2* sizeof( double));

  return mat;
}

int mat_append( matrix_s *mat, int i, int j, double val) {
  if( mat->lennz == mat->nz) {
    mat->I = realloc( mat->I, 2*mat->lennz * sizeof( int));
    mat->J = realloc( mat->J, 2*mat->lennz * sizeof( int));
    mat->val = realloc( mat->val, 2*mat->lennz * sizeof( double));
    mat->lennz *= 2;
  }
  mat->I[mat->nz] = i;
  mat->J[mat->nz] = j;
  mat->val[mat->nz] = val;
  return mat->nz++;
}

int write2mmfile( FILE *fp, mesh_s *mesh) {
  fprintf( fp, "%%%%MatrixMarket weightedmatrix coordinate pattern general\n");
  int nhypernets = mesh->nverts;
  int nhyperverts = mesh->ntris;
  fprintf( fp, "%d %d %d %d\n", nhypernets, nhyperverts, mesh->hypergraph->nz, 3);

  for( int i = 0; i < mesh->hypergraph->nz; i++) {
    fprintf( fp, "%d %d\n", mesh->hypergraph->I[i]+1, mesh->hypergraph->J[i]+1);
  }

  for( int i = 0; i < nhyperverts; i++) {
    fprintf( fp, "1\n"); //TODO whats the correct weight for this?
  }

  for( int i = 0; i < nhypernets; i++) {
    fprintf( fp, "%d\n", 1-(int)mesh->b[i]); //can sigma(e) be zero as well?
  }

  return 0;
}

int readfromvfile( FILE *fp, mesh_s *mesh) {
  int ncols, P, col, s;
  fscanf( fp, "%d %d\n", &ncols, &P);
  if( ncols != mesh->ntris) exit( 1);

  for( int i = 0; i < ncols; i++) {
    fscanf( fp, "%d %d\n", &col, &s);
    if( col != i+1) exit( 1);
    mesh->d[i] = s-1;
  }

  return 0;
}

int main( void) {
  const int nverts = 5;
  const int ntris = 4;
  const int P = 2; //number of processors
  double x[nverts] = {0, 1, 1, 0, 0.5};            //x-values for vertices
  double y[nverts] = {0, 0, 1, 1, 0.5};            //y-values for vertices
  int b[nverts] = {1, 1, 1, 1, 0};            //boolean boundary for vertices
  int t[ntris*3] = {4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0};     //triangles
  int d[ntris] = {-1};

  mesh_s mesh = {
    .nverts = nverts,
    .ntris = ntris,
    .x = x,
    .y = y,
    .b = b,
    .t = t,
    .P = P,
    .d = d
  };

  //create hypergraph matrix
  matrix_s *hypergraph = create( ntris, nverts);
  for( int i = 0; i < nverts; i++) {
    for( int j = 0; j < ntris; j++) {
      for( int k = 0; k < 3; k++) {
        if( t[3*j+k] == i)
          append( hypergraph, i, j, 1);
      }
    }
  }

  mesh.hypergraph = hypergraph;

  write2mmfile( stdout, &mesh);
  FILE *vfile = fopen( "bla.mtx-v2", "r");
  readfromvfile( vfile, &mesh);
  for( int i = 0; i < mesh.P; i++) {
    //printf("%d %d\n", i, mesh.d[i]);
  }
  return 0;
}
