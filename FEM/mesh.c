#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mesh.h"

/**
 * Workflow:
 *  1. Create some mesh, load it into the struct
 *  2. Save the hypergraph to file using write2mtx()
 *  3. Run Mondriaan on this file to get the triangle distribution
 *  4. Load its result using readfromvfile()
 *  5. Write the resulting mesh and its distribution to file using write2meshfile()
 */

matrix_s mat_init(int n, int m) {
	matrix_s mat;
  mat.nz = 0;
  mat.n = n;
  mat.m = m;
  mat.lennz = 2;
  mat.I = malloc( 2* sizeof( int));
  mat.J = malloc( 2* sizeof( int));
  mat.val = malloc( 2* sizeof( double));
  return mat;
}

matrix_s * mat_create( int n, int m) {
  matrix_s *mat = malloc( sizeof( matrix_s));
	*mat = mat_init(n,m);
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

int write2meshfile( FILE *fp, mesh_s *mesh, int distributed) {
  if( distributed)
    fprintf( fp, "%%distributed\n");
  else
    fprintf( fp, "%%sequential\n");
  fprintf( fp, "%d %d", mesh->nverts, mesh->ntris);
  if( distributed)
    fprintf( fp, " %d", mesh->P);
  fprintf( fp, "\n");

  for( int i = 0; i < mesh->nverts; i++) {
    fprintf( fp, "%g %g %d\n", mesh->x[i], mesh->y[i], mesh->b[i]);
  }
  for( int i = 0; i < mesh->ntris; i++) {
    fprintf( fp, "%d %d %d", mesh->t[3*i+0], mesh->t[3*i+1], mesh->t[3*i+2]);
    if( distributed)
      fprintf( fp, " %d", mesh->d[i]);
    fprintf( fp, "\n");
  }

  return 0;
}

mesh_s *readfrommeshfile( FILE *fp) {
  int distributed = 0;
  mesh_s *mesh = malloc( sizeof( mesh_s));
  char buf[256];
  fscanf( fp, "%%%s\n", buf);
  if( strstr( buf, "distributed") != NULL)
    distributed = 1;
  fscanf( fp, "%d %d", &mesh->nverts, &mesh->ntris);
  if( distributed)
    fscanf( fp, " %d", &mesh->P);
  fscanf( fp, "\n");
  mesh->x = malloc( mesh->nverts * sizeof( double));
  mesh->y = malloc( mesh->nverts * sizeof( double));
  mesh->b = malloc( mesh->nverts * sizeof( int));

  mesh->t = malloc( 3 * mesh->ntris * sizeof( int));
  mesh->d = malloc( mesh->ntris * sizeof( int));
  
  for( int i = 0; i < mesh->nverts; i++) {
    fscanf( fp, "%lg %lg %d\n", &mesh->x[i], &mesh->y[i], &mesh->b[i]);
  }

  for( int i = 0; i < mesh->ntris; i++) {
    fscanf( fp, "%d %d %d", &mesh->t[3*i+0], &mesh->t[3*i+1], &mesh->t[3*i+2]);
    if( distributed)
      fscanf( fp, " %d", &mesh->d[i]);
    fscanf( fp, "\n");
  }

  return mesh;
}

int write2mtx( FILE *fp, mesh_s *mesh) {
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

void create_hypergraph_from_mesh( mesh_s *mesh) {
  matrix_s *hypergraph = mat_create( mesh->ntris, mesh->nverts);
  for( int i = 0; i < mesh->nverts; i++) {
    for( int j = 0; j < mesh->ntris; j++) {
      for( int k = 0; k < 3; k++) {
        if( mesh->t[3*j+k] == i)
          mat_append( hypergraph, i, j, 1);
      }
    }
  }

  mesh->hypergraph = hypergraph;
}

#ifdef FOKJOE
int main( void) {
  FILE *mfile = fopen( "lshaped.m", "r");
  FILE *mtxfile = fopen( "lshaped.mtx", "w");
  FILE *vfile = fopen( "lshaped.mtx-v2", "r");
  FILE *dmfile = fopen( "lshaped.dm", "w");

  mesh_s *mesh = readfrommeshfile( mfile);
  create_hypergraph_from_mesh( mesh);
  write2mtx( mtxfile, mesh);

  //now run Mondriaan on this mtx file
  exit(0);

  mesh->P = 2;
  readfromvfile( vfile, mesh);

  write2meshfile( dmfile, mesh, 1); //write distributed m file
  mesh_s *new_mesh = readfrommeshfile( dmfile);
  return 0;
}
#endif
