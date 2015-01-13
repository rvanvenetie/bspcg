#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mesh.h"

int main( int argc, char **argv) {
  char mfn[1024], vfn[1024];
  int P;
  FILE *m, *v;
  if( argc == 1) {
    strcpy( mfn, "lshaped.m");
    P = 2;
  } else {
    if( argc < 2) {
      printf("Usage: %s [mfile P [vfile]]\n", argv[0]);
      exit(1);
    }
    strcpy( mfn, argv[1]);
    P = atoi( argv[2]);
  }
  if( argc > 3) {
    strcpy( vfn, argv[3]);
  } else {
    sprintf( vfn, "%stx-v%d", mfn, P);
  }

  m = fopen( mfn, "r");
  v = fopen( vfn, "r");

  if( !m|| !v) {
    printf("file not found!\n");
    exit(1);
  }

  mesh_dist mesh = readfrommeshfile( m);
  mesh.P = P;
  readfromvfile( v, &mesh);

  write2meshfile( stdout, &mesh, 1);
  return 0;
}
