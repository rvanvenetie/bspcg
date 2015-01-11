#include <stdlib.h>
#include <stdio.h>
#include "mesh.h"

int main( int argc, char **argv) {
  FILE *fp;
  if( argc == 1) {
    fp = fopen( "lshaped.m", "r");
  } else {
    fp = fopen( argv[1], "r");
  }

  if( !fp) {
    printf("file not found!\n");
    exit(1);
  }

  mesh_dist mesh = readfrommeshfile( fp);
  matrix_s *hg = create_hypergraph_from_mesh( &mesh);

  write2mtx( stdout, &mesh, hg);
  return 0;
}
