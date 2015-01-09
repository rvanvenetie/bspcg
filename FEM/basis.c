#include "basis.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

basis *basis_load( void) {
  basis *b = malloc( sizeof( basis));
  int n;

  char buffer[50];
  sprintf( buffer, "%idof/basis_0.txt", DOF);
  FILE *fp = fopen( buffer, "r");
  if( fp == NULL) perror( "error opening file");
  fscanf( fp, "%d\n", &n);
  b->n = n;

  b->A = malloc( 4*n * n * sizeof( double));
  b->B = malloc( 4*n * n * sizeof( double));
  b->C = malloc( 4*n * n * sizeof( double));

  b->C1 = malloc( 4*n * n * sizeof( double));
  b->C1inv = malloc( 4*n * n * sizeof( double));
  b->C2 = malloc( 4*n * n * sizeof( double));
  b->C2inv = malloc( 4*n * n * sizeof( double));

  b->f = malloc( 4*n * sizeof( double));
  for( int k = 0; k < 4; k++) {
    sprintf( buffer, "%idof/basis_%i.txt", DOF, k);
    fp = fopen( buffer, "r");
    if( fp == NULL) perror( "error opening file");
    fscanf( fp, "%d\n", &n);

    int i, j;
    for( i = 0; i < n; i++) {
      for( j = 0; j < n-1; j++) {
        fscanf( fp, "%lg,", &(b->A)[n*n*k + n*i+j]);
      }
      fscanf( fp, "%lg\n", &(b->A)[n*n*k + n*(i+1)-1]);
    }

    fscanf( fp, "\n");

    for( i = 0; i < n; i++) {
      for( j = 0; j < n-1; j++) {
        fscanf( fp, "%lg,", &(b->B)[n*n*k + n*i+j]);
      }
      fscanf( fp, "%lg\n", &(b->B)[n*n*k + n*(i+1)-1]);
    }

    fscanf( fp, "\n");

    for( i = 0; i < n; i++) {
      for( j = 0; j < n-1; j++) {
        fscanf( fp, "%lg,", &(b->C)[n*n*k + n*i+j]);
      }
      fscanf( fp, "%lg\n", &(b->C)[n*n*k + n*(i+1)-1]);
    }

    fclose( fp);

    // load rhs from file
    sprintf( buffer, "%idof/rhs_%i.txt", DOF, k);
    fp = fopen( buffer, "r");
    if( fp == NULL) perror( "error opening file");
    fscanf( fp, "%d\n", &n);
    assert( n == b->n);

    for( int i = 0; i < n; i++) {
      fscanf( fp, "%lg\n", &(b->f)[n*k + i]);
    }

    fclose( fp);

    // load C's from file
    sprintf( buffer, "%idof/C_%i.txt", DOF, k);
    fp = fopen( buffer, "r");
    if( fp == NULL) perror( "error opening file");
    fscanf( fp, "%d\n", &n);
    assert( n == b->n);

    for( i = 0; i < n; i++) {
      for( j = 0; j < n-1; j++) {
        fscanf( fp, "%lg,", &(b->C1)[n*n*k + n*i+j]);
      }
      fscanf( fp, "%lg\n", &(b->C1)[n*n*k + n*(i+1)-1]);
    }
    fscanf( fp, "\n");

    for( i = 0; i < n; i++) {
      for( j = 0; j < n-1; j++) {
        fscanf( fp, "%lg,", &(b->C1inv)[n*n*k + n*i+j]);
      }
      fscanf( fp, "%lg\n", &(b->C1inv)[n*n*k + n*(i+1)-1]);
    }
    fscanf( fp, "\n");

    for( i = 0; i < n; i++) {
      for( j = 0; j < n-1; j++) {
        fscanf( fp, "%lg,", &(b->C2)[n*n*k + n*i+j]);
      }
      fscanf( fp, "%lg\n", &(b->C2)[n*n*k + n*(i+1)-1]);
    }
    fclose( fp);

    for( i = 0; i < n; i++) {
      for( j = 0; j < n-1; j++) {
        fscanf( fp, "%lg,", &(b->C2inv)[n*n*k + n*i+j]);
      }
      fscanf( fp, "%lg\n", &(b->C2inv)[n*n*k + n*(i+1)-1]);
    }
    fclose( fp);
  }

  /*
  for( int k = 0; k < 4; k++) {
    for( int i = 0; i < n; i++) {
      for( int j = 0; j < n; j++) {
        printf("%g ", b->A[k*n*n + n*i + j]);
      }
    }
    printf("\n");
  }
  */

  return b;
}
