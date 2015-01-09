#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "edge.h"
#include "tritype.h"
#include "helper.h"
#include "tree.h"
#include "types.h"

typedef struct i2v {
  int npoints;
  int *points;
  int index;
} i2v;

int get_index_from_vertex( i2v *p, int i) {
  assert( i < p->npoints);
  if( p->points[i] == -1) {
    p->points[i] = p->index++;
  }
  return p->points[i];
}

elnode *generate_elnode( workspace *w) {
  assert( w->is_conform);

  elnode *eln = malloc( sizeof( elnode));
  int bn = w->b->n;
  eln->e = malloc( w->nleaves * bn * sizeof( int));
  for( int i = 0; i < w->nleaves; i++) {
    for( int j = 0; j < w->b->n; j++) {
      eln->e[w->b->n*i + j] = -1;
    }
  }
  int points[w->npoints];

  for( int n = 0; n < w->npoints; n++) {
    points[n] = -1;
  }

  i2v index_helper = { .points = points, .npoints = w->npoints, .index = 0};
  for( int l = 0; l < w->nleaves; l++) {
    int nn = find_n( w->leaves[l]->r);
    int dof = (nn+1)*(nn+2)/2;
    assert( dof <= bn);
    tri *t = w->leaves[l]->tt;
    for( int n = 0; n < 3; n++) {
      if( w->boundary_points[t->p[n]]) {
        continue;
      }
      eln->e[l*bn+n] = get_index_from_vertex( &index_helper, t->p[n]);
      fprintf(stderr, "vert node index %i given on index %i,%i\n", eln->e[l*bn+n], l, 0);
    }
    if( nn >= 2) {
      for( int k = 1; k <= nn; k++) {
        //add face nodes
        if( k >= 3) {
          for( int r1 = 0; r1 <= k-3; r1++) {
            eln->e[l*bn+(k+2)*(k+1)/2 + r1 - (k-2)] = index_helper.index++;
            fprintf(stderr,"face node index %i given on index %i,%i (k=%i)\n", eln->e[l*bn+(k+2)*(k+1)/2 + r1 - (k-2)], l, (k+2)*(k+1)/2 + r1 - (k-2), k);
          }
        }
        if( k < nn) {
          //add edge nodes
          for( int n = 0; n < 3; n++) {
            tree *neighbour = edge_get( w->edges, t->p[(n+1)%3], t->p[n]);
            int neighbour_index = -1;
            int neighbour_edge = -1;

            //find leaf across edge (n, (n+1) mod 3)
            for( int m = 0; m < w->nleaves; m++) {
              if( neighbour_edge > -1) break;
              tri *o_t = w->leaves[m]->tt;
              for( int M = 0; M < 3; M++) {
                if( o_t->p[(M+1)%3] == t->p[n] && o_t->p[M] == t->p[(n+1)%3]) {
                  //found neighbour
                  neighbour = w->leaves[m];
                  neighbour_index = m;
                  neighbour_edge = M;
                  break;
                }
              }
            }

            //we're on the boundary
            if( neighbour == NULL) {
              continue;
              eln->e[l*bn+(k+2)*(k+1)/2 + n] = index_helper.index++;
              fprintf(stderr, "bndr node index %i given on index %i,%i\n", index_helper.index-1, l, (k+2)*(k+1)/2 + n);
            } 

            //degree of neighbour is not high enough
            if( find_n( neighbour->r) < k) {
              fprintf(stderr,"No DOF given\n");
              continue;
            }

            assert( neighbour_index > -1);

            if( neighbour_index > l) {
              eln->e[l*bn+(k+2)*(k+1)/2 + n] = index_helper.index++;
              fprintf(stderr,"edge node index %i given on index %i,%i (k=%i)\n", index_helper.index-1, l, (k+2)*(k+1)/2 + n, k);
            } else {
              eln->e[l*bn+(k+2)*(k+1)/2 + n] = eln->e[neighbour_index * bn + (k+2)*(k+1)/2 + neighbour_edge];
              fprintf(stderr,"edge node index %i taken on index %i,%i (k=%i)\n", eln->e[neighbour_index * bn + (k+2)*(k+1)/2 + neighbour_edge], l, (k+2)*(k+1)/2 + n, k);
            }
          }
        }
      }
    }
  }
  eln->M = index_helper.index;
  return eln;
}

fem_system *generate_fem_system( workspace *w, elnode *eln) {
  int n = w->b->n;

  fem_system *sys = malloc( sizeof( fem_system));
  sys->A = malloc( eln->M * eln->M * sizeof( double));
  sys->f = malloc( eln->M * sizeof( double));
  sys->M = eln->M;

  for( int i = 0; i < sys->M; i++) {
    for( int j = 0; j < sys->M; j++) {
      sys->A[i*sys->M + j] = 0.0;
    }
    sys->f[i] = 0.0;
  }

  for( int l = 0; l < w->nleaves; l++) {
    tri *t = w->leaves[l]->tt;
    point a = w->points[t->p[0]];
    point b = w->points[t->p[1]];
    point c = w->points[t->p[2]];
    double f11 = b.x - a.x;
    double f12 = c.x - a.x;
    double f21 = b.y - a.y;
    double f22 = c.y - a.y;
    double D = f11*f22 - f12*f21;
    double E1 = f12*f12 + f22*f22;
    double E2 = f11*f12 + f21*f22;
    double E3 = f11*f11 + f21*f21;

    int nn = find_n( w->leaves[l]->r);
    int dof = (nn+1)*(nn+2)/2;

    for( int r = 0; r < dof; r++) {
      int i = eln->e[l*n + r];
      if( i == -1) continue;
      for( int s = 0; s <= r; s++) {
        int j = eln->e[l*n + s];
        if( j == -1) continue;
        double krs = 1/fabs( D) * (E1*w->b->A[tritype_to_rtype(t)*n*n + r*n+s] 
                                 + E2*w->b->B[tritype_to_rtype(t)*n*n + r*n+s] 
                                 + E3*w->b->C[tritype_to_rtype(t)*n*n + r*n+s]);
        sys->A[i*sys->M + j] += krs;
        if( j != i) {
          sys->A[j*sys->M + i] += krs;
        }
      }

      sys->f[i] += fabs(D) * w->b->f[tritype_to_rtype(t)*n + r];
    }
  }

  return sys;
}
