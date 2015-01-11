typedef struct {
  int n, nz, nzA;
  int *inc;
  int nrows, ncols, *rowindex, *colindex;
  double *val;
} distributed_matrix;

typedef struct {
  int n, nv;
  int *vindex;
} vector_distribution;

distributed_matrix load_symm_distributed_matrix_from_file( const char *filename, const int p, const int s);
vector_distribution load_vector_distribution_from_file( const char *disfile, const char *valfile, const int p, const int s, double **pval);
