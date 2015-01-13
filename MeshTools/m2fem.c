#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <libgen.h> 
#include <cblas.h>
#include "mesh.h"
#include "bspedupack.h"
#include "fem.h"


char * mesh_file;
FILE * output_mat;
FILE * output_rhs;
void  hack_fem_mat() {
	FILE * mesh = fopen(mesh_file, "r");
	mesh_dist mesh_total = readfrommeshfile(mesh);
	mesh_total.P = 1;
	mesh_total.p = calloc(sizeof(int), mesh_total.n_tri);
	bsp_fem_data hack = bsp_fem_init(0,1,&mesh_total);
	mat2mtx(output_mat, &hack.mat);
	//Output RHS in corresponding file format
	fprintf(output_rhs, "%d\n", hack.dof);
	int cntr = 1;
	for (int i = 0; i < hack.dof; i++)
		fprintf(output_rhs, "%d %lf\n", cntr++, hack.rhs[i]);
}

void bspexportmat() {
  bsp_begin(1);
	hack_fem_mat();
	bsp_end();
}
int main(int argc, char *argv[]) {
  bsp_init(bspexportmat,argc, argv); //Execute after we start with initalization
	if (argc != 2)
	{
		fprintf(stderr, "Invalid argument, give mesh file\n");
		exit(1);
	}
	mesh_file = argv[1];
	output_mat = stdout;
	output_rhs = stderr;
	bspexportmat();
  exit(0);
}
