#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    struct sparsematrix A;
    long P, k, q, n, i, j, t, nzp;
    int *Nprocs;

    printf("Test InitNprocs: ");
    P = 8; /* number of non-empty processors */
    k= 10*P;
    q= 4;
    n= q*k; /* n by n matrix A with nonzeros in positions
                A[i,j] with i mod q = j mod q = 0 */
    A.m = n;
    A.n = n;
    A.NrNzElts = k*k; 
    A.NrProcs = 2*P; /* we add P empty processors, to test some pathology */ 

    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.Pstart = (long *) malloc((2*P+1)* sizeof(long));
    Nprocs = (int *) malloc(n* sizeof(int));

    if (A.i == NULL || A.j  == NULL || A.Pstart == NULL || Nprocs == NULL) {
        printf("Error\n");
        exit(1);
    }

    /* Fill matrix with k*k nonzeros:
       k nonempty rows, each with k nonzeros */ 
    t= 0;
    for (i=0; i<k; i++) {
        for (j=0; j<k; j++) {
            A.i[t] = i*q;
            A.j[t] = j*q;
            t++;
        }
    }

    nzp = k*k/P;
    /* Procs 0, 2, ... , 2P-2 are empty.
       Procs 1, 3, ... , 2P-1 have nzp nonzeros each. */
    for (i=0; i<P; i++) {
        A.Pstart[2*i] = i*nzp;
        A.Pstart[2*i+1] = A.Pstart[2*i];
    }
    A.Pstart[2*P] = k*k;

    /* Initialise number of processors */
    if (!InitNprocs(&A, ROW, Nprocs)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result values and corresponding indices */
    for (j=0; j<n; j++) { 
        if(j%q == 0) {
            /* column j is nonempty and has P processors */
            if (Nprocs[j] != P) {
                printf("Error\n");
                exit(1);
            }
        } else {
            /* column j is empty and has 0 processors */ 
            if (Nprocs[j] != 0) {  
                printf("Error\n");
                exit(1);
            } 
        }
    }

    printf("OK\n");
    exit(0);

} /* end main */
