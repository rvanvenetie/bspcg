#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "bspedupack.h"

#define MODE GOLDBACH
#define BENCHMARK 0

#define TWIN 0
#define GOLDBACH 1
#define JUST_PRIMES 2
#define JUST_ONCE 1

/*
 * Usage: ./main <p> <n> to generate all primes below n with p nodes. In 
 * benchmark mode, only use first argument.
 */

typedef char * sieve_array;
#define PRIME 1
#define COMPOSITE 0

void print_sieve_primes(int length, int offset,  sieve_array sieve) {
  for (int i = 0; i < length; i++) {
    //printf("%d ", i + offset);
    if (sieve[i] == PRIME) {
      printf("%d\t", i + offset);
    }  
  }
}

void print_twin_primes(sieve_array sieve, int n) {
  for (int i = 2; i <= n - 2; i++)
    if ((sieve[i] == PRIME) && (sieve[i+2] == PRIME)) {
      printf("(%d,%d) ", i,i+2);
    }  
}

/*
 * Index in the sieve_array corresponds to the actual number
 * Array result must be initalized.
 * 
 * Does the sieve_array indexing for all primes < sqrt_N
 * This means that index 0 and 1 are undefined.
 */
void sieve_primes_sqrt(int n, sieve_array result) {
  int sqrtN = (int) sqrt(n);
  if (MODE == TWIN)
    sqrtN += 2;
    
  for (int i = 2; i <= sqrtN; i++) {
    if (result[i] == PRIME) {
      for (int j = i*i; j <= sqrtN; j += i) {
        result[j] = COMPOSITE;
      }
    }
  }
}

/*
 * Calculates the blocksize for processor s. Also returns the starting number + end number.
 * Distributes the numbers floor(sqrt(n))..n over p processors.
 */
void sieve_primes_dist(int p, int s, int n, int * start, int * end, int * blocksize) {
  int first_num  = (int) sqrt(n) + 1;
  int n_left = n - first_num + 1;
  *blocksize = (int) (n_left) / p;
  
  *start = (*blocksize) * s + first_num; 

  //The last block needs to do more work
  if (s == p - 1)  
    *blocksize += n_left - p * (*blocksize);
  else if (MODE == TWIN) //Let the other blocks calculate two more numbers 
    *blocksize += 2; 
    
  *end = *start + (*blocksize) - 1; 
}

/*
 * Cross out all multiples of primes in inital primes (the first 0 .. sqrt N).
 * Use easy block structure. Returns pointer to start of the block in all_primes.
 */
sieve_array sieve_primes_block(int p,int s,int n,sieve_array all_primes) {
  int start,end,blocksize;
  int sqrtN = (int) sqrt(n);
  
  sieve_primes_dist(p,s,n,&start,&end,&blocksize);
  sieve_array A = all_primes + start;
  //memset(A, PRIME, blocksize); //Initialize
  
  for (int i = 2; i <= sqrtN; i++) {
    if (all_primes[i] == PRIME) {
      int j = start;
    
      if (j % i) //j is not a multiple of i, take next multiple of i.
        j += i - j % i;
      if( j < i*i)
        j = i*i;
      
      for (; j <= end; j += i) {
          A[j-start] = COMPOSITE;
      }
    }
  }
  return A;  
}

void twin_primes_block(sieve_array primes, int offset, int size) {
  for (int i = 0; i < size - 2; i++) 
    if (!((primes[i + offset] == PRIME) && (primes[i + offset +2] == PRIME)))
      primes[i+offset] = COMPOSITE;  
}

int P,N; //Global stuff

void bspprimes() {
  int p,s,n;
  double time_start,time_init_primes,  time_initial, time_sieve, time_end;
  sieve_array block_primes, all_primes;
  
  bsp_begin(P);
  p= bsp_nprocs(); /* p = number of processors obtained */
  s= bsp_pid();    /* s = processor number */
    
  #if BENCHMARK == 1
  if (s == 0)
    N = 1;
  for (int i_bench = 0; i_bench < 8; i_bench++) {
    if (s == 0)
      N = N * 10;
  #endif

  n = N;
  bsp_sync();
  time_start = bsp_time();

  //First register n. Registration needs a sync step
  bsp_push_reg(&n,SZINT);
  bsp_sync();

  //All get the value of n from the first processor.
  bsp_get(0,&n,0,&n,SZINT);
  bsp_sync();

  //We now have the value of n
  //We do not need the linked version of n anymore, delete superstep
  bsp_pop_reg(&n);
  
  //Allocate memory array for all primes and bind register
  all_primes = malloc(n+1);
  memset(all_primes, PRIME, n+1);
  bsp_push_reg(all_primes, n+1);
  
  time_initial = bsp_time();
  
  //Let this processor get all the initial primes
  sieve_primes_sqrt(n, all_primes);
  time_init_primes = bsp_time();
  
  //Evenly distribute blocks across processors and let each processor
  //cross out the primes in his block.
  int start,end,blocksize;
  sieve_primes_dist(p,s,n,&start,&end,&blocksize);  
  block_primes = sieve_primes_block(p,s,n, all_primes);
  
  //Possibly calculate twin primes
  if (MODE == TWIN) {
    if (s == 0)
      twin_primes_block(all_primes, 0 , (int) sqrt(n) + 2);
  
    twin_primes_block(all_primes, start, blocksize);
  }
  time_sieve = bsp_time();
  
  /* Start communicating all primes */
  bsp_sync();
  
  #if JUST_ONCE == 1
  //Communicate primes in a block
  for (int t= 0; t < p; t++)
    bsp_put(t,block_primes,all_primes, start, blocksize);
  #else
  //Communicate primes separately
  for (int t = 0; t < p; t++) 
    for (int j = start; j <= end; j++)
      if (all_primes[j] == PRIME)
        bsp_put(t,all_primes + j, all_primes, j, 1); 
  #endif
  bsp_sync();
   
  time_end = bsp_time();
  /* Now every processor should have all the primes in all_primes */
  if (BENCHMARK && (s ==0)) {
    printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n",
       p,
       n,
       time_initial - time_start,
       time_init_primes - time_initial,
       time_sieve - time_init_primes,
       time_end - time_sieve,
       time_end - time_start);
       
  } else if (s == 0) {
    printf("***** TIMING REPORT N = %d FOR PROCESSOR %d OF  %d *****\n",n, s + 1, p);
    printf("Begin block: %d \t Block size: %d\n", start, blocksize);
    printf("Initial communication:   %f sec\n", time_initial - time_start);
    printf("Initial prime gathering: %f sec\n", time_init_primes - time_initial); 
    printf("Sieving prime block:     %f sec\n", time_sieve - time_init_primes);
    printf("Communication primes:    %f sec\n", time_end - time_sieve);
    printf("Total:                   %f sec\n", time_end - time_start);
    printf("***** END TIMING REPORT *****\n");
    //print_sieve_primes(n+1,0,all_primes);
  }
  
  /*
   * Goldbach mode? We now have all primes. Let us test the goldbach conjecture on
   * a dataset. We use cyclic distribution here, as we expect this to evenly
   * distribute the computing cost. 
   */
  if (MODE == GOLDBACH){
    printf("aasdf\n");
    double time_gold_start = bsp_time();
    int goldbach_fails;
    for (int i = 4 + s*2; i <= n; i += P*2) {
      goldbach_fails = 1;
      for (int j = 2; j <= i/2; j++)
        if ((all_primes[j] == PRIME) && (all_primes[i-j] == PRIME))
          goldbach_fails = 0;
      if( goldbach_fails) {
        printf("Goldbach for processor %d:\n  Failed %d\n", s,goldbach_fails);
      }
    }
    double time_gold_end = bsp_time();
    printf("Processor %d Golbach Time: %f\n", s, time_gold_end - time_gold_start);
  }
  bsp_pop_reg(all_primes);
  free(all_primes);
  #if BENCHMARK == 1
    bsp_sync();
  }

  #endif
  bsp_end();
}


int main(int argc, char *argv[]) {
  bsp_init(bspprimes,argc, argv); //Execute bspprimes after we start with initalization
  if (argc == 3) {
    P = atoi(argv[1]);
    N = atoi(argv[2]);
  } else if (BENCHMARK) {
    P = bsp_nprocs();
    N = 10000;
  } else {
    printf("How many processors do you want to use?\n"); 
    fflush(stdout);
    scanf("%d",&P);
    printf("Please enter n:\n"); fflush(stdout);
    scanf("%d",&N);
  }
  if(N<0)
    bsp_abort("Error in input: n is negative");
  if (P > bsp_nprocs()){
      printf("Sorry, not enough processors available.\n");
      exit(1);
  }
  bspprimes();
  exit(0);
}
