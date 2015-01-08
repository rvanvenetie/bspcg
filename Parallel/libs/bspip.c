#include "bspip.h"
#include "bspedupack.h"

/*
 * Edited version of bspip that uses vectors that are stored
 * according to a certain distribution (both must have the same).

 * Input:
 * p - number of processors.
 * s - current processor number.
 * x - locally stored vector elements of vector x
 * y - locally stored vector elements of vector y
 * nv - amount of vector elements locally stored.
 */
double bspip_dist(int p, int s,  //Processor information  
    double * x, double * y,      //Vector data
    int nv)                      //Distribution information
{
  /*
   * Every processor calculates its own local inner product.
   * It then puts this value in all the other processors.
   * After which every processor calculates the total inner
   * product
   */
  //Initialize array that holds all partial sums
  double inprod_local;
  double * inprod;
  inprod = vecallocd(p);
  bsp_push_reg(inprod, p * SZDBL);
  bsp_sync();

  //Calculate local inproduct
  inprod_local = 0;
  for (int i = 0; i < nv; i++)
    inprod_local += x[i] * y[i];

  //Put the local value in all the processors
  for (int t = 0; t < p; t++) 
    bsp_put(t, &inprod_local, inprod, s * SZDBL, SZDBL);

  bsp_sync();
  //Calculate the actual inner product
  inprod_local = 0;
  for (int t = 0; t < p; t++)
    inprod_local += inprod[t];

  //Free data stuff
  bsp_pop_reg(inprod);
  vecfreed(inprod);
	bsp_sync();

  return inprod_local;
} 
