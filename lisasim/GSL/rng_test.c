#include <stdio.h>
#include "gsl_rng.h"
#include <sys/time.h>

int main (void) {
  gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlux);
  struct timeval tv;

  int i, n = 10;

  gettimeofday(&tv,0);
  printf("seconds is: %ld, microseconds is: %ld\n",tv.tv_sec,tv.tv_usec);
  printf("seed is: %ld\n",tv.tv_sec+tv.tv_usec);
  
  gsl_rng_set (r, tv.tv_usec); 

  for (i = 0; i < n; i++) {
    double u = gsl_rng_uniform (r);
    printf ("%.16g\n", u);
  }

  gsl_rng_free (r);

  return 0;
}
