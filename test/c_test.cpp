#include <stdio.h>

extern "C" void c_test (int& i, double& r, int* i1, double* r1) {

  printf ("%d  %d  %d\n", i, i1[0], i1[1]);
  printf ("%f  %f  %f\n", r, r1[0], r1[1]);

} 
