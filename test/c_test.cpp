#include <stdio.h>

struct c_struct {
  int i, j;
};

extern "C" void my_c (c_struct& cs) {

  printf ("%d  %d\n", cs.i, cs.j);
  cs.j = -72;

} 
