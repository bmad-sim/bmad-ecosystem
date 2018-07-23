/*========================================================================
  PISA  (www.tik.ee.ethz.ch/pisa/; www.lepp.cornell.edu/~ib38/apisa/)
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich

  Cornell University
  Ithaca, NY 14853
  ========================================================================
  NSGA2

  Implements data exchange trough files.
  
  file: nsga2_io.cpp
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by:  Stefan Bleuler, bleuler@tik.ee.ethz.ch
  rewritten by: Ivan Bazarov, bazarov@cornell.edu
  last change:  July 22, 2004
  ========================================================================
*/


#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>

#include "nsga2.hpp"

int read_pop(char *filename, pop *pp, int size, int dim, int con)
/* Reads individuals from file into pop */
{
  int entries = 0;
  char tag[4];
  FILE *fp;
  int result;
  
  assert(dim >= 0);
  assert(pp != NULL);
  
  fp = fopen(filename, "r");
  assert(fp != NULL);
  
  fscanf(fp, "%d", &entries);
  assert(entries == size * (dim + con + 1));
  
  for(int j = 0; j < size; j++) {
    /* reading index of individual */
    result = fscanf(fp, "%d", &(pp->ind_array[j]->index));
    
    for(int i = 0; i < dim; i++) {
      /* reading objective values of ind */
      result = fscanf(fp, "%le", &(pp->ind_array[j]->f[i]));
      if(result == EOF) { /* file not completely written */
	fclose(fp);
	return(1); /* signalling that reading failed */
      }
    }
    for(int i = 0; i < con; ++i) {
      /* reading constraint values of ind */
      result = fscanf(fp, "%le", &(pp->ind_array[j]->c[i]));
      if(result == EOF) { /* file not completely written */
	fclose(fp);
	return(1); /* signalling that reading failed */
      }
    }
  }
  
  /* after all data elements: "END" expected */
  fscanf(fp, "%s", tag);
  if (strcmp(tag, "END") != 0) {
    fclose(fp);
    return(1);  /* signalling that reading failed */
  }
  else { /* "END" ok */
    fclose(fp);
    /* delete file content if reading successful */
    fp = fopen(filename, "w");
    assert(fp != NULL);
    fprintf(fp, "0");
    fclose(fp);
    
    return(0);  /* signalling that reading was successful */
  }
}


void write_pop(char* filename, pop* pp, int size)
/* Writes a pop or PISA_to a given filename. */
{
  int i;
  FILE *fp;
  
  assert(0 <= size && size <= pp->size);
  
  fp = fopen(filename, "w");
  assert(fp != NULL);
  
  fprintf(fp, "%d\n", size); /* number of elements */
  
  for(i = 0; i < size; i++)
    fprintf(fp, "%d\n", pp->ind_array[i]->index);
     
  fprintf(fp, "END");
  fclose(fp);
}


int check_file(char* filename) {
  int control_element = 1;
  
  FILE *fp;
  
  fp = fopen(filename, "r");
  assert(fp != NULL);
  fscanf(fp, "%d", &control_element);
  fclose(fp);
  
  if(0 == control_element) return(0); /* file is ready for writing */
  else return(1); /* file is not ready for writing */
}
