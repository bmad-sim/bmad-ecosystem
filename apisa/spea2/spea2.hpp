/*========================================================================
  APISA  (www.tik.ee.ethz.ch/pisa/; www.lepp.cornell.edu/~ib38/apisa/)
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich

  Cornell University
  Ithaca, NY 14853
  ========================================================================
  SPEA2 - Strength Pareto EA 2

  Implementation in C++ for the selector side.
  
  Header file.
  
  file: spea2.hpp
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by:  Stefan Bleuler, bleuler@tik.ee.ethz.ch
  rewritten by: Ivan Bazarov, bazarov@cornell.edu
  last change:  July 22, 2004
  ========================================================================
*/

#ifndef SPEA2_HPP
#define SPEA2_HPP

/*-----------------------| specify Operating System |------------------*/
/* necessary for wait() */

/* #define PISA_WIN */
#define PISA_UNIX

/*----------------------------| macro |----------------------------------*/

#define PISA_ERROR(x) fprintf(stderr, "\nError: " x "\n"), fflush(stderr), exit(EXIT_FAILURE)

/*---------------------------| constants |-------------------------------*/
#define FILE_NAME_LENGTH 128 /* maximal length of filenames */
#define CFG_ENTRY_LENGTH 128 /* maximal length of entries in cfg file */
#define PISA_MAXDOUBLE 1E99  /* Internal maximal value for double */

/*----------------------------| structs |--------------------------------*/

typedef struct ind_st  /* an individual */
{
  int index;
  double *f; /* objective vector */
  double* c; /* constraint vector */
  double fitness;
} ind;

typedef struct pop_st  /* a population */
{
  int size;
  int maxsize;
  ind** ind_array;
} pop;

/*-------------| functions for control flow (in spea2.cpp) |------------*/

void write_flag(char *filename, int flag);
int read_flag(char *filename);
void wait(double sec);

/*-------------| auxiliary functions |------------*/

int msg_to_file(const char *file, const char *str);

inline void error(const char* v) {
  fprintf(stderr, "%s\nprogram stops.\n", v);
  exit(1);
  return;
}

/*---------| initialization function (in spea2_functions.cpp) |---------*/

void initialize(char *paramfile, char *filenamebase);

/*--------| memory allocation functions (in spea2_functions.cpp) |------*/

void* chk_malloc(size_t size);
pop* create_pop(int size);
ind* create_ind();

void free_memory();
void free_pop(pop *pp);
void free_ind(ind *p_ind);

/*-----| functions implementing the selection (spea2_functions.cpp) |---*/

void selection();
void mergeOffspring();
void calcFitnesses();
void calcDistances();
int getNN(int index, int k);
double getNNd(int index, int k);
void environmentalSelection();
void truncate_nondominated();
void truncate_dominated();
void matingSelection();

void select_initial();
void select_normal();
bool is_feasible(ind* i1);
bool con_dominates(ind* i1, ind* i2);
bool dominates(ind *p_ind_a, ind *p_ind_b);
int is_equal(ind *p_ind_a, ind *p_ind_b);
double calcDistance(ind *p_ind_a, ind *p_ind_b);
int irand(int range);

/*--------------------| data exchange functions |------------------------*/

/* in spea2_functions.cpp */

int read_ini();
int read_var();
void write_sel();
void write_arc();
int check_sel();
int check_arc();

/* in spea2_io.cpp */

int read_pop(char *filename, pop *pp, int size, int dim, int con);
void write_pop(char *filename, pop *pp, int size);
int check_file(char *filename);

#endif /* SPEA2_HPP */
