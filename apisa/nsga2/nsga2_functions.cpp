/*========================================================================
  APISA  (www.tik.ee.ethz.ch/pisa/; www.lepp.cornell.edu/~ib38/apisa/)
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich

  Cornell University
  Ithaca, NY 14853
  ========================================================================
  NSGA2
  
  Implements most functions.
  
  file: nsga2_functions.c
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
#include <cmath>

#include "nsga2.hpp"

/* common parameters */
int alpha;  /* number of individuals in initial population */
int mu;     /* number of individuals selected as parents */
int lambda; /* number of offspring individuals */
int dim;    /* number of objectives */
int con;    /* number of constraints */

/* local parameters from paramfile*/
int seed;   /* seed for random number generator */
int tournament;  /* parameter for tournament selection */

/* other variables */
char cfgfile[FILE_NAME_LENGTH];  /* 'cfg' file (common parameters) */
char inifile[FILE_NAME_LENGTH];  /* 'ini' file (initial population) */
char selfile[FILE_NAME_LENGTH];  /* 'sel' file (parents) */
char arcfile[FILE_NAME_LENGTH];  /* 'arc' file (archive) */
char varfile[FILE_NAME_LENGTH];  /* 'var' file (offspring) */

/* population containers */
pop* pp_all;
pop* pp_new;
pop* pp_sel;

/* NSGA2 internal global variables */
int* copies;
int** front;
double* dist;

bool verbose; /* whether to print out diagnostics messages */
char const* diag_file = "nsga2_diag.log";

/*-----------------------| initialization |------------------------------*/

void initialize(char *paramfile, char *filenamebase)
/* Performs the necessary initialization to start in state 0. */
{
  FILE *fp;
  int result;
  char str[CFG_ENTRY_LENGTH];
  
  /* reading parameter file with parameters for selection */
  fp = fopen(paramfile, "r");
  assert(fp != NULL);
  
  fscanf(fp, "%s", str);
  assert(strcmp(str, "seed") == 0);
  result = fscanf(fp, "%d", &seed);
  
  fscanf(fp, "%s", str);
  assert(strcmp(str, "tournament") == 0);
  fscanf(fp, "%d", &tournament); /* fscanf() returns EOF
				    if reading fails. */
  char ans[4];
  fscanf(fp, "%s", str);
  assert(strcmp(str, "verbose") == 0);
  result = fscanf(fp, "%s", ans);
  assert(!(strcmp(ans, "YES") && strcmp(ans, "NO")));
  verbose = !strcmp(ans, "YES") ? true : false;
  
  assert(result != EOF); /* no EOF, parameters correctly read */
  
  fclose(fp);
  
  srand(seed); /* seeding random number generator */
  
  sprintf(varfile, "%svar", filenamebase);
  sprintf(selfile, "%ssel", filenamebase);
  sprintf(cfgfile, "%scfg", filenamebase);
  sprintf(inifile, "%sini", filenamebase);
  sprintf(arcfile, "%sarc", filenamebase);
  
  /* reading cfg file with common configurations for both parts */
  fp = fopen(cfgfile, "r");
  assert(fp != NULL);
  
  fscanf(fp, "%s", str);
  assert(strcmp(str, "initial_population_size") == 0);
  fscanf(fp, "%d", &alpha);
  assert(alpha > 0);
  
  fscanf(fp, "%s", str);
  assert(strcmp(str, "parent_set_size") == 0);
  fscanf(fp, "%d", &mu);
  assert(mu > 0);
  
  fscanf(fp, "%s", str);
  assert(strcmp(str, "offspring_set_size") == 0);
  fscanf(fp, "%d", &lambda);
  assert(lambda > 0);
  
  fscanf(fp, "%s", str);
  assert(strcmp(str, "objectives") == 0);
  fscanf(fp, "%d", &dim);
  assert(dim > 0);

  fscanf(fp, "%s", str);
  assert(strcmp(str, "constraints") == 0);
  result = fscanf(fp, "%d", &con);
  
  assert(result != EOF); /* no EOF, 'dim' correctly read */
  
  fclose(fp);
  
  /* create individual and archive pop */
  pp_all = create_pop(alpha + lambda);
  pp_sel = create_pop(mu);    
}


/*-------------------| memory allocation functions |---------------------*/

void* chk_malloc(size_t size)
/* Wrapper function for malloc(). Checks for failed allocations. */
{
  void* return_value = malloc(size);
  if(return_value == NULL)
    PISA_ERROR("Selector: Out of memory.");
  return(return_value);
}

pop* create_pop(int maxsize)
/* Allocates memory for a population. */
{
  pop* pp;
  
  assert(maxsize >= 0);
  
  pp = (pop*) chk_malloc(sizeof(pop));
  pp->size = 0;
  pp->maxsize = maxsize;
  pp->ind_array = (ind**) chk_malloc(maxsize * sizeof(ind*));
  
  for(int i = 0; i < maxsize; i++) pp->ind_array[i] = NULL;
  
  return(pp);
}


ind* create_ind()
/* Allocates memory for one individual. */
{
  ind* p_ind;
  
  assert(dim > 0);
  assert(con >= 0);
  
  p_ind = (ind*) chk_malloc(sizeof(ind));
  
  p_ind->index = -1;
  p_ind->fitness = -1;
  p_ind->f = (double*) chk_malloc(dim * sizeof(double));
  p_ind->c = (double*) chk_malloc(con * sizeof(double));
  return(p_ind);
}


void free_memory()
/* Frees all memory. */
{
  free_pop(pp_sel);
  free_pop(pp_all);               
}


void free_pop(pop *pp)
/* Frees memory for given population. */
{
  assert(pp != NULL);
  
  free(pp->ind_array);
  free(pp);
}


void free_ind(ind *p_ind)
/* Frees memory for given individual. */
{
  assert(p_ind != NULL);
  
  free(p_ind->f);
  free(p_ind->c);
  free(p_ind);
}


/*-----------------------| selection functions|--------------------------*/

void selection()
{
  int i;
  int size;
  
  /* Join offspring individuals from variator to population */
  mergeOffspring();
  
  size = pp_all->size;
  
  /* Create internal data structures for selection process */
  /* Vectors */
  copies = (int*) chk_malloc(size * sizeof(int));
  dist = (double*) chk_malloc(size * sizeof(double));
  
  /* Matrices */
  front = (int**) chk_malloc(size * sizeof(int*));
  for(i = 0; i < size; i++) {
    front[i] = (int*) chk_malloc(size * sizeof(int));
  }
  
  /* Calculates NSGA2 fitness values for all individuals */
  calcFitnesses();

  /* Calculates distance cuboids */
  calcDistances();

  /* Performs environmental selection
     (truncates 'pp_all' to size 'alpha') */
  environmentalSelection();
  
  /* Performs mating selection
     (fills mating pool / offspring population pp_sel */
  matingSelection();
  
  /* Frees memory of internal data structures */    
  free(copies);
  free(dist);
  for(i = 0; i < size; i++) free(front[i]);
  free(front);
  
  return;
}


void mergeOffspring()
{
  int i;
  
  assert(pp_all->size + pp_new->size <= pp_all->maxsize);
  
  for(i = 0; i < pp_new->size; i++)
    pp_all->ind_array[pp_all->size + i] = pp_new->ind_array[i];
  
  pp_all->size += pp_new->size;
  
  free_pop(pp_new);
}


void calcFitnesses()
{
  int i, j, l;
  int size;
  int num;
  int *d;
  int *f;
  
  size = pp_all->size;
  d = (int*) chk_malloc(size * sizeof(int));
  f = (int*) chk_malloc(size * sizeof(int));
  
  /* initialize fitness and strength values */
  for (i = 0; i < size; i++) {
    pp_all->ind_array[i]->fitness = 0;
    d[i] = 1;
    f[i] = 1;
    copies[i] = 0;
  }
  
  /* calculate strength values */
  num = size;
  for(l = 0; l < size; l++) {
    /* find next front */
    for(i = 0; i < size; i++) {
      d[i] = 0;
      if(f[i] != 0) {
	for(j = 0; j < i && d[i] == 0; j++)
	  if(f[j] != 0) 
	    if(con_dominates(pp_all->ind_array[j], pp_all->ind_array[i])) d[i] = 1;
	for(j = i+1; j < size && d[i] == 0; j++)
	  if(f[j] != 0)
	    if(con_dominates(pp_all->ind_array[j], pp_all->ind_array[i])) d[i] = 1;
      }
    }
    
    /* extract front */
    for(i = 0; i < size; i++) {
      if(f[i] != 0 && d[i] == 0) {
	pp_all->ind_array[i]->fitness = l;
	f[i] = 0;
	num--;
	front[l][copies[l]] = i;
	copies[l] += 1;
      }
    }
    
    if (num == 0) break;
  }
  
  free(d);
  free(f);
  return;
}

void calcDistances()
{
  int i, j, l, d;
  int size = pp_all->size;
  double dmax = PISA_MAXDOUBLE / (dim + 1);
  
  for(i = 0; i < size; i++) dist[i] = 1;
  
  for(l = 0; l < size; l++) {
    for (d = 0; d < dim; d++) {
      /* sort accorting to d-th objective */
      for(i = 0; i < copies[l]; i++) {
	int min_index = -1;
	int min = i;
	for (j = i + 1; j < copies[l]; j++) {
	  if (pp_all->ind_array[front[l][j]]->f[d] <
	      pp_all->ind_array[front[l][min]]->f[d])
	    min = j;
	}
	min_index = front[l][min];
	front[l][min] = front[l][i];
	front[l][i] = min_index;
      }
      
      /* add distances */
      for (i = 0; i < copies[l]; i++) {
	if (i == 0 || i == copies[l] - 1)
	  dist[front[l][i]] += dmax;
	else {
	  dist[front[l][i]] +=
	    pp_all->ind_array[front[l][i+1]]->f[d] -
	    pp_all->ind_array[front[l][i-1]]->f[d];
	}
      }
    }
  }
}


void environmentalSelection()
{
  int i, j;
  int size = pp_all->size;
  
  static int nc = 0;
  ++nc;
  
  if(verbose) {
    char tmp [1024];
    sprintf(tmp, " nsga2: #%d population size (%d)\n", nc, size);
    msg_to_file(diag_file, tmp);
  }
  
  for(i = 0; i < size; i++)
    pp_all->ind_array[i]->fitness += 1.0 / dist[i];
  
  for(i = 0; i < alpha; i++) {
    ind *p_min;
    int min = i;
    for(j = i + 1; j < size; j++) {
      if(pp_all->ind_array[j]->fitness <
	 pp_all->ind_array[min]->fitness)
	min = j;
    }
    p_min = pp_all->ind_array[min];
    pp_all->ind_array[min] = pp_all->ind_array[i];
    pp_all->ind_array[i] = p_min;
  }
    
  for(i = alpha; i < size; i++) {
    free(pp_all->ind_array[i]);
    pp_all->ind_array[i] = NULL;
  }
  
  pp_all->size = alpha;
  
  return;
}


void matingSelection()
/* Fills mating pool 'pp_sel' */
{
  int i, j;
  
  for(i = 0; i < mu; i++) {
    int winner = irand(pp_all->size);
    
    for(j = 1; j < tournament; j++) {
      int opponent = irand(pp_all->size);
      if (pp_all->ind_array[opponent]->fitness
	  < pp_all->ind_array[winner]->fitness || winner == opponent) {
	winner = opponent;
      }
    }  
    pp_sel->ind_array[i] = pp_all->ind_array[winner];
  }
  pp_sel->size = mu;
}


void select_initial()
/* Performs initial selection. */
{
  selection();
  return;
}


void select_normal()
/* Performs normal selection.*/
{
  selection();
  return;
}


bool is_feasible(ind* i1)
/* Determines if an individual is feasible.
   Constraints are of the form g(x) >= 0. */
{
  for(int j = 0; j < con; ++j) if(i1->c[j] < 0) return false;
  return(true);
}


bool con_dominates(ind* i1, ind* i2)
/* Determines if one individual constrain-dominates another.
   Minimizing fitness values.
   Constraints are of the form g(x) >= 0. */
{
  bool i1_feasible = is_feasible(i1);
  bool i2_feasible = is_feasible(i2);
  
  if(i1_feasible)
    if(i2_feasible) return(dominates(i1,i2));
    else return(true);
  else if(i2_feasible) return(false);
  
  /* here if both individuals are infeasible */
  bool is_better = false;
  for(int j = 0; j < con; ++j) {
    if((i1->c[j] >= 0) && (i2->c[j] >= 0)) continue; 
    if(i1->c[j] < i2->c[j]) return(false);
    else if(i1->c[j] > i2->c[j]) is_better = true;
  }
  return(is_better);
}


bool dominates(ind* p_ind_a, ind* p_ind_b)
/* Determines if one individual dominates another.
   Minimizing fitness values. */
{
  bool a_is_worse = false;
  bool equal = true;
  
  for(int i = 0; i < dim && !a_is_worse; i++) {
    a_is_worse = p_ind_a->f[i] > p_ind_b->f[i];
    equal = (p_ind_a->f[i] == p_ind_b->f[i]) && equal;
  }
  
  return(!equal && !a_is_worse);
}


int irand(int range)
/* Generate a random integer. */
{
  int j;
  j=(int) ((double)range * (double) rand() / (RAND_MAX+1.0));
  return (j);
}


/*--------------------| data exchange functions |------------------------*/

int read_ini() {
  pp_new = create_pop(alpha);
  
  for(int i = 0; i < alpha; i++)
    pp_new->ind_array[i] = create_ind();
  pp_new->size = alpha;
  
  return(read_pop(inifile, pp_new, alpha, dim, con));                    
}


int read_var() {
  pp_new = create_pop(lambda);
  
  for(int i = 0; i < lambda; i++)
    pp_new->ind_array[i] = create_ind();
  
  pp_new->size = lambda;
  return (read_pop(varfile, pp_new, lambda, dim, con));
}


void write_sel() {
  write_pop(selfile, pp_sel, mu);
  return;
}


void write_arc() {
  write_pop(arcfile, pp_all, pp_all->size);
  return;
}


int check_sel() {
  return(check_file(selfile));
}


int check_arc() {
  return (check_file(arcfile));
}
