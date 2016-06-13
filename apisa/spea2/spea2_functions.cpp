/*========================================================================
  APISA  (www.tik.ee.ethz.ch/pisa/; www.lepp.cornell.edu/~ib38/apisa/)
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich

  Cornell University
  Ithaca, NY 14853
  ========================================================================
  SPEA2 - Strength Pareto EA 2
  
  Implements most functions.
  
  file: spea2_functions.cpp
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

#include "spea2.hpp"

/* common parameters */
int alpha;  /* number of individuals in initial population */
int mu;     /* number of individuals selected as parents */
int lambda; /* number of offspring individuals */
int dim;    /* number of objectives */
int con;    /* number of constraints */

/* local parameters from paramfile*/
int seed;   /* seed for random number generator */
int tournament;  /* parameter for tournament selection */
int kth; /* k-th neighbor in niching operator */
bool kth_is_sqrt; /* k-th is equal to sqrt of merged population if true */
bool verbose;

char const* diag_file = "spea2_diag.log";

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


/* SPEA2 internal global variables */
int* fitness_bucket;
int* fitness_bucket_mod;
int* copies;
int* old_index;
double** dist; /* distance between individuals in combined population */
int** NN;



/*-----------------------| initialization |------------------------------*/

void initialize(char* paramfile, char* filenamebase)
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
  fscanf(fp, "%d", &seed);
  
  fscanf(fp, "%s", str);
  assert(strcmp(str, "tournament") == 0);
  fscanf(fp, "%d", &tournament);

  char ans[5];
  fscanf(fp, "%s", str);
  assert(strcmp(str, "k_neighbor") == 0);
  fscanf(fp, "%s", ans);
  if(!strcmp(ans, "SQRT")) kth_is_sqrt = true;
  else kth_is_sqrt = false;
  
  if(!kth_is_sqrt) {
    kth = atoi(ans);
    assert(kth > 0);
  }

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
  void *return_value = malloc(size);
  if(return_value == NULL) PISA_ERROR("Selector: Out of memory.");
  return (return_value);
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
  /* Join offspring individuals from variator to population */
  mergeOffspring();
  
  int size = pp_all->size;

  /* set k-th neighbor number */
  if(kth_is_sqrt) kth = (int) sqrt(size);
  assert(kth > 0);

  /* Create internal data structures for selection process */
  /* Vectors */
  fitness_bucket = (int*) chk_malloc(size * size * sizeof(int));
  fitness_bucket_mod = (int*) chk_malloc(size * sizeof(int));
  copies = (int*) chk_malloc(size * sizeof(int));
  old_index = (int*) chk_malloc(size * sizeof(int));
  
  /* Matrices */
  dist = (double**) chk_malloc(size * sizeof(double*));
  NN = (int**) chk_malloc(size * sizeof(int*));
  for (int i = 0; i < size; i++) {
    dist[i] = (double*) chk_malloc(size * sizeof(double));	  
    NN[i] = (int*) chk_malloc(size * sizeof(int));
  }

  /* Calculates SPEA2 fitness values for all individuals */
  calcFitnesses();

  /* Calculates distance matrix dist[][] */
  calcDistances();
  
  /* Performs environmental selection
     (truncates 'pp_all' to size 'alpha') */
  environmentalSelection();

  /* Performs mating selection
     (fills mating pool / offspring population pp_sel */
  matingSelection();

  /* Frees memory of internal data structures */    
  free(fitness_bucket);
  free(fitness_bucket_mod);
  free(copies);
  free(old_index);
  for(int i = 0; i < size; i++) {
    free(dist[i]);
    free(NN[i]);
  }
  free(dist);
  free(NN);
  
  return;
}


void mergeOffspring()
{
  assert(pp_all->size + pp_new->size <= pp_all->maxsize);
  
  for(int i = 0; i < pp_new->size; i++)
    pp_all->ind_array[pp_all->size + i] = pp_new->ind_array[i];
  
  pp_all->size += pp_new->size;
  
  free_pop(pp_new);
}


void calcFitnesses()
{  
  int size = pp_all->size;
  int* strength = (int*) chk_malloc(size * sizeof(int));
  
  /* initialize fitness and strength values */
  for(int i = 0; i < size; i++) {
    pp_all->ind_array[i]->fitness = 0;
    strength[i] = 0;
    fitness_bucket[i] = 0;
    fitness_bucket_mod[i] = 0;	
    for(int j = 0; j < size; j++) fitness_bucket[i * size + j] = 0;
  }
    
  /* calculate strength values */
  for(int i = 0; i < size; i++)
    for(int j = 0; j < size; j++)
      if(con_dominates(pp_all->ind_array[i], pp_all->ind_array[j])) strength[i]++;
  
  /* Fitness values =  sum of strength values of dominators */
  for(int i = 0; i < size; i++) {
    int sum = 0;
    for(int j = 0; j < size; j++)
      if(con_dominates(pp_all->ind_array[j], pp_all->ind_array[i])) sum += strength[j];

    pp_all->ind_array[i]->fitness = sum;
    fitness_bucket[sum]++;
    fitness_bucket_mod[(sum / size)]++;
  }
  
  free(strength);
  return;
}


void calcDistances()
{
  int size = pp_all->size;
  
  /* initialize copies[] vector and NN[][] matrix */
  for(int i = 0; i < size; i++) {
    copies[i] = 1;
    for(int j = 0; j < size; j++) NN[i][j] = -1;
  }
    
  /* calculate distances */
  for(int i = 0; i < size; i++) {
    NN[i][0] = i;
    for(int j = i + 1; j < size; j++) {
      dist[i][j] = calcDistance(pp_all->ind_array[i], pp_all->ind_array[j]);
      assert(dist[i][j] < PISA_MAXDOUBLE);
      dist[j][i] = dist[i][j];
      if(dist[i][j] == 0) {
	NN[i][copies[i]] = j;
	NN[j][copies[j]] = i;
	copies[i]++;
	copies[j]++;
      }
    }
    dist[i][i] = 0;
  }
}


int getNN(int index, int k)
/* lazy evaluation of the k-th nearest neighbor
   pre-condition: (k-1)-th nearest neigbor is known already */
{
  assert(index >= 0);
  assert(k >= 0);
  assert(copies[index] > 0);
    
  if(NN[index][k] < 0) {
    double min_dist = PISA_MAXDOUBLE;
    int min_index = -1;
    int prev_min_index = NN[index][k-1];
    double prev_min_dist = dist[index][prev_min_index];
    assert(prev_min_dist >= 0);
    
    for(int i = 0; i < pp_all->size; i++) {
      double my_dist = dist[index][i];
	    
      if(my_dist < min_dist && index != i) {
	if(my_dist > prev_min_dist ||
	   (my_dist == prev_min_dist && i > prev_min_index)) {
	  min_dist = my_dist;
	  min_index = i;
	}
      }
    }
    
    NN[index][k] = min_index;
  }
  
  return(NN[index][k]);
}


double getNNd(int index, int k)
/* Returns the distance to the k-th nearest neigbor
   if this individual is still in the population.
   For for already deleted individuals, returns -1 */
{
  int neighbor_index = getNN(index, k);
    
  if(copies[neighbor_index] == 0) return(-1);
  else return(dist[index][neighbor_index]);
}


void environmentalSelection()
{
  int new_size = 0;

  static int nc = 0; /* call number for diagnostics purposes */
  ++nc;
    
  if(verbose) {
    char tmp [1024];
    sprintf(tmp, " spea2: #%d 1st front size (%d", nc, fitness_bucket[0]);
    if(fitness_bucket[0] > alpha) sprintf(tmp, "%s*)\n", tmp);
    else sprintf(tmp, "%s)\n", tmp);
    
    msg_to_file(diag_file, tmp);
  }

  if(fitness_bucket[0] > alpha) truncate_nondominated();
  else if(pp_all->size > alpha)	truncate_dominated();

  /* Move remaining individuals to top of array in 'pp_all' */
  for(int i = 0; i < pp_all->size; i++) {
    ind* temp_ind = pp_all->ind_array[i];
    if(temp_ind != NULL) {
      assert(copies[i] > 0);	    
      pp_all->ind_array[i] = NULL;
      pp_all->ind_array[new_size] = temp_ind;
      old_index[new_size] = i;
      new_size++;    
    }
  }
  assert(new_size <= alpha);
  pp_all->size = new_size;
    
  return;
}


void truncate_nondominated()
/* truncate from nondominated individuals (if too many) */
{
  /* delete all dominated individuals */
  for(int i = 0; i < pp_all->size; i++) {
    if(pp_all->ind_array[i]->fitness > 0) {
      free_ind(pp_all->ind_array[i]);
      pp_all->ind_array[i] = NULL;
      copies[i] = 0;
    }
  }
    
  /* truncate from non-dominated individuals */
  while(fitness_bucket[0] > alpha) {
    int* marked;
    int max_copies = 0;
    int count = 0;
    int delete_index;

    marked = (int*) chk_malloc(pp_all->size * sizeof(int));

    /* compute inds with maximal copies */
    for(int i = 0; i < pp_all->size; i++) {
      if(copies[i] > max_copies) {
	count = 0;
	max_copies = copies[i];
      }
      if(copies[i] == max_copies) {
	marked[count] = i;
	count++;
      }
    }
    
    //assert(count >= max_copies); 
    
    if(count > max_copies) {    
      int* neighbor;
      neighbor = (int*) chk_malloc(count * sizeof(int));
      for(int i = 0; i < count; i++)
	neighbor[i] = kth; /* k-th neighbor */
      
      while(count > max_copies) {
	double min_dist = PISA_MAXDOUBLE;
	int count2 = 0;
	
	for(int i = 0; i < count; i++) {
	  double my_dist = -1;
	  while(my_dist == -1 && neighbor[i] < pp_all->size) {
	    my_dist = getNNd(marked[i],neighbor[i]);
	    neighbor[i]++;
	  }
	  
	  if(my_dist < min_dist) {
	    count2 = 0;
	    min_dist = my_dist;
	  }
	  if(my_dist == min_dist) {
	    marked[count2] = marked[i];
	    neighbor[count2] = neighbor[i];
	    count2++;
	  }
	}
	
	count = count2;
	if(min_dist == -1) break; /* all have equal distances */
      }
      
      free(neighbor);
    }
    
    /* remove individual from population */
    delete_index = marked[irand(count)];
    free_ind(pp_all->ind_array[delete_index]);
    pp_all->ind_array[delete_index] = NULL;
    for(int i = 0; i < count; i++)
      if(dist[delete_index][marked[i]] == 0) copies[marked[i]]--;

    copies[delete_index] = 0; /* Indicates that this index is empty */
    fitness_bucket[0]--;
    fitness_bucket_mod[0]--;
    free(marked);
  }
  
  return;
}


void truncate_dominated()
/* truncate from dominated individuals */
{
  int size = pp_all->size;
  int num = 0;
  
  int i = -1;
  while (num < alpha) {
    i++;
    num += fitness_bucket_mod[i];
  }
  
  int j = i * size;
  num = num - fitness_bucket_mod[i] + fitness_bucket[j];
  while(num < alpha) {
    j++;
    num += fitness_bucket[j];
  }
    
  if(num == alpha) {
    for(i = 0; i < size; i++) {
      if(pp_all->ind_array[i]->fitness > j) {
	free_ind(pp_all->ind_array[i]);
	pp_all->ind_array[i] = NULL;
      }
    }
  }
  else { /* if not all fit into the next generation */
    int k;
    int free_spaces;
    int fill_level = 0;
    int* best;

    free_spaces = alpha - (num - fitness_bucket[j]);
    best = (int*) chk_malloc(free_spaces * sizeof(int));
    for (i = 0; i < size; i++) {
      if(pp_all->ind_array[i]->fitness > j) {
	free_ind(pp_all->ind_array[i]);
	pp_all->ind_array[i] = NULL;
      }
      else if(pp_all->ind_array[i]->fitness == j) {
	if(fill_level < free_spaces) {
	  best[fill_level] = i;
	  fill_level++;
	  for(k = fill_level - 1; k > 0; k--) {
	    int temp;
	    if(getNNd(best[k], kth) <= getNNd(best[k - 1], kth)) /* k-th neighbor */
	      break;
	    temp = best[k];
	    best[k] = best[k-1];
	    best[k-1] = temp;
	  }
	}
	else {
	  if(getNNd(i, kth) <= getNNd(best[free_spaces - 1], kth)) { /* k-th neighbor */
	    free_ind(pp_all->ind_array[i]);
	    pp_all->ind_array[i] = NULL;
	  }
	  else {
	    free_ind(pp_all->ind_array[best[free_spaces - 1]]);
	    pp_all->ind_array[best[free_spaces - 1]] = NULL;
	    best[free_spaces - 1] = i;
	    for(k = fill_level - 1; k > 0; k--) {
	      int temp;
	      if(getNNd(best[k], kth) <= getNNd(best[k - 1], kth)) /* k-th neighbor */
		break;
	      temp = best[k];
	      best[k] = best[k-1];
	      best[k-1] = temp;
	    }
	  }
	}
      }
    }
  }
  return;
}


void matingSelection()
/* Fills mating pool 'pp_sel' */
{
  for(int i = 0; i < mu; i++) {
    int winner = irand(pp_all->size);
	
    for(int j = 1; j < tournament; j++) {
      int opponent = irand(pp_all->size);
      if(pp_all->ind_array[opponent]->fitness
	 < pp_all->ind_array[winner]->fitness || winner == opponent)
	winner = opponent;
      else if(pp_all->ind_array[opponent]->fitness
	      == pp_all->ind_array[winner]->fitness) {
	if(dist[old_index[opponent]][getNN(old_index[opponent], kth)] /* k-th neighbor */
	   > dist[old_index[winner]][getNN(old_index[winner], kth)]) /* k-th neighbor */
	  winner = opponent;
      }
    }  
    pp_sel->ind_array[i] = pp_all->ind_array[winner];
  }
  pp_sel->size = mu;
  return;
}


void select_initial()
/* Performs initial selection. */
{
  selection();
}


void select_normal()
/* Performs normal selection.*/
{
  selection();
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


bool dominates(ind *p_ind_a, ind *p_ind_b)
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


int is_equal(ind *p_ind_a, ind *p_ind_b)
/* Determines if two individuals are equal in all objective values. */
{
  int equal = 1;
  
  for(int i = 0; i < dim; i++)
    equal = (p_ind_a->f[i] == p_ind_b->f[i]) && equal;
  
  return(equal);
}


double calcDistance(ind *p_ind_a, ind *p_ind_b)
{
  double distance = 0;
  
  if(is_equal(p_ind_a, p_ind_b) == 1) return(0);
  
  for (int i = 0; i < dim; i++)
    distance += pow(p_ind_a->f[i] - p_ind_b->f[i],2);
  
  return(sqrt(distance));
}


int irand(int range)
/* Generate a random integer. */
{
  int j = (int) ((double)range * (double) rand() / (RAND_MAX+1.0));
  return(j);
}


/*--------------------| data exchange functions |------------------------*/

int read_ini()
{
  pp_new = create_pop(alpha);
  
  for(int i = 0; i < alpha; i++)
    pp_new->ind_array[i] = create_ind();
  pp_new->size = alpha;
  
  return(read_pop(inifile, pp_new, alpha, dim, con));                    
}


int read_var()
{
  pp_new = create_pop(lambda);
  
  for(int i = 0; i < lambda; i++)
    pp_new->ind_array[i] = create_ind();
  
  pp_new->size = lambda;
  return(read_pop(varfile, pp_new, lambda, dim, con));
}


void write_sel() {
  write_pop(selfile, pp_sel, mu);
  return;
}


void write_arc()
{
  write_pop(arcfile, pp_all, pp_all->size);
  return;
}


int check_sel()
{
  return(check_file(selfile));
}


int check_arc()
{
  return(check_file(arcfile));
}
