/*

Function opti (vec, gen, pop, n_vars, merit, v0, v_del) 

Subroutines to minimize a merit function.

Code adapted from Dr. Dobbs Journal, April 1997 with modifications
Article: Differential Evolution by Kenneth Price and Rainer Storn
Source Code: page 78
Added end_flag so merit function can terminate search.    -- DCS

Note: To call opti from C use:
  opti_(vec, gen, pop, n_vars, merit, v0, v_del) 

Input:
  gen           -- Integer: Number of "generations" to let solutions evolve.
  pop           -- Integer: Number of "strands" (solutions) to use.
  n_vars        -- Integer: Number of "genes" (variables).
  v0(n_vars)    -- Real(8): Initial strand guess.
  v_del(n_vars) -- Real(8): Typical variation of the genes.

Output:
  vec(n_vars)   -- Real(8): Optimum solution.
  opti          -- Real(8): Figure of merit of the optimum solution.

The Merit function prototype is:

function merit (vec, end_flag) result (this_merit)
  real(8) vec(*)     ! Input: trial strand
  integer end_flag   ! Output: Set = 1 to terminate optimization.
  real(8) this_merit ! Output: Merit value corresponting to vec.
  ...
end function

*/



//######################################################################
//########        START OPTI.H      ####################################
//######################################################################
// Use genetic algorithm to optimize functions
#ifndef _OPTI_H_
#define _OPTI_H_

#include <cmath>
#include <cstdlib>
#include <iostream>

// #elif defined(CESR_WINCVF)
// #include <math.h>
// #include <stdlib.h>
// #include <iostream>
// #endif

using namespace std;

typedef double (*EvalFunc) (double *, int&);

inline int RandInt(int limit) {return rand()%limit;}
inline double UniRand() {return double(rand())/RAND_MAX;}

class Strand;
class GenePool;
class Opti;
class Mini;
class Maxi;


class Strand{
public:
 double *gene;
protected:
 int length;
 double cost;
public:
 Strand() : gene(0), length(0) {}
 Strand(int len) : length(len) {gene = new double[length];}
 Strand(Strand &s0) : length(s0.length){
  gene = new double[length];
  this->operator=(s0);
 }
 Strand(double *v, int len) : length(len) {
  gene = new double[length];
  for(int i=0; i<len; i++) gene[i] = v[i];
 }
 ~Strand() {if(gene) delete [] gene;}
 int Length(){ return length; }
 void Length(int len){
  if(length>0) delete [] gene;
  length = len;
  gene = new double[length];
 }
 double& Cost() {return cost;}
 void Init(double delta){
  for(int i=0; i<length; i++) gene[i] = 2.0*delta*(UniRand()-0.5);
 }
 void Init(double offset, double delta){
  for(int i=0; i<length; i++) gene[i] = 2.0*delta*(UniRand()-0.5) + offset;
 }
 void Init(Strand& s, double delta){
  for(int i=0; i<length; i++) gene[i] = s.gene[i] + 2.0*delta*(UniRand()-0.5);
 }
 void Init(double v[], double delta){
  for(int i=0; i<length; i++) gene[i] = v[i] + 2.0*delta*(UniRand()-0.5);
 }
 void Init(double v0[], double vdel[]){
  for(int i=0; i<length; i++)
   gene[i] = v0[i] + 2.0*vdel[i]*(UniRand()-0.5);
 }
 Strand& operator=(Strand& s){
  cost = s.cost;
  for(int i=0; i<length; i++) gene[i] = s.gene[i];
  return *this;
 }
 double operator()(int index){return gene[index];}
 friend class GenePool;
 friend class Opti;
 friend class Mini;
 friend class Maxi;
};    


class GenePool{
public:
 Strand *strand;
protected:
 int population, length;
public:
 GenePool() : strand(0), population(0), length(0) {}
 GenePool(int pop) : population(pop), length(0){
  strand = new Strand[population];
 }
 GenePool(int pop, int len) : population(pop){
  strand = new Strand[population];
  Length(len);
 }
 virtual ~GenePool() {if (strand) delete [] strand;}
 int Length() {return length;}
 void Length(int len){
  length = len;
  if(strand)
   for(int i=0; i<population; i++) strand[i].Length(length);
 }
 int Population() {return population;}
 void Population(int pop){
  if (strand) delete [] strand;
  population = pop;
  strand = new Strand[population];
 }
 void Init(double delta){
  if(strand && length>0)
   for(int i=0; i<population; i++) strand[i].Init(delta);
 }
 void Init(double delta, double offset){
  if(strand && length>0)
   for(int i=0; i<population; i++) strand[i].Init(delta, offset);
 }
 void Init(Strand& s, double delta){
  if(strand && length>0)
   for(int i=0; i<population; i++) strand[i].Init(s, delta);
 }
 void Init(double v[], double delta){
  if(strand && length>0)
   for(int i=0; i<population; i++) strand[i].Init(v, delta);
 }
 void Init(double v0[], double vdel[]){
  strand[0].Init(v0, 0.0);
  if(strand && length>0){
   for(int i=1; i<population; i++) strand[i].Init(v0, vdel);
  }
 }
};


class Opti : public GenePool{
protected:
 int generations;
 Strand trial;
 EvalFunc f;
 virtual int Better(int) = 0;
 virtual int Better(int, int) = 0;
 virtual int Best(){
     int index, theMin;
     for(index=1,theMin=0; index<population; index++)
  	  theMin = Better(index, theMin) ? index : theMin;
     return theMin;
 }
public:
 Opti(int gen, int pop, int len, EvalFunc func)
   : GenePool(pop, len), generations(gen), trial(len) { f = func; }
 void Setf(EvalFunc func) {f = func;}
 virtual Strand& Optimize();
 virtual Strand& Optimize(int numgen){generations = numgen; return Optimize();}
};


class Mini : public Opti{
public:
 Mini(int gen, int pop, int len, EvalFunc func)
   : Opti(gen, pop, len, func) {}
 virtual int Better(int index){
   	return trial.cost < strand[index].cost;
 }
 virtual int Better(int ind1, int ind2){
   	return strand[ind1].cost < strand[ind2].cost;
 }
};


class Maxi : public Opti{
public:
 Maxi(int gen, int pop, int len, EvalFunc func)
   : Opti(gen, pop, len, func) {}
 virtual int Better(int index){
   	return trial.cost > strand[index].cost;
 }
 virtual int Better(int ind1, int ind2){
   	return strand[ind1].cost > strand[ind2].cost;
 }
};

double MiniVec(double solution[], const int& gen, const int& pop, 
	const int& len, EvalFunc func, double v0[], 
	const double& delta);
double MiniOff(double solution[], const int& gen, const int& pop,
	const int& len, EvalFunc func, const double& offset,
	const double& delta);
extern "C" double minidel(double solution[], const int& gen, const int& pop,
	const int& len, EvalFunc func, double v0[], double vdel[]);
double MaxiVec(double solution[], const int& gen, const int& pop, 
	const int& len, EvalFunc func, double v0[], 
	const double& delta);
double MaxiDel(double solution[], const int& gen, const int& pop,
	const int& len, EvalFunc func, double v0[], double vdel[]);
      
#endif

//######################################################################
//########          END OPTI.H      ####################################
//######################################################################

#if defined(CESR_WINCVF)
extern "C" void ERR_EXIT(); 
#else
extern "C" void err_exit_();
#endif

Strand& Opti::Optimize(){
 int a, b, c, j, end_flag;
 const double CR = 0.8, frac = 1.4;
 double rfrac, flex_cr;

 flex_cr = CR * 10.0/length;
 if ( 0.1 > flex_cr ) flex_cr = 0.1;
 if ( CR < flex_cr ) flex_cr = CR;

 for(int i=0; i<population; i++) strand[i].cost = f(strand[i].gene, end_flag);

 int numgen = generations;
 while(numgen--){
  for(int i=0; i<population; i++){
   do {a = RandInt(population);} while(a==i);
   do {b = RandInt(population);} while(b==i || b==a);
   do {c = RandInt(population);} while(c==i || c==b || c==a);
   j = RandInt(length);
   rfrac = frac * UniRand();
   for(int k=1; k<=length; k++){
     if(UniRand() < flex_cr || k == length){
      trial.gene[j] = strand[a].gene[j] +
	rfrac * (strand[b].gene[j] - strand[c].gene[j]);
     }
     else trial.gene[j] = strand[a].gene[j];
     j = (j+1)%length;
   }
   end_flag = 0;
   trial.cost = f(trial.gene, end_flag);
   if (Better(i)) strand[i] = trial;
   if (end_flag == 1) return strand[Best()];
  }
 }
 return strand[Best()];
}



double MiniOff(double solution[], const int& gen, const int& pop,
	const int& len, EvalFunc func, const double& offset,
	const double& delta){

 Mini xxx(gen, pop, len, func);
 xxx.Init(offset, delta);
 Strand sol = xxx.Optimize();
 for(int index=0; index<len; index++) solution[index] = sol(index);
 return sol.Cost();
}


double MiniVec(double solution[], const int& gen, const int& pop, 
	const int& len, EvalFunc func, double v0[], 
	const double& delta){

 Mini xxx(gen, pop, len, func);
 xxx.Init(v0, delta);
 Strand sol = xxx.Optimize();
 for(int index=0; index<len; index++) solution[index] = sol(index);
 return sol.Cost();
}

#if defined(CESR_WINCVF)
extern "C" double OPTI(double solution[], const int& gen, const int& pop, 
	const int& len, EvalFunc func, double v0[], 
	double vdel[]){
#else
extern "C" double opti_(double solution[], const int& gen, const int& pop, 
	const int& len, EvalFunc func, double v0[], 
	double vdel[]){

#endif

 if (pop < 4) {
   cout << "ERROR IN OPTI: POPULATION MUST BE AT LEAST 4!\n";

#if defined(CESR_WINCVF)
   ERR_EXIT(); 
#else
   err_exit_();
#endif
 }

 Mini xxx(gen, pop, len, func);
 xxx.Init(v0, vdel);
 Strand sol = xxx.Optimize();
 for(int index=0; index<len; index++) solution[index] = sol(index);
 return sol.Cost();
}


double MaxiVec(double solution[], const int& gen, const int& pop, 
	const int& len, EvalFunc func, double v0[], const double& delta){

 Maxi xxx(gen, pop, len, func);
 xxx.Init(v0, delta);
 Strand sol = xxx.Optimize();
 for(int index=0; index<len; index++) solution[index] = sol(index);
 return sol.Cost();
}


double MaxiDel(double solution[], const int& gen, const int& pop, 
	const int& len, EvalFunc func, double v0[], double vdel[]){

 Maxi xxx(gen, pop, len, func);
 xxx.Init(v0, vdel);
 Strand sol = xxx.Optimize();
 for(int index=0; index<len; index++) solution[index] = sol(index);
 return sol.Cost();
}


