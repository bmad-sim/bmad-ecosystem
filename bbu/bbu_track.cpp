#include <cstdio>
#include <limits>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#define MAXIT 100
#define UNUSED (-1.11e30)
#define accuracy (1e-3)
#define LARGEN 1e300
#define SMALLN 1e-300 


const int dimension = 4;
const int turns = 2;
const double pi = 4*atan(1.0);
const double clight = 299792458.0;
const double mc2 = 511000;


using namespace std;



struct beam {
      double coor[dimension];
      struct beam *next;
};

struct cavity {
      int c_ind;
      struct beam *head;
      double time;
      double dt;
};

/* Transport phase space coordinates from one element to another */

void transport(double* mat,
	       int cav_ind,
	       const double *x1, 
	       double* x2)
{
   int i,j;
   double u;

   for (i=0;i<dimension;i++) x2[i]=0;  
     
   for (j=0;j<dimension;j++)
      for (i=0;i<dimension;i++)
	 x2[i] += mat[cav_ind*dimension*dimension + i*dimension + j] * x1[j];
}

void nrerror(char error_text[]){

/* Numerical Recipes standard error handler */

   printf("Numerical Recipes run-time error...\n");
   printf("%s\n",error_text);
   printf("...now exiting to system...\n");
   exit(1);
}

inline double SIGN(double x1, double x2){
   return (x1*x2<0 ? -x1 : x1);
}

inline double rnd_uniform(double ival, double rng) {
  return ival + (rand() / (RAND_MAX + 1.0) - 0.5) * rng;
}

inline int RTCI (double val) {
  double tmp = val - int(val);
  if (tmp > 0.5) return int(val + 1);
  else if (tmp < -0.5) return int(val - 1);
  else return int(val);
}

double theta2(double* v, int s, double& check) {
  
  double x = (s-1.0) * s /2.0;
  double base = s * (s-1.0) * s * (2.0 * s - 1.0)/6.0  - x * x;
  double xy = 0;
  double y = 0;
  double yy = 0;
  for (int i = 0; i < s; ++i) {
    xy += v[i+s] * i;
    y += v[i+s];
    yy += v[i+s]*v[i+s];
  }
  double ssxx = (s-1.0) * s * (2.0 * s - 1.0)/6.0 - x*x/s;
  double ssyy = yy - y*y/s;
  double ssxy = xy - x*y/s;
  double res = (s * xy - x * y) / base;
  check = (ssxy * ssxy)/(ssxx * ssyy);
  cout << "theta: " << res << "  correlation:  " << check <<endl;
  return res;
}



double tbbu(
      const double current,     // cur
	    int    hom_all,           // cav_num
      int    hom_total,         // hom_total
	    int*   hom,               // hom
	    int*   start,             // start
	    int    bunch_num,         // csbunch_num
	    double bfreq,             // csbfreq
	    int    N,                 // N
	    int    n_prt,             // n_prt
	    double d_amp,             // csd_amp
	    double pos_err_amp,       // cspos_err_amp
	    double t_noise,           // cst_noise
	    double*  hom_time,        // cavtime
	    double*  mat,             // cavmat
	    double*  homQ,            // homQ
	    double*  homRoverQ,       // homRoQ
	    double*  homf,            // homfreq
	    double* cosang,           // cosang
	    double* sinang,           // sinang
	    complex<double>* pwi,     // cspw
	    double* cpower,           // power
	    double* pcoor,            // pcoor
	    double* hompw,            // hompw
	    int& calls,               // cscalls
	    double& crl)              // crl
{
   
  int hom_uniq = hom_all;
  double* homa;
  double* homb;
  double qbunch;
  complex<double>* homb2;
  complex<double>* pw;
  int k,l,i,j,n,it,jt,init;
  double t,dt,dx,dy,power,TB,av,r;
  int outi = 0;
  int uniqidx = 0;
  int index;
  double coor[dimension];
  double progress, old_progress; 

  double*  dplx    = new double [hom_uniq];
  double*  dply    = new double [hom_uniq];
  
  old_progress = 0;
  
  homa = new double [hom_total];
  homb = new double [hom_total];
  homb2 = new complex<double> [hom_total * turns];
  pw   = new complex<double> [hom_total]; 

/* Calculate the bunch charge */

  if (current) qbunch = current / (bfreq *1000);
  
/* Calculate coefficients used in HOM wake function */
  
  for (i = 0; i < hom_total; i++){
	homa[i] = homRoverQ[i]*homf[i]*homf[i]*qbunch/(2.0*clight);
	homb[i] = -homf[i]/(2.0*homQ[i]);
        pw[i] = pwi[i];
     }
  
/* Use seq to track which cavity's HOM power needs to be updated by which bunch*/
/* Use bunch to represent the bunch train */

  struct beam* bunch;
  struct beam* p;
  struct cavity* seq; 
  seq = new struct cavity[hom_all*turns];
  bunch = new struct beam[bunch_num];

  TB = 1.0/bfreq;
  for (i=1; i< bunch_num;i++){
     bunch[i].next = &bunch[i-1];
     for (j=0;j<dimension;j++) bunch[i].coor[j] = 0; // Initializing the bunches' phase space coordinates
  }
  bunch[0].next = &bunch[bunch_num-1];
  for (j=0;j<dimension;j++) bunch[0].coor[j] = 0;


  k=0;
  for (i=1; i<hom_uniq*turns; i++){    
     while(k/bfreq <= hom_time[i] && k < bunch_num){
	seq[i].head = &bunch[k];
	k+=1;
     }
     k-=1;
     seq[i].c_ind = i;
     seq[i].time = hom_time[i]-k*TB;
  } 
  seq[0].c_ind = 0;
  seq[0].head = &bunch[0];
  seq[0].time = 0;
  
   
  for (i=0;i<hom_uniq;i++){
     j = i + hom_uniq;
     if (seq[i].time >= seq[j].time) seq[i].dt = seq[i].time - seq[j].time;
     else seq[i].dt = TB + seq[i].time - seq[j].time;
     seq[j].dt = TB - seq[i].dt;
  }

/*homb2[] is used to optimize the speed*/
  for (k=0; k < turns; k++)
     for (i=0; i < hom_uniq; i++){
	dt = seq[k*hom_uniq + i].dt;
	init = start[i];
     for(j=0; j < hom[i]; j++)
	homb2[k*hom_total + init + j] = polar(exp(homb[init+j]* dt),homf[init+j]* dt);
  }

/*Sort the array in the time order*/
  for (i=0; i<hom_uniq*turns-1; i++)
     for (j=hom_uniq*turns-1; j>i; j--){
	if (seq[j].time < seq[j-1].time){
	   k = seq[j-1].c_ind;
	   p = seq[j-1].head;
	   t = seq[j-1].time;
	   dt = seq[j-1].dt;
	   seq[j-1].c_ind = seq[j].c_ind;
	   seq[j-1].head = seq[j].head;
	   seq[j-1].time = seq[j].time;
	   seq[j-1].dt = seq[j].dt;
	   seq[j].c_ind = k;
	   seq[j].head = p;
	   seq[j].time = t;
	   seq[j].dt = dt;
       }
    }

  
  for (int i = 0; i < hom_uniq; ++i){
    dplx[i] = rnd_uniform(0, pos_err_amp);
    dply[i] = rnd_uniform(0, pos_err_amp);
  }

   // generate homs stuff
 

  ifstream srandseed("rndsd");
  if (! srandseed) srand(1);
  else {
    srandseed >> i;
    srand(i);
  }
  
  

  
  av = 1.0/(hom_uniq*turns);
  r = 1.0/N;
/*Tracking through the cavity*/
  for (n=0; n <= N ; n++){
     power = 0.0;
     for (i=0; i<hom_uniq*turns; i++){
        dt = seq[i].dt;
        uniqidx = seq[i].c_ind % hom_uniq;
	k = seq[i].c_ind / hom_uniq;
        init = start[uniqidx];
	dx = (*seq[i].head).coor[0]-dplx[uniqidx];
	dy = (*seq[i].head).coor[2]-dply[uniqidx];
        for (j = 0; j < hom[uniqidx]; ++j){
	   pw[init+j] = homb2[hom_total*k + init + j] * pw[init+j] +  homa[init+j] * (dx * cosang[init+j] + dy * sinang[init+j]);
	   power+= pw[init+j].real()*pw[init+j].real()+pw[init+j].imag()*pw[init+j].imag();
	   (*seq[i].head).coor[1] += pw[init+j].imag() * cosang[init+j];
	   (*seq[i].head).coor[3] += pw[init+j].imag() * sinang[init+j];
	}
        
       
        if (seq[i].c_ind != hom_uniq*turns-1){ 
	   transport(mat,seq[i].c_ind,(*seq[i].head).coor, coor);	   
	   for (l=0; l<dimension; l++) (*seq[i].head).coor[l] = coor[l];
	}
	else if(n%n_prt == 0){
           for (l=0; l<dimension; l++) pcoor[outi*dimension+l] = (*seq[i].head).coor[l];
	}

	  seq[i].head = (*seq[i].head).next;
     }

     for (l=0;l<dimension;l++) (*seq[0].head).coor[l] = 0;
   
     if (n % n_prt == 0) {
       cpower[outi] = power*av;
       outi++;
       progress = double(n)*100.0*r;
       if (progress > old_progress + 2) {
         printf("progress: % 5.1f%%\r", progress);
         fflush(stdout);
         old_progress = progress;
       }
     }
  }

     


      
    
  for (i = 0; i < outi; i++) {
    if (cpower[i] <= SMALLN) {
           printf("The trial current %g is in super stable region and will be increased\n", current);
	   return (-0.05);
	   break;
	}
	else if(cpower[i] >= LARGEN){ 
           printf("The trial current %g is in super unstable region and will be decreased\n", current);
	   return (0.05);
	   break;
	}
    cpower[i] = log(cpower[i]);

  }

  for (i = 0; i < hom_total; i++) {
     hompw[i] = sqrt(pw[i].real()*pw[i].real()+pw[i].imag()*pw[i].imag()) ;   
  }

 
  delete[] pw;
  delete[] homa;
  delete[] homb;
  delete[] homb2;
  delete[] dplx;
  delete[] dply;
  delete[] seq;
  delete[] bunch;

  calls++;
  cout << "call " << calls << ", trial current:" << current << "mA\n";
  return theta2(cpower, N/(2*n_prt), crl);
}
  


double root(double x1, double fx1, double x2, double fx2, double xacc,
	    int    hom_all,
            int    hom_total,
	    int*   hom,
	    int*   start,
	    int    bunch_num,
	    double bfreq,
	    int    N,
	    int    n_prt,
	    double d_amp,
	    double pos_err_amp,
	    double t_noise,
	    double*   hom_time,
	    double* mat,
	    double*  homQ,
	    double*  homRoverQ,
	    double*  homf,
	    double* cosang,
	    double* sinang,
	    complex<double>* pw,
	    double* cpower,
	    double* pcoor,
            double* hompw,
	    int& calls,
            double& crl)

{
   int j;
   double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
   fl = fx1;
   fh = fx2;

  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl=x1;
    xh=x2;                                          
    ans=UNUSED;
                                          
    for (j=1;j<=MAXIT;j++) {
        xm=0.5*(xl+xh);                                          
        fm=tbbu(xm, hom_all, hom_total, hom, start,bunch_num, bfreq, N, n_prt, d_amp, pos_err_amp, t_noise,
		   hom_time, mat, homQ, homRoverQ, homf, cosang, sinang, pw, cpower, pcoor, hompw,calls, crl);
        s=sqrt(fm*fm-fl*fh);
        if (s == 0.0) return ans;
        xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
        if (fabs((xh-xl)/xl) <= xacc)  return ans;	
        ans=xnew;
        fnew=tbbu(ans, hom_all, hom_total, hom,start, bunch_num, bfreq, N, n_prt, d_amp, pos_err_amp, t_noise,
		   hom_time, mat, homQ, homRoverQ, homf, cosang, sinang, pw, cpower, pcoor, hompw, calls, crl);
        if (fnew == 0.0) return ans;
        if (SIGN(fm,fnew) != fm) {
            xl=xm;
            fl=fm;
            xh=ans;
            fh=fnew;
        } else if (SIGN(fl,fnew) != fl) {
            xh=ans;
            fh=fnew;
        } else if (SIGN(fh,fnew) != fh) {
            xl=ans;
            fl=fnew;
        } else nrerror("never get here.");
        if (fabs((xh-xl)/xl) <= xacc){	 
           return ans;
	}
     }
    nrerror("zriddr exceeds maximum iterations");
}
   else {
      if (fabs(fl) == 0.0) return x1;
      if (fabs(fh) == 0.0) return x2;
      nrerror("root must be bracketed in zriddr.");
   }
  return 0.0;

}

 
extern "C" void bbu_track_(
double& cur,             // cur
int& cav_num,            // bbu%num_unique_cav
int& hom_total,          // bbu%hom_num
double* cavmat,          // bbu%mat_c
double* cavtime,         // bbu%time_c
int* hom,                // bbu%hom
double* homQ,            // bbu%Q_c
double* homRoQ,          // bbu%RoQ_c
double* cshomfreq,       // bbu%freq_c
double* homang,          // bbu%angle_c
double* hompw,           // bbu%power_c
double* power,           // power_t
double* pcoor,           // coor_t_c
double& csbfreq,         // bbu%b_freq
int& N,                  // bbu%b_time
int& n_prt,              // bbu%p_time
double& csd_amp,         // bbu%n_amp
double& cspos_err_amp,   // bbu%p_amp
double& cst_noise,       // bbu%n_on
int& check)              // 1

{ 
  double cur_low = 0.0;
  double cur_high = 0.0;
  double theta = 0.0;
  double crl = 0;
  double xcrl = 0;
  double theta_high = 0.0;
  double theta_low = 0.0;
  double alpha = 0.0;
  int cscalls = 0;
  int i,k;
  string bbuin;
  int*    start;
  double* homfreq; 
  double* cosang;
  double* sinang;
  complex<double>* cspw;


  int n_in = int(cavtime[cav_num*turns-1] *csbfreq); // Number of bunches in the lat at any one time.
 
 
  if (n_prt < 1) n_prt = 1;
  if (n_prt >= N) {
     cout << "printout time is too long\n";
     cout << "program stops";
     exit(1);
  }

		      
  if (N < n_in) {
     cout << "The number of bunches" << N << "<" << n_prt <<endl;
     cout << "beam does not fill the whole structure, increase run time\n";
     cout << "program stops";
     exit(1);
  }

  
  start   = new int    [cav_num];
  homfreq = new double [hom_total];
  cosang  = new double [hom_total];
  sinang  = new double [hom_total];
  cspw    = new complex<double> [hom_total];
 
  start[0] = 0;
  for (i=1; i < cav_num; i++) start[i] = start[i-1] + hom[i-1]; 
  for (i = 0; i < hom_total; i++) {
     homfreq[i] = 2*pi*cshomfreq[i];
     cosang[i] = cos(homang[i]);
     sinang[i] = sin(homang[i]);
     cspw[i] = hompw[i] * polar(1.0,pi*0.5);   
    }
 
 

 int csbunch_num=int(cavtime[cav_num*turns-1]*csbfreq); 
  
  
 theta = tbbu(cur, cav_num, hom_total, hom, start, csbunch_num, csbfreq, N,
	      n_prt, csd_amp, cspos_err_amp, cst_noise, cavtime, cavmat, homQ, homRoQ,
	      homfreq, cosang, sinang, cspw, power, pcoor, hompw, cscalls, crl);

 if (check == 0) {
    for (i=0;i <= N/n_prt;i++) power[i]=exp(power[i]);
    delete[] cspw;
    delete[] cosang;
    delete[] sinang;
    delete[] homfreq;
    delete[] start;
    return;
 }
 else {

    if (theta < 0) {
       cur_high = cur;
       theta_high = theta;
       while  (theta_high <0) {
     
	  cur_low = cur_high;
	  theta_low = theta_high;
	  cur_high *= 2;
  
    
	  theta_high = tbbu(cur_high, cav_num, hom_total, hom, start, csbunch_num, csbfreq, N, 
			    n_prt, csd_amp, cspos_err_amp, cst_noise, cavtime, cavmat, homQ, homRoQ,
			    homfreq, cosang, sinang, cspw, power, pcoor, hompw, cscalls, crl);
       }
    }
    else {
       cur_low = cur;
       theta_low = theta;
       while  (theta_low >0) {
     
	  cur_high = cur_low;
	  theta_high = theta_low;    
	  cur_low *= 0.5;
    
	  theta_low = tbbu(cur_low, cav_num, hom_total, hom, start, csbunch_num, csbfreq, N,
			   n_prt, csd_amp, cspos_err_amp, cst_noise, cavtime, cavmat, homQ, homRoQ,
			   homfreq, cosang, sinang, cspw, power, pcoor, hompw, cscalls, crl);
       }
    }
    cur = root(cur_low, theta_low, cur_high, theta_high, accuracy, cav_num, hom_total, hom, start, csbunch_num, csbfreq, N,
	       n_prt, csd_amp, cspos_err_amp, cst_noise, cavtime, cavmat, homQ, homRoQ,
	       homfreq, cosang, sinang, cspw, power, pcoor, hompw, cscalls, crl);

    theta = theta2(power, N/(2*n_prt), xcrl);
    printf("The threshold current is %g mA\n",cur);
    printf("The accuracy is %f%\n",accuracy * 100.0);
    printf("The residual power is %g\n", exp(power[N/n_prt-1]));
    printf("The ratio of residual power to the threshold current is %g\n", exp(power[N/n_prt-1])/cur);
    printf("The slope theta =%g and the correlation factor r*r = %g\n", theta, xcrl);
    cout << "number of calls: " << cscalls << endl;
    for (i=0;i<=N/n_prt;i++) power[i]=exp(power[i]);
 }
 
 

 delete[] cspw;
 delete[] cosang;
 delete[] sinang;
 delete[] homfreq;
 delete[] start;
 
}
