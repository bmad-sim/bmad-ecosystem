#include <cstdio>
#include <cstdlib>
#include <cmath>


extern "C" void check_result2_(double* var, double* thresh, int& N)
{
   char buffer1 [50] ;
   char buffer2 [50] ;
   int i;
   FILE *out_file;
   out_file=fopen("threshdep","w");
   for (i=0;i<N;i++) fprintf(out_file,"%20.16e   %20.16e\n", var[i], thresh[i]);
   fclose(out_file);
  
   
  
   sprintf(buffer1, "gnuplot gnu2");
   system(buffer1);
   puts ("Trying to generate the plot of the threshold current");
   sprintf(buffer2, "gv thresh.ps &");
   i = system (buffer2);
   if (i==-1) puts ("Error in generating the plot");
   else puts ("Plot successfully generated");
}

