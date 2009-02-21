#include <cstdio>
#include <cstdlib>
#include <cmath>


extern "C" void check_result1_(double* power_t, int& N)
{
   char buffer1 [50] ;
   char buffer2 [50] ;
   char buffer3 [50] ;
   int i;
   FILE *out_file;
   out_file=fopen("wout","w");
   for (i=0;i<N;i++) fprintf(out_file,"%d   %20.16e\n", i, power_t[i]);
   fclose(out_file);
     
   sprintf(buffer1, "gnuplot gnu1");
   system(buffer1);
   puts ("Trying to generate the plot of the Power V.S. t");
   sprintf(buffer2, "gv wt.ps &");
   sprintf(buffer3, "gv logwt.ps &");
   i = system (buffer2);
   system(buffer3);
   if (i==-1) puts ("Error in generating the plot");
   else puts ("Plot successfully generated");
}

