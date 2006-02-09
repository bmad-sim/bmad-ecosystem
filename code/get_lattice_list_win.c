/* int get_lattice_list_win_(char *lat_list, int *num_lats, char *directory,
		long int list_len, long int dir_len)

  Windows version of get_lattice_list_unix_ .
  Function to get the names of the lattices of the form:
       directory\BMAD_*.*

   Input:
     directory -- *char : directory to search for lattice files

   Output:
     lat_list -- *char : list of lattice names
     num_lats -- *int  : number of lattices found


	 When called from Fortran, the trailing underscore is omitted. */

#include "CESR_platform.h"

#if defined(CESR_WINCVF)

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

void trim(char *str, int len);
int get_lattice_list_win_(char *lat_list, int *num_lats, char *directory,
			   long int list_len, long int dir_len);
/*  this main() for testing only

void main(){
  char directory[40]="C:\\user_local\\cesrulib", lat_list[500];
  long int list_len, dir_len;
  int num_lats, error;
  dir_len=23;
  list_len=500;
  error=get_lattice_list_win_(lat_list, &num_lats, directory,list_len,dir_len);
  if (error == 0){  printf(lat_list);}
  else { printf("error = %d\n",error);};
}
*/

int get_lattice_list_win_(char *lat_list, int *num_lats, char *directory,
		long int list_len, long int dir_len)
{
  FILE *pf;
  int err;
  char tmp[80],comline[160];

  /* creating DOS command line "dir /B directory\bmad_*.* > file.tmp" */  
  trim(directory, dir_len);
  *comline='\0';
  strcat(comline,"dir /B ");
  strcat(comline,directory);
  strcat(comline,"\\bmad_*.*");
  strcat(comline," > xxx0123xxx.tmp"); 
  //  printf("comline : %s \n",comline);

  err = system(comline); /* execute command line */
  if (err !=0) {return -1;};
 
  pf=fopen("xxx0123xxx.tmp","r");
  
  *num_lats=0;
  *lat_list='\0';
  while(!feof(pf)){
    (*num_lats)++;
    fscanf(pf,"%s\n",tmp);
    strcat(tmp," ");
    strcat(lat_list,tmp);
  }
  
  fclose(pf);
  
  err = system("erase xxx0123xxx.tmp");
  if (err != 0){
    printf("get_lattice_list_win err =%d\n",err);
    return -2;
  }
  return 0;
 }
  

/* Trim trailing spaces from FORTRAN string */
void trim(char *str, int len)
{
  int i;

  for (i=0; i<len; i++) {
    if (isspace(str[i])) {
      str[i]='\0';
      break;
    }
  }
}
#endif


