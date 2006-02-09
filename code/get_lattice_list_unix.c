/* int get_lattice_list_unix_(char *lat_list, int *num_lats, char *directory,
		long int list_len, long int dir_len)

   Function to get the names of the lattices of the form:
       directory/BMAD_*.*

   Input:
     directory -- *char : directory to search for lattice files

   Output:
     lat_list -- *char : list of lattice names
     num_lats -- *int  : number of lattices found

   Note: This is a UNIX compatible version of get_lattice_list.  Its
         functionality is identical to the VMS version, but in order to
	 pass arguments to/from Fortran, its inner workings are quite
	 different.  This should only be an issue for someone who wishes to
	 modify it, not someone who wishes to use it in place of its VMS
	 counterpart.

	 When called from Fortran, the trailing underscore is omitted.
	 */
/*----------------------------------------------------------
 $Id$

 $Log$
 Revision 1.6  2006/02/09 21:31:06  cesrulib
 windows porting

 Revision 1.5  2004/01/13 20:03:46  cesrulib
 Eliminate cast problem for gcc compiler.

 Revision 1.4  2002/02/23 20:32:16  dcs
 Double/Single Real toggle added

 Revision 1.3  2002/01/11 17:39:37  palmer
 Avoid module clashes on VMS with the symbol match_wild

 Revision 1.2  2001/10/08 17:18:14  rwh24
 DCS changes to f90 files.
 Bug fixes to c file.

 Revision 1.1  2001/09/27 18:33:14  rwh24
 UNIX compatibility updates

*----------------------------------------------------------*/

#include "CESR_platform.h"

#if defined(CESR_UNIX) || defined(CESR_LINUX)
#include <dirent.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

void trim(char *str, int len);
void pad_to_length(char *str, int len);
int c_match_wild(const char *pattern, const char *str);

int get_lattice_list_unix_(char *lat_list, int *num_lats, char *directory,
		long int list_len, long int dir_len)
{
  DIR *dir_pointer;
  struct dirent *dp;
  char tmp[80];
  
  /* set filename mask here using '*' and '?' */
  const char *mask = "bmad_*.*";

  trim(directory, dir_len);
  dir_pointer = opendir(directory);
  if (!dir_pointer) {
    printf("Error opening directory\n");
    return -1;
  }

  *num_lats=0;
  *lat_list='\0';
  for (dp = readdir(dir_pointer); dp != NULL; dp = readdir(dir_pointer)) {
    strcpy(tmp, dp->d_name);
    if (c_match_wild(mask, tmp)) {
      pad_to_length(tmp, list_len);
      strcat(lat_list, tmp);
      (*num_lats)++;
    }
  }

  return closedir(dir_pointer);
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


/* Pad strings with spaces */
void pad_to_length(char *str, int len)
{
  int i;
  while (strlen(str)<len) strcat(str, " ");
}


/* Check for wildcard match (could be improved) */
int c_match_wild(const char *pattern, const char *str)
{
  char c;
  const char *s;

  for(;;) {
    switch(c=*pattern++) {
    case 0:
      if (!*str)
        return 1;
      return 0;
    case '?':
      ++str;
      break;
    case '*':
      if (!*pattern)
        return 1;
      s=str;
      while (*s) {
        if (*s == *pattern && c_match_wild(pattern, s))
          return 1;
        ++s;
      }
      break;
    default:
      if (*str++ != c)
        return 0;
      break;
    }
  }
  /* Should never get here */
  return -1;
}

#endif
