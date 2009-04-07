#include "CESR_platform.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef CESR_VMS
#ifndef CESR_WINCVF
#include <readline/readline.h>
#include <readline/history.h>
#endif
#endif
#ifndef CESR_WINCVF
void read_line_(char* tag, char* str, int tag_len, int str_len) {
#ifndef CESR_VMS
  printf ("\n");
  char* str2 = readline (tag);
  strcpy (str, str2);
  int ii;
  for (ii = strlen(str2); ii < str_len; ii++) 
    str[ii] = ' ';

/* If the line has any text in it, save it on the history stack */

  if (str2 && *str2) add_history (str2);

  free (str2);
#endif
}
#endif
#ifdef CESR_WINCVF
 void READ_LINE(char* tag, char* str, int tag_len, int str_len) {
 }
#endif
