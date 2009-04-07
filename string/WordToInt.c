/*
* int WordToInt (char *str, int *err_flag) {
*
* Function to convert the first word in a string to an integer. 
* The only character allowed just after the word is a blanks, tabs, 
* or line_feed. 
*
* See: StringToInt as an alternative.
*
* #include "[cesr.dcs.lib]dcslib.h"
*/


#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <errno.h>

int WordToInt (char *str, int *err_flag) {

  int retval;
  char *endptr;

  errno = 0;
  retval = strtol (str, &endptr, 0);

  if (endptr[0] != 32 && endptr[0] != 0 &&
      endptr[0] != 10 && endptr[0] != 9) 
    *err_flag = 1;
  else 
    *err_flag = errno;

  return retval;

}
