/*
* double WordToReal (char *str, int *err_flag) {
*
* Function to convert the first word in a string to an double. 
* The only character allowed just after the word is a blanks, tabs, 
* or line_feed. 
*
* See: StringToReal as an alternative.
*
* #include "[cesr.dcs.lib]dcslib.h"
*/


#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <errno.h>

double WordToReal (char *str, int *err_flag) {

  double retval;
  char *endptr;

  errno = 0;
  retval = strtod (str, &endptr);

  if (endptr[0] != 32 && endptr[0] != 0 &&
      endptr[0] != 10 && endptr[0] != 9) 
    *err_flag = 1;
  else 
    *err_flag = errno;

  return retval;

}
