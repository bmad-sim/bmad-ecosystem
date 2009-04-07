/*
* double StringToReal (char *str, int *err_flag) {
*
* Function to convert a string to an double. No trailing characters allowed
* except blanks, tabs, and line_feeds
*
* See: WordToReal as an alternative.
*
* #include "[cesr.dcs.lib]dcslib.h"
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

double StringToReal (char *str, int *err_flag) {

  int i;
  char *endptr;
  double retval;

  errno = 0;
  retval = strtod (str, &endptr);

  *err_flag = errno;

  if (errno != 0) return retval;

  for (i=0; i<40; i++) {
    if (endptr[i] == 0) return retval;
    if (endptr[i] != 32 && endptr[i] != 0 &&
        endptr[i] != 10 && endptr[i] != 9) {
      *err_flag = 1;
      return retval;
    }
  }

  return retval;

}
