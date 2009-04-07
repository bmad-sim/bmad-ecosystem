/*
 * This is called internally by re_match.  It should probably not be called
 * by the user.
 */

#include "CESR_platform.h"

#if defined(CESR_UNIX) || defined(CESR_WINCVF)
#include <regex.h>
#include <stdlib.h>

/*
 * Match string against the extended regular expression in
 * pattern, treating errors as no match.
 *
 * Return 1 for match, 0 for no match.
 */


int match_reg_c_(const char *string, char *pattern)
{
  int status;
  regex_t re;
  
  if (regcomp(&re, pattern, REG_EXTENDED|REG_NOSUB) != 0) {
    return(0);      /* Report error. */
  }
  status = regexec(&re, string, (size_t) 0, NULL, 0);
  regfree(&re);
  if (status != 0) {
    return(0);      /* Report error. */
  }
  return(1);
}
#endif

#if defined(CESR_VMS) 
int match_reg_c(const char *string, char *pattern)
{
   printf("match_reg_c is not available under OpenVMS");
   return(0);
}
#endif
