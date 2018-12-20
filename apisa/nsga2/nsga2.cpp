/*========================================================================
  APISA  (www.tik.ee.ethz.ch/pisa/; www.lepp.cornell.edu/~ib38/apisa/)
  ========================================================================
  Computer Engineering (TIK)
  ETH Zurich

  Cornell University
  Ithaca, NY 14853
  ========================================================================
  NSGA2

  Implementation in C++ for the selector side.
  
  Implements Petri net.
  
  file: nsga2.cpp
  author: Marco Laumanns, laumanns@tik.ee.ethz.ch

  revision by:  Stefan Bleuler, bleuler@tik.ee.ethz.ch
  rewritten by: Ivan Bazarov, bazarov@cornell.edu
  last change:  July 22, 2004
  ========================================================================
*/

/* CAUTION: <unistd.h> is not standard C
   It is used only for sleep() and usleep() in wait().
   In Windows use Sleep() in <windows.h> or implement busy waiting.
*/

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <ctime>
#include <cstring>

#include "nsga2.hpp"

#ifdef PISA_UNIX
#include <unistd.h>
#endif

#ifdef PISA_WIN
#include <windows.h>
#endif

/*------------------------------| main() |-------------------------------*/

int main(int argc, char* argv[]) {
  /* command line parameters */
  char paramfile[FILE_NAME_LENGTH];     /* file with local parameters */
  char filenamebase[FILE_NAME_LENGTH];  /* filename base,
                                              e.g., "dir/test." */
  double poll = 1.0;                    /* polling interval in seconds */
  
  /* other variables */
  int state = -1;
  char statefile[FILE_NAME_LENGTH];
  int result;
     
  /* reading command line parameters */
  if (argc != 4)
    PISA_ERROR("Selector: wrong number of arguments");
  sscanf(argv[1], "%s", paramfile);
  sscanf(argv[2], "%s", filenamebase);
  sscanf(argv[3], "%lf", &poll);
  
  /* generate name of statefile */
  sprintf(statefile, "%ssta", filenamebase);
  
  /* main loop */
  while (state != 6) { /* stop state for selector */
                       /* Caution: if reading of the statefile fails
                          (e.g. no permission) this is an infinite loop */
    state = read_flag(statefile);
    
    if(state == 1) { /* inital selection */
      initialize(paramfile, filenamebase);
      result = read_ini();   /* read ini file */
      if (result == 0) {     /* reading ini file successful */
	select_initial();    /* do selection */
	write_arc();         /* write arc file (all individuals
				that could ever be used again) */
	write_sel();         /* write sel file */
	state = 2;
	write_flag(statefile, state);
      }                      /* else don't do anything and wait again */
    }
    
    
    else if(state == 3) {    /* selection */
      if(check_arc() == 0 && check_sel() == 0) {
	result = read_var();  /* read var file */
	if (result == 0) {    /*reading var file successful */
	  select_normal(); /* do selection */
	  write_arc();     /* write arc file (all individuals
			      that could ever be used again) */
	  write_sel();     /* write sel file */
	  state = 2;
	  write_flag(statefile, state);
	} /* else don't do anything and wait again */
      } /* else don't do anything and wait again */
    }
          
    else if(state == 5) { /* variator just terminated,
			     here you can do what you want */
      state = 6;          /* e.g., terminate too */
      write_flag(statefile, state);
    }

    else if (state == 9) { /* variator ready for reset,
			      here you can do what you want */
      state = 10;          /* e.g., get ready for reset too */
      write_flag(statefile, state);
    }
          
    else if (state == 10) { /* reset */
      free_memory();
      state = 11;
      write_flag(statefile, state);
    }
          
    else wait(poll); /* state == -1 (reading failed) or state concerns variator */
  } /* state == 6 (stop) */
  
  free_memory();
  state = 7;
  write_flag(statefile, state);
  return (0);
}

/*--------------------| functions for control flow |---------------------*/

void write_flag(char* filename, int flag)
/* Write the state flag to given file. */
{
  FILE *fp;
  
  assert(0 <= flag && flag <= 11);
  
  fp = fopen(filename, "w");
  assert(fp != NULL);
  fprintf(fp, "%d", flag);
  fclose(fp);
}


int read_flag(char* filename)
/* Read state flag from given file. */
{
  int result;
  int flag = -1;
  FILE *fp;
  fp = fopen(filename, "r");
  if(fp != NULL) {
    result = fscanf(fp, "%d", &flag);
    fclose(fp);
    if (result == 1) { /* excatly one element read */
      if (flag < 0 || flag > 11)
	PISA_ERROR("Selector: Invalid state read from file.");
    }
  }
  return (flag);
}


void wait(double sec)
/* Makes the calling process sleep. */
/* pre: sec >= 0.01 */
{
#ifdef PISA_UNIX
  unsigned int int_sec;
  unsigned int usec;
  
  assert(sec > 0);
  
  int_sec = (unsigned int) floor(sec);
  usec = (unsigned int) floor((sec - floor(sec)) * 1e6);
  /* split it up, usleep can fail if argument greater than 1e6 */

     
  /* two asserts to make sure your file server doesn't break down */
  assert(!((int_sec == 0) && (usec == 0))); /* should never be 0 */
  assert((int_sec * 1e6) + usec >= 10000);  /* you might change this one
					       if you know what you are
					       doing */
  
  sleep(int_sec);
  usleep(usec);
#endif

#ifdef PISA_WIN
  unsigned int msec;
  assert(sec > 0);
  msec = (unsigned int) floor(sec * 1e3);
  assert(msec >= 10); /* making sure we are really sleeping for some time*/
  Sleep(msec);
#endif
  
}


int msg_to_file(const char* file, const char* str)
/* can be used for debugging purpose. Writes date and time to a 'file'
   on the first call, then appends string 'str' to the same file. */
{
  FILE* fp;
  static bool first_call = true;
  
  if(file != NULL) {

    struct tm* curr_date;
    char date_string[64];
    time_t now = time(0);
    curr_date = (struct tm*)localtime(&now);
    if(curr_date != NULL)
      strftime(date_string, sizeof (date_string), "%d.%m.%Y %H:%M:%S",
	       curr_date);
    else strcpy(date_string, "no date"); /* no data available */
      
    if(first_call) {
      first_call = false;
      fp = fopen(file, "w");
      assert(fp != NULL);
    }
    else {
      fp = fopen(file, "a");
      assert(fp != NULL);
    }
    
    fprintf(fp, "(%s) %s", date_string, str);
    fclose(fp);
    
  }
  else error("invalid file name");

  return(0);
}
