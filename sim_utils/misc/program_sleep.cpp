//+
// extern "C" void program_sleep_(const double& delay) 
// #include "cesr_utils.h"
//
// This routine pauses the program for a time (in seconds) set by delay. 
//-

#include <time.h>

extern "C" void program_sleep_(const double& delay) {

  struct timespec ts;
  ts.tv_sec = int(delay);
  ts.tv_nsec = int(1e9 * (delay - ts.tv_sec));
  nanosleep (&ts, NULL);

}
