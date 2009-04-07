#include "CESR_platform.h"

#ifndef CESR_WINCVF
#ifndef CESR_VMS
#include <stdio.h> 
#include <ctype.h> 
#include <errno.h> 
#include <termios.h>    
#include <sys/termios.h>   
#include <sys/fcntl.h>  
#include <sys/ioctl.h> 
#include <unistd.h>

static int oldf; 
static struct termios torig, tnocanon; 

void init_tty_char_() { 

   /* get the original file flags */ 
   oldf = fcntl(0, F_GETFL, 0); 

   /* and the original tty settings */ 
   tcgetattr(0, &torig); 
   tnocanon = torig; 

   /* set for non-canonical io */ 
   tnocanon.c_lflag &= ~ICANON; 
} 

/* -------------------------------------------------------- */ 
/* first arg is character returned, second arg is 0 for non-blocking, 
    non-zero for blocking io, third non-zero to flush typeahead */ 

void get_tty_char_c_(int* ii, int* block, int* flush) { 
   char ch; 
   int n, ret; 

   if (!*block) { 
      fcntl(0, F_SETFL, oldf | O_NONBLOCK); 
   } 
   tcsetattr(0, TCSANOW, &tnocanon); 

   /* 
#ifdef DEBUG 
   printf("type any character\n"); 
   fflush(stdin); 
   sleep(1); 
#endif 
*/ 

/* EWOULDBLOCK means no characters waiting to be read */ 

   ret = read(0, &ch, 1); 
   if (ret < 1) { 
      if (ret < 0 && errno != EWOULDBLOCK) perror("read"); 
      *ii = 0; 
   } else { 
      *ii = ch; 
   } 

/* flush typeahead if requested */ 

   if (*flush) { 
      tcflush(0, TCIFLUSH); 
   } 

/* restore original modes */ 
   tcsetattr(0, TCSANOW, &torig); 
   if (!*block) { 
      fcntl(0, F_SETFL, oldf); 
   } 
}
#else
void init_tty_char_() {
} 
#endif
#else
#ifdef CESR_WINCVF
void INIT_TTY_CHAR() {
} 
void GET_TTY_CHAR_C(int* ii, int* block, int* flush) { 
}
#endif
#endif























