/* Called about every 20 ms while get_tty_char_c waits for a key. Set by the GTK plot window
   (qp_gtk_window.c) so the window keeps repainting while Tao waits in single mode. */

void (*acc_tty_idle_hook) (void) = 0;

#ifndef CESR_WINCVF

#include <stdio.h> 
#include <ctype.h> 
#include <errno.h> 
#include <termios.h>    
#include <sys/termios.h>   
#include <sys/fcntl.h>  
#include <sys/ioctl.h> 
#include <unistd.h>
#include <poll.h>

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

/* While blocking, call the idle hook (if set) every 20 ms. */

   if (*block && acc_tty_idle_hook) {
      struct pollfd pfd = {0, POLLIN, 0};
      while (acc_tty_idle_hook) {
         int np = poll(&pfd, 1, 20);
         if (np > 0 || (np < 0 && errno != EINTR)) break;
         acc_tty_idle_hook();
      }
   }

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

#else                  // CESR_WINCVF
#include <conio.h>

void init_tty_char_() {}

void get_tty_char_c_(int* ii, int* block, int* flush) {
   if (!*block && !kbhit()) {
       *ii = 0;
       return;
   }
   if (flush) {
       *ii = getche();
   }
   else {
       *ii = getch();
   }
}
#endif























