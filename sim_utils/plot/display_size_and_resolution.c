//+
// void display_size_and_res(int ix_screen, double* x_size, double* y_size, double* x_res, double* y_res)
//
// Routine to return the size and resolution of a display screen.
// Note: All output numbers will be set to zero if there is a problem obtaining the display information.
// 
// Input:
//    ix_screen -- Screen index. Generally this should be 0.
//
// Output:
//    x_size  -- Horizontal screen size in mm.
//    y_size  -- Vertical screen size in mm.
//    x_res   -- Horizontal resolution in pixels / mm.
//    y_res   -- Vertical resolution in pixels / mm.
//-

// No X11 case --------------------------------------------------------------

#if defined (CESR_NOPLOT)

void display_size_and_res(int ix_screen, double* x_size, double* y_size, double* x_res, double* y_res)
{
  *x_size = 0;
  *y_size = 0;
  *x_res = 0;
  *y_res = 0;
}

// Have X11 case --------------------------------------------------------------

#else

#if defined(CESR_WINCVF)
#include <X11/Xwindows.h>
#endif

// Fix for OSX X11 issue with conda-forge build
#include <X11/Xfuncproto.h>
#ifndef _X_SENTINEL
# define _X_SENTINEL(x)
#endif
#ifndef _X_DEPRECATED
# define _X_DEPRECATED
#endif

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>

void display_size_and_res(int ix_screen, double* x_size, double* y_size, double* x_res, double* y_res)
{

  *x_size = 0;
  *y_size = 0;
  *x_res = 0;
  *y_res = 0;

  Display *dpy; 
  char *displayname = NULL;            /* server to contact */
  dpy = XOpenDisplay (displayname);
  if (!dpy) return;
  if (ix_screen >= 0 && ix_screen < ScreenCount(dpy)) {

    //  printf ("name of display:    %s\n", DisplayString (dpy));
    //  printf ("default screen number:    %d\n", DefaultScreen (dpy));
    //  printf ("number of screens:    %d\n", ScreenCount (dpy));

    *x_size = XDisplayWidthMM(dpy, ix_screen);
    *y_size = XDisplayHeightMM (dpy, ix_screen);
    *x_res = (double) XDisplayWidth(dpy, ix_screen) / XDisplayWidthMM(dpy, ix_screen);
    *y_res = (double) XDisplayHeight(dpy, ix_screen) / XDisplayHeightMM(dpy, ix_screen);
  }

  XCloseDisplay(dpy);
}

#endif
