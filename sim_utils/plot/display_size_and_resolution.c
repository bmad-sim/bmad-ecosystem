//+
// void display_size_and_res(int ix_screen, double* x_size, double* y_size, double* x_res, double* y_res,
//                           int* x_pix, int* y_pix)
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
//    x_pix   -- Horizontal screen size in pixels. Set even if the physical size is not known.
//    y_pix   -- Vertical screen size in pixels. Set even if the physical size is not known.
//-

// No X11 case --------------------------------------------------------------

#if defined (CESR_NOPLOT)

void display_size_and_res(int ix_screen, double* x_size, double* y_size, double* x_res, double* y_res,
                          int* x_pix, int* y_pix)
{
  *x_size = 0;
  *y_size = 0;
  *x_res = 0;
  *y_res = 0;
  *x_pix = 0;
  *y_pix = 0;
}

// GTK case (ACC_PLOT_DISPLAY_TYPE = "GTK"): No X11 needed ---------------------

#elif defined (CESR_GTK)

#include "qp_gtk_window.h"

void display_size_and_res(int ix_screen, double* x_size, double* y_size, double* x_res, double* y_res,
                          int* x_pix, int* y_pix)
{
  int scale;
  // GDK resolution is in logical pixels, the units of the GTK plot window.
  qpg_display_info(x_size, y_size, x_res, y_res, x_pix, y_pix, &scale);
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

void display_size_and_res(int ix_screen, double* x_size, double* y_size, double* x_res, double* y_res,
                          int* x_pix, int* y_pix)
{

  *x_size = 0;
  *y_size = 0;
  *x_res = 0;
  *y_res = 0;
  *x_pix = 0;
  *y_pix = 0;

  Display *dpy; 
  char *displayname = NULL;            /* server to contact */
  dpy = XOpenDisplay (displayname);
  if (!dpy) return;
  if (ix_screen >= 0 && ix_screen < ScreenCount(dpy)) {

    //  printf ("name of display:    %s\n", DisplayString (dpy));
    //  printf ("default screen number:    %d\n", DefaultScreen (dpy));
    //  printf ("number of screens:    %d\n", ScreenCount (dpy));

    *x_pix = XDisplayWidth(dpy, ix_screen);
    *y_pix = XDisplayHeight(dpy, ix_screen);
    *x_size = XDisplayWidthMM(dpy, ix_screen);
    *y_size = XDisplayHeightMM (dpy, ix_screen);
    if (*x_size > 0 && *y_size > 0) {
      *x_res = (double) *x_pix / *x_size;
      *y_res = (double) *y_pix / *y_size;
    }
  }

  XCloseDisplay(dpy);
}

#endif
