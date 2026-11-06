/*
 * qp_gtk_prototype.c
 *
 * Standalone test of the GTK plot window with a readline command loop, standing in for Tao.
 * Checks that the window stays responsive while the program waits at the prompt (readline's
 * rl_event_hook pumps GTK events) and that drawing is sharp on HiDPI screens.
 *
 * Commands:
 *   draw [n]    Draw a test plot with n sine periods (default 3).
 *   busy <s>    Compute for s seconds without processing window events (like a Tao optimizer).
 *   info        Print display size and resolution from GDK.
 *   png <file>  Write the off-screen picture (device resolution) to a PNG file.
 *   close       Close the window (the next "draw" reopens it).
 *   quit
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <readline/readline.h>
#include <readline/history.h>
#include "qp_gtk_window.h"

static int periods = 3;

/* Stand-in for a PLplot page: frame, grid, a curve, symbols, and text with Greek and subscripts. */

static void render (int w, int h) {
  cairo_t *cr = qpg_begin_draw(w, h);
  double x0 = 60, y0 = 30, x1 = w - 30, y1 = h - 50;
  int i;

  if (!cr) {printf("(no display: nothing drawn)\n"); return;}

  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_paint(cr);

  /* Hairline grid: shows whether the screen resolution is used (sharp vs. blurred). */
  cairo_set_source_rgb(cr, 0.75, 0.75, 0.75);
  cairo_set_line_width(cr, 0.5);
  for (i = 1; i < 10; i++) {
    double x = x0 + i * (x1 - x0) / 10, y = y0 + i * (y1 - y0) / 10;
    cairo_move_to(cr, x, y0); cairo_line_to(cr, x, y1);
    cairo_move_to(cr, x0, y); cairo_line_to(cr, x1, y);
  }
  cairo_stroke(cr);

  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, 1);
  cairo_rectangle(cr, x0, y0, x1 - x0, y1 - y0);
  cairo_stroke(cr);

  /* Curve */
  cairo_set_source_rgb(cr, 0, 0, 0.8);
  cairo_set_line_width(cr, 1.5);
  for (i = 0; i <= 400; i++) {
    double t = (double)i / 400;
    double x = x0 + t * (x1 - x0), y = 0.5*(y0 + y1) - 0.4 * (y1 - y0) * sin(2 * M_PI * periods * t);
    if (i == 0) cairo_move_to(cr, x, y); else cairo_line_to(cr, x, y);
  }
  cairo_stroke(cr);

  /* Symbols */
  cairo_set_source_rgb(cr, 0.85, 0, 0);
  for (i = 0; i <= 20; i++) {
    double t = (double)i / 20;
    double x = x0 + t * (x1 - x0), y = 0.5*(y0 + y1) - 0.4 * (y1 - y0) * sin(2 * M_PI * periods * t);
    cairo_new_sub_path(cr);
    cairo_arc(cr, x, y, 3, 0, 2 * M_PI);
  }
  cairo_set_line_width(cr, 1);
  cairo_stroke(cr);

  /* Text */
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_select_font_face(cr, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
  cairo_set_font_size(cr, 14);
  cairo_move_to(cr, x0, y0 - 10);
  cairo_show_text(cr, "GTK + Cairo prototype:  β-function,  σ = 1.5 mm,  Å");
  cairo_set_font_size(cr, 11);
  cairo_move_to(cr, 0.5*(x0 + x1) - 30, h - 20);
  cairo_show_text(cr, "s [m]");
  cairo_move_to(cr, x0, h - 20);
  {
    char buf[80];
    snprintf(buf, sizeof buf, "%d x %d logical px", w, h);
    cairo_show_text(cr, buf);
  }

  qpg_end_draw(cr);
}

static void on_resize (int w, int h, void *data) {
  render(w, h);
}

/* readline calls this repeatedly while waiting for a key. */
static int event_hook (void) {
  qpg_pump();
  return 0;
}

int main (void) {
  char *line;
  int w = 600, h = 400;

  if (qpg_init() != 0) printf("No display available: window commands will do nothing.\n");

  qpg_set_resize_callback(on_resize, NULL);
  qpg_open(w, h, "qp_gtk_prototype");
  render(w, h);

  rl_event_hook = event_hook;
  rl_set_keyboard_input_timeout(20000);   /* microseconds: pump events about 50 times a second */

  while ((line = readline("proto> ")) != NULL) {
    char cmd[64] = "", arg[512] = "";
    if (*line) add_history(line);
    sscanf(line, "%63s %511s", cmd, arg);

    if (!strcmp(cmd, "draw")) {
      if (arg[0]) periods = atoi(arg);
      render(w, h);
    } else if (!strcmp(cmd, "busy")) {
      double secs = arg[0] ? atof(arg) : 5, x = 0;
      clock_t end = clock() + (clock_t)(secs * CLOCKS_PER_SEC);
      while (clock() < end) x += sin(x);   /* No qpg_pump() here, like a long Tao calculation. */
      printf("done (%g)\n", x);
    } else if (!strcmp(cmd, "info")) {
      double wmm, hmm, xr, yr;
      int sc;
      if (qpg_display_info(&wmm, &hmm, &xr, &yr, &sc) == 0)
        printf("monitor %.0f x %.0f mm, %.3f x %.3f logical px/mm, scale factor %d\n", wmm, hmm, xr, yr, sc);
      else
        printf("display info not available\n");
    } else if (!strcmp(cmd, "png")) {
      printf(qpg_write_png(arg[0] ? arg : "qp_gtk_prototype.png") == 0 ? "written\n" : "write failed\n");
    } else if (!strcmp(cmd, "close")) {
      qpg_close();
    } else if (!strcmp(cmd, "quit") || !strcmp(cmd, "exit")) {
      free(line);
      break;
    } else if (cmd[0]) {
      printf("commands: draw [n], busy <s>, info, png <file>, close, quit\n");
    }
    free(line);
  }

  qpg_close();
  return 0;
}
