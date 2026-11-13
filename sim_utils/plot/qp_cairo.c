/*
 * qp_cairo.c
 *
 * Page handling for the quick_plot Cairo backend (ACC_PLOT_PACKAGE = "cairo").
 * Called from cairo_interface.f90.
 *
 * Every drawing call is encoded as an op (see qp_cairo_render.h) and:
 *   - Hardcopy pages (PS, PDF, SVG, PNG): the op is drawn right away onto a Cairo surface.
 *   - Screen page: the op is added to the page's display list and drawn into the GTK plot window
 *     (qp_gtk_window.c). When the window is resized the whole list is redrawn to fill the window.
 */

#if defined (CESR_CAIRO)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <cairo.h>
#include <cairo-pdf.h>
#include <cairo-ps.h>
#include <cairo-svg.h>
#include "qp_cairo_render.h"
#include "qp_gtk_window.h"

enum {QPC_SCREEN = 0, QPC_PS = 1, QPC_PNG = 2, QPC_SVG = 3, QPC_PDF = 4};

typedef struct {
  int kind;
  char file[1024];
  double w_in, h_in;
  cairo_surface_t *surf;
  cairo_t *cr;
  qpr_renderer r;             /* Hardcopy: draws. Screen: state tracking only. */
  /* Screen only */
  uint8_t *buf;               /* Display list. */
  size_t len, cap;
  size_t drawn;               /* Part of the list drawn in the window. */
  int need_redraw;            /* Redraw the whole list (after a clear or resize). */
  int view_w, view_h;         /* Window drawing size [logical pixels]. */
  qpr_renderer sr;            /* Window renderer (keeps state between partial draws). */
} qpc_dev;

#define QPC_MAX_DEV 16
static qpc_dev *devs[QPC_MAX_DEV];
static int ndev = 0;

/* Scratch buffer used to build one op. */
static uint8_t *ob = NULL;
static size_t oblen = 0, obcap = 0;

static qpc_dev *cur_dev (void) {return ndev > 0 ? devs[ndev-1] : NULL;}

/*-------------------------------------------------------------------------*/

static void grow (uint8_t **b, size_t *cap, size_t need) {
  if (need <= *cap) return;
  size_t nc = *cap ? *cap : 4096;
  while (nc < need) nc *= 2;
  *b = (uint8_t *)realloc(*b, nc);
  if (!*b) {fprintf(stderr, "qp_cairo: out of memory\n"); exit(1);}
  *cap = nc;
}

static void ob_put (const void *p, size_t n) {
  grow(&ob, &obcap, oblen + n);
  memcpy(ob + oblen, p, n);
  oblen += n;
}

static void ob_begin (uint8_t op) {
  uint32_t zero = 0;
  oblen = 0;
  ob_put(&op, 1);
  ob_put(&zero, 4);
}

static void ob_d (double v) {ob_put(&v, sizeof v);}
static void ob_i (int32_t v) {ob_put(&v, sizeof v);}

/* Finish the op: execute it on the page's renderer and, for the screen, add it to the display list. */

static void ob_end_ex (int local, int record) {
  qpc_dev *d = cur_dev();
  uint32_t plen = (uint32_t)(oblen - 5);
  memcpy(ob + 1, &plen, 4);
  if (!d) return;
  if (local) qpr_exec(&d->r, ob, oblen);
  if (record && d->kind == QPC_SCREEN) {
    grow(&d->buf, &d->cap, d->len + oblen);
    memcpy(d->buf + d->len, ob, oblen);
    d->len += oblen;
  }
}

static void ob_end (void) {ob_end_ex(1, 1);}

/*-------------------------------------------------------------------------*/
/* Screen window. */

static double now (void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

/* Draw the not yet drawn part of the display list (or all of it) into the window. */

static void draw_screen (qpc_dev *d) {
  cairo_t *cr;
  if (d->kind != QPC_SCREEN) return;
  if (!d->need_redraw && d->drawn == d->len) return;

  cr = qpg_begin_draw(d->view_w, d->view_h);
  if (!cr) return;    /* No display. */

  if (d->need_redraw) {
    /* Stretch the page to fill the window: positions are scaled
     * separately in x and y, while text, symbols and line widths are scaled by the smaller of
     * the two so they keep their shape and never grow past what fits. */
    double sx = d->view_w / d->w_in, sy = d->view_h / d->h_in;
    double sc = sx < sy ? sx : sy;
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);
    qpr_init(&d->sr, cr, sc, 0, 0, 1);
    qpr_set_xy_scale(&d->sr, sx, sy);
    d->drawn = 0;
    d->need_redraw = 0;
  }

  d->sr.cr = cr;
  qpr_exec(&d->sr, d->buf + d->drawn, d->len - d->drawn);
  d->drawn = d->len;
  d->sr.cr = NULL;
  qpg_end_draw(cr);

  /* Testing aid: save what the window shows. */
  if (getenv("QP_GTK_SNAPSHOT")) qpg_write_png(getenv("QP_GTK_SNAPSHOT"));
}

/* Called by the window when its size changes. */

static void on_resize (int width, int height, void *data) {
  qpc_dev *d = (qpc_dev *)data;
  d->view_w = width;
  d->view_h = height;
  d->need_redraw = 1;
  draw_screen(d);
}

/*-------------------------------------------------------------------------*/
/* Page management. */

/* Open a page. Sizes are in inches. dpi: pixels per inch for PNG and for the screen window. */

int qpc_open (int kind, const char *file, int nfile, double w_in, double h_in, double dpi) {
  qpc_dev *d;
  double scale = 72;

  if (ndev >= QPC_MAX_DEV) return -1;
  d = (qpc_dev *)calloc(1, sizeof(qpc_dev));
  d->kind = kind;
  d->w_in = w_in;
  d->h_in = h_in;
  if (nfile > (int)sizeof(d->file) - 1) nfile = sizeof(d->file) - 1;
  memcpy(d->file, file, nfile);
  d->file[nfile] = 0;

  switch (kind) {
  case QPC_PS:
    d->surf = cairo_ps_surface_create(d->file, w_in * 72, h_in * 72);
    break;
  case QPC_PDF:
    d->surf = cairo_pdf_surface_create(d->file, w_in * 72, h_in * 72);
    break;
  case QPC_SVG:
    d->surf = cairo_svg_surface_create(d->file, w_in * 72, h_in * 72);
    break;
  case QPC_PNG:
    scale = dpi;
    d->surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, (int)(w_in * dpi + 0.5), (int)(h_in * dpi + 0.5));
    break;
  default:
    d->kind = QPC_SCREEN;
    d->view_w = (int)(w_in * dpi + 0.5);
    d->view_h = (int)(h_in * dpi + 0.5);
    if (qpg_screen_open(d->view_w, d->view_h, "quick_plot") != 0)
      fprintf(stderr, "quick_plot: No display available: plots will not be shown on screen.\n");
    qpg_set_resize_callback(on_resize, d);
    break;
  }

  if (d->surf && cairo_surface_status(d->surf) != CAIRO_STATUS_SUCCESS) {
    fprintf(stderr, "qp_cairo: Cannot create surface for %s: %s\n", d->file,
            cairo_status_to_string(cairo_surface_status(d->surf)));
    cairo_surface_destroy(d->surf);
    free(d);
    return -1;
  }

  if (d->surf) d->cr = cairo_create(d->surf);
  qpr_init(&d->r, d->cr, scale, 0, 0, d->kind != QPC_SCREEN);
  devs[ndev++] = d;

  ob_begin(QPR_OP_PAGE); ob_d(w_in); ob_d(h_in); ob_end();
  ob_begin(QPR_OP_CLEAR); ob_end();
  d->need_redraw = 1;
  return ndev - 1;
}

void qpc_close (void) {
  qpc_dev *d = cur_dev();
  if (!d) return;

  if (d->kind == QPC_SCREEN) {
    qpg_set_resize_callback(NULL, NULL);
    qpg_screen_close();
    free(d->buf);
  } else if (d->kind == QPC_PNG) {
    cairo_status_t st = cairo_surface_write_to_png(d->surf, d->file);
    if (st != CAIRO_STATUS_SUCCESS)
      fprintf(stderr, "qp_cairo: Error writing %s: %s\n", d->file, cairo_status_to_string(st));
  } else {
    cairo_surface_finish(d->surf);
  }

  if (d->cr) {
    if (cairo_status(d->cr) != CAIRO_STATUS_SUCCESS)
      fprintf(stderr, "qp_cairo: Drawing error: %s\n", cairo_status_to_string(cairo_status(d->cr)));
    cairo_destroy(d->cr);
  }
  if (d->surf) cairo_surface_destroy(d->surf);
  free(d);
  devs[--ndev] = NULL;
}

void qpc_end (void) {
  while (ndev > 0) qpc_close();
}

/* Show what has been drawn. Window events are processed at most every 20 ms here unless force
 * is nonzero: the input hooks (qp_gtk_window.c) keep the window up to date while the program waits for
 * input. quick_plot forces processing when it finishes a plot (qp_wait_to_flush(.false.)) so the
 * finished plot is always shown even when the program does not wait for input in between, EG
 * Tao single mode with a key held down. */

void qpc_flush (int force) {
  static double last_pump = 0;
  qpc_dev *d = cur_dev();
  double t;
  if (!d || d->kind != QPC_SCREEN) return;
  draw_screen(d);
  t = now();
  if (force || t - last_pump > 0.02) {
    qpg_pump();
    last_pump = t;
  }
}

/* Clear the page. For the screen the display list is restarted, then the current state
 * stack is replayed into it so ops that follow the clear draw the same as before. The window
 * keeps showing the old plot until the next flush, so it never shows a blank page while a new
 * plot is being drawn. */

void qpc_clear (void) {
  qpc_dev *d = cur_dev();
  int i;
  if (!d) return;

  if (d->kind != QPC_SCREEN) {
    ob_begin(QPR_OP_CLEAR); ob_end();
    return;
  }

  d->len = 0;
  d->need_redraw = 1;
  ob_begin(QPR_OP_PAGE); ob_d(d->w_in); ob_d(d->h_in); ob_end_ex(0, 1);
  ob_begin(QPR_OP_CLEAR); ob_end_ex(0, 1);
  for (i = 0; i <= d->r.depth; i++) {
    ob_begin(QPR_OP_STATE); ob_put(&d->r.st[i], sizeof(qpr_state)); ob_end_ex(0, 1);
    if (i < d->r.depth) {ob_begin(QPR_OP_PUSH); ob_end_ex(0, 1);}
  }
}

/*-------------------------------------------------------------------------*/
/* State and drawing. */

void qpc_push (void) {ob_begin(QPR_OP_PUSH); ob_end();}
void qpc_pop (void)  {ob_begin(QPR_OP_POP); ob_end();}

void qpc_color (double r, double g, double b) {
  ob_begin(QPR_OP_COLOR); ob_d(r); ob_d(g); ob_d(b); ob_end();
}

void qpc_line_width (double w)   {ob_begin(QPR_OP_LINE_WIDTH); ob_d(w); ob_end();}
void qpc_line_pattern (int p)    {ob_begin(QPR_OP_LINE_PATTERN); ob_i(p); ob_end();}
void qpc_clip (int on)           {ob_begin(QPR_OP_CLIP); ob_i(on); ob_end();}
void qpc_char_size (double h)    {ob_begin(QPR_OP_CHAR_SIZE); ob_d(h); ob_end();}
void qpc_sym_size (double h)     {ob_begin(QPR_OP_SYM_SIZE); ob_d(h); ob_end();}

void qpc_graph_pos (double x1, double x2, double y1, double y2) {
  ob_begin(QPR_OP_GRAPH_POS); ob_d(x1); ob_d(x2); ob_d(y1); ob_d(y2); ob_end();
}

void qpc_polyline (int n, const double *x, const double *y) {
  ob_begin(QPR_OP_POLYLINE); ob_i(n);
  ob_put(x, n * sizeof(double));
  ob_put(y, n * sizeof(double));
  ob_end();
}

void qpc_rect (double x1, double x2, double y1, double y2, int fill) {
  ob_begin(QPR_OP_RECT); ob_d(x1); ob_d(x2); ob_d(y1); ob_d(y2); ob_i(fill); ob_end();
}

void qpc_symbol (double x, double y, int type) {
  ob_begin(QPR_OP_SYMBOL); ob_d(x); ob_d(y); ob_i(type); ob_end();
}

void qpc_text (double x, double y, double angle, double justify, const char *text, int n) {
  ob_begin(QPR_OP_TEXT); ob_d(x); ob_d(y); ob_d(angle); ob_d(justify); ob_i(n);
  ob_put(text, n);
  ob_end();
}

void qpc_arrow (double x1, double y1, double x2, double y2, double head_angle,
                double head_barb, double head_size, int filled) {
  ob_begin(QPR_OP_ARROW);
  ob_d(x1); ob_d(y1); ob_d(x2); ob_d(y2); ob_d(head_angle); ob_d(head_barb); ob_d(head_size);
  ob_i(filled);
  ob_end();
}

double qpc_text_len (const char *text, int n) {
  qpc_dev *d = cur_dev();
  if (!d) return 0;
  return qpr_text_len(&d->r, text, n);
}

#endif
