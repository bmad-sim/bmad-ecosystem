/*
 * qp_cairo.c
 *
 * quick_plot process side of the Cairo backend. Called from cairo_interface.f90.
 *
 * Every drawing call is encoded as an op (see qp_cairo_render.h) and:
 *   - Hardcopy pages (PS/PDF/SVG/PNG): the op is rendered immediately onto a Cairo surface.
 *   - Screen page ("X"): the op is appended to the page's display list, which is sent to a
 *     viewer process on flush. Only the ops added since the last flush are sent, except after
 *     a page clear when the whole list is resent.
 *
 * The viewer connection:
 *   - If QP_CAIRO_DISPLAY is set to "fd:N", ops are written to the already open file descriptor N.
 *     This is the hook for a GUI front end that launches Tao and draws the plots itself.
 *   - Otherwise the viewer program is spawned with popen. Its name is taken from QP_CAIRO_VIEWER,
 *     default is $ACC_EXE/qp_cairo_viewer if that exists, else "qp_cairo_viewer" on PATH.
 *     If the viewer window is closed the viewer is restarted on the next plot update.
 */

#if defined (CESR_CAIRO)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cairo.h>
#include <cairo-pdf.h>
#include <cairo-ps.h>
#include <cairo-svg.h>
#include "qp_cairo_render.h"

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#define popen  _popen
#define pclose _pclose
#else
#include <signal.h>
#endif

enum {QPC_SCREEN = 0, QPC_PS = 1, QPC_PNG = 2, QPC_SVG = 3, QPC_PDF = 4};

typedef struct {
  int kind;
  char file[1024];
  double w_in, h_in;
  cairo_surface_t *surf;
  cairo_t *cr;
  qpr_renderer r;
  /* Screen only */
  uint8_t *buf;
  size_t len, cap, sent;
  int need_reset;
  FILE *pipe;
  int is_popen;
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

/* Finish the op. local: also execute on the local renderer. record: append to the screen display list. */

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
/* Viewer transport. */

static int open_transport (qpc_dev *d) {
  const char *disp = getenv("QP_CAIRO_DISPLAY");
  const char *viewer = getenv("QP_CAIRO_VIEWER");
  char cmd[2048];

#ifndef _WIN32
  /* A closed viewer must not kill the program with SIGPIPE. Only change the default disposition. */
  struct sigaction sa;
  if (sigaction(SIGPIPE, NULL, &sa) == 0 && sa.sa_handler == SIG_DFL) {
    sa.sa_handler = SIG_IGN;
    sigaction(SIGPIPE, &sa, NULL);
  }
#endif

  if (disp && strncmp(disp, "fd:", 3) == 0) {
    d->pipe = fdopen(atoi(disp + 3), "wb");
    d->is_popen = 0;
  } else {
    char exe_path[1024];
    if (!viewer || !viewer[0]) {
      /* Prefer the viewer installed next to the Bmad programs. Else rely on PATH. */
      const char *acc_exe = getenv("ACC_EXE");
      FILE *f = NULL;
      viewer = "qp_cairo_viewer";
      if (acc_exe && acc_exe[0]) {
        snprintf(exe_path, sizeof exe_path, "%s/qp_cairo_viewer", acc_exe);
        if ((f = fopen(exe_path, "rb")) != NULL) {fclose(f); viewer = exe_path;}
      }
    }
    snprintf(cmd, sizeof cmd, "\"%s\" --width %d --height %d --title quick_plot",
             viewer, (int)(d->w_in * 72 + 0.5), (int)(d->h_in * 72 + 0.5));
    d->pipe = popen(cmd, "w");
    d->is_popen = 1;
  }

  if (!d->pipe) {
    fprintf(stderr, "qp_cairo: Cannot start plot viewer: %s\n", viewer ? viewer : disp);
    return -1;
  }
#ifdef _WIN32
  _setmode(_fileno(d->pipe), _O_BINARY);
#endif
  d->need_reset = 1;
  return 0;
}

static void close_transport (qpc_dev *d) {
  if (!d->pipe) return;
  if (d->is_popen) pclose(d->pipe);
  else fclose(d->pipe);
  d->pipe = NULL;
}

/* Set QP_CAIRO_DEBUG to trace viewer connection events on stderr. */
static void debug_msg (const char *what, qpc_dev *d) {
  if (getenv("QP_CAIRO_DEBUG"))
    fprintf(stderr, "qp_cairo: %s (list %lu bytes, sent %lu, reset %d)\n", what,
            (unsigned long)d->len, (unsigned long)d->sent, d->need_reset);
}

static int send_msg (qpc_dev *d, uint32_t type, const uint8_t *p, size_t n) {
  qpr_msg_header h = {QPR_MAGIC, QPR_VERSION, type, 0, (uint64_t)n};
  if (fwrite(&h, sizeof h, 1, d->pipe) != 1) return -1;
  if (n > 0 && fwrite(p, 1, n, d->pipe) != n) return -1;
  if (fflush(d->pipe) != 0) return -1;
  return 0;
}

static void flush_screen (qpc_dev *d) {
  int err = -1, attempt;
  if (d->kind != QPC_SCREEN) return;
  if (!d->need_reset && d->sent == d->len) return;

  /* If the viewer has gone away (EG the user closed its window) start a new one and resend
   * the whole page. Only one retry so a viewer that cannot start does not cause a loop. */
  for (attempt = 0; attempt < 2; attempt++) {
    if (!d->pipe) {
      if (open_transport(d) != 0) return;
      debug_msg("viewer started", d);
    }
    if (d->need_reset) err = send_msg(d, QPR_MSG_RESET, d->buf, d->len);
    else               err = send_msg(d, QPR_MSG_APPEND, d->buf + d->sent, d->len - d->sent);
    if (!err) break;
    debug_msg("viewer write failed", d);
    close_transport(d);
    d->need_reset = 1;
    if (!d->is_popen) return;   /* A GUI connection is not ours to restart. */
  }

  if (err) return;
  d->sent = d->len;
  d->need_reset = 0;
}

/*-------------------------------------------------------------------------*/
/* Page management. */

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
    cairo_ps_surface_set_eps(d->surf, 1);
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
  default:   /* Screen: a tiny surface only used for measuring text. */
    d->kind = QPC_SCREEN;
    d->surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 1, 1);
    break;
  }

  if (cairo_surface_status(d->surf) != CAIRO_STATUS_SUCCESS) {
    fprintf(stderr, "qp_cairo: Cannot create surface for %s: %s\n", d->file,
            cairo_status_to_string(cairo_surface_status(d->surf)));
    cairo_surface_destroy(d->surf);
    free(d);
    return -1;
  }

  d->cr = cairo_create(d->surf);
  qpr_init(&d->r, d->cr, scale, 0, 0, d->kind != QPC_SCREEN);
  if (kind == QPC_PNG) d->r.min_line_width = 1;
  devs[ndev++] = d;

  ob_begin(QPR_OP_PAGE); ob_d(w_in); ob_d(h_in); ob_end();
  ob_begin(QPR_OP_CLEAR); ob_end();
  d->need_reset = 1;
  return ndev - 1;
}

void qpc_close (void) {
  qpc_dev *d = cur_dev();
  if (!d) return;

  if (d->kind == QPC_SCREEN) {
    if (d->pipe) {
      send_msg(d, QPR_MSG_CLOSE, NULL, 0);
      close_transport(d);
    }
    free(d->buf);
  } else if (d->kind == QPC_PNG) {
    cairo_status_t st = cairo_surface_write_to_png(d->surf, d->file);
    if (st != CAIRO_STATUS_SUCCESS)
      fprintf(stderr, "qp_cairo: Error writing %s: %s\n", d->file, cairo_status_to_string(st));
  } else {
    cairo_surface_finish(d->surf);
  }

  if (cairo_status(d->cr) != CAIRO_STATUS_SUCCESS)
    fprintf(stderr, "qp_cairo: Drawing error: %s\n", cairo_status_to_string(cairo_status(d->cr)));

  cairo_destroy(d->cr);
  cairo_surface_destroy(d->surf);
  free(d);
  devs[--ndev] = NULL;
}

void qpc_end (void) {
  while (ndev > 0) qpc_close();
}

void qpc_flush (void) {
  qpc_dev *d = cur_dev();
  if (d) flush_screen(d);
}

/* Clear the page. For the screen the display list is restarted, then the current state
 * stack is replayed into it so ops that follow the clear draw the same as before. */

void qpc_clear (void) {
  qpc_dev *d = cur_dev();
  int i;
  if (!d) return;

  if (d->kind != QPC_SCREEN) {
    ob_begin(QPR_OP_CLEAR); ob_end();
    return;
  }

  d->len = 0;
  d->sent = 0;
  d->need_reset = 1;
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
void qpc_sym_fill (int f)        {ob_begin(QPR_OP_SYM_FILL); ob_i(f); ob_end();}

void qpc_graph_pos (double x1, double x2, double y1, double y2) {
  ob_begin(QPR_OP_GRAPH_POS); ob_d(x1); ob_d(x2); ob_d(y1); ob_d(y2); ob_end();
}

void qpc_text_bg (int opaque, double r, double g, double b) {
  ob_begin(QPR_OP_TEXT_BG); ob_i(opaque); ob_d(r); ob_d(g); ob_d(b); ob_end();
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
