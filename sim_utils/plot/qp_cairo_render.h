/*
 * qp_cairo_render.h
 *
 * Drawing engine for the quick_plot Cairo backend (ACC_PLOT_PACKAGE = "cairo").
 *
 * quick_plot drawing calls are encoded as a byte stream of "ops" (a display list). Hardcopy pages
 * (PDF, PS, SVG, PNG) execute each op as it arrives. The screen page keeps the list so the plot can
 * be redrawn at any size when the window is resized.
 *
 * The drawing reproduces what quick_plot drew with PLplot and its Cairo drivers: Pango text with
 * the same font sizes, super/subscript rules and placement, symbols drawn as the same Unicode
 * glyphs, and the same line, dash and fill rules.
 *
 * Units in the stream: positions are inches from the lower left corner of the page.
 * Text heights are quick_plot points. Symbol heights are quick_plot points times any size factor.
 * Line widths are quick_plot units (1 unit = 0.005 inch, as in PGPLOT).
 *
 * Stream format: each op is
 *   uint8  opcode
 *   uint32 payload length in bytes
 *   payload (native-endian doubles / int32s)
 */

#ifndef QP_CAIRO_RENDER_H
#define QP_CAIRO_RENDER_H

#include <stddef.h>
#include <stdint.h>
#include <cairo.h>

enum qpr_opcode {
  QPR_OP_PAGE         = 1,   /* page_w, page_h [inch]                                      */
  QPR_OP_CLEAR        = 2,   /* Paint the page white.                                      */
  QPR_OP_STATE        = 3,   /* Full qpr_state (used to resync after a clear).             */
  QPR_OP_PUSH         = 4,
  QPR_OP_POP          = 5,
  QPR_OP_COLOR        = 6,   /* r, g, b in [0, 1]                                          */
  QPR_OP_LINE_WIDTH   = 7,   /* width [qp units]                                           */
  QPR_OP_LINE_PATTERN = 8,   /* int32 pattern                                              */
  QPR_OP_CLIP         = 9,   /* int32 on/off                                               */
  QPR_OP_GRAPH_POS    = 10,  /* x1, x2, y1, y2 [inch]                                      */
  QPR_OP_CHAR_SIZE    = 11,  /* height [pt]                                                */
  QPR_OP_SYM_SIZE     = 12,  /* height [pt]                                                */
  QPR_OP_POLYLINE     = 15,  /* int32 n, x[n], y[n]                                        */
  QPR_OP_RECT         = 16,  /* x1, x2, y1, y2, int32 fill pattern                         */
  QPR_OP_SYMBOL       = 17,  /* x, y, int32 symbol type                                    */
  QPR_OP_TEXT         = 18,  /* x, y, angle [deg], justify, int32 n, chars[n]              */
  QPR_OP_ARROW        = 19   /* x1, y1, x2, y2, head_angle, head_barb, head_size, int32 filled */
};

/* Line patterns, fill patterns and symbol types use the quick_plot enum values. */

enum { QPR_SOLID = 1, QPR_DASHED = 2, QPR_DASH_DOT = 3, QPR_DOTTED = 4, QPR_DASH_DOT3 = 5 };
enum { QPR_SOLID_FILL = 1, QPR_NO_FILL = 2, QPR_HATCHED = 3, QPR_CROSS_HATCHED = 4 };

typedef struct {
  double rgb[3];
  double line_width;          /* qp units */
  int32_t line_pattern;
  int32_t clip;
  double graph[4];            /* x1, x2, y1, y2 [inch] */
  double char_size;           /* pt */
  double sym_size;            /* pt */
} qpr_state;

#define QPR_MAX_DEPTH 64

typedef struct {
  cairo_t *cr;
  double scale;               /* device units per inch */
  double ox, oy;              /* device offset of the page's upper left corner */
  double page_w, page_h;      /* inch */
  int draw;                   /* 0 => only track state */
  int depth;
  qpr_state st[QPR_MAX_DEPTH];
} qpr_renderer;

/* Initialize a renderer drawing onto cr. cr may be NULL if draw = 0. */
void qpr_init (qpr_renderer *r, cairo_t *cr, double scale, double ox, double oy, int draw);
void qpr_set_view (qpr_renderer *r, double scale, double ox, double oy);

/* Execute all ops in buf. Returns 0 on success, -1 if the stream is truncated/corrupt. */
int qpr_exec (qpr_renderer *r, const uint8_t *buf, size_t len);

/* Approximate length of text in inches at the current char size (same estimate as the PLplot backend). */
double qpr_text_len (qpr_renderer *r, const char *text, int n);

#endif
