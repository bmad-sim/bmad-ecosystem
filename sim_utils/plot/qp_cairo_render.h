/*
 * qp_cairo_render.h
 *
 * Shared drawing engine for the quick_plot Cairo backend.
 *
 * quick_plot drawing calls are encoded as a byte stream of "ops" (a display list).
 * The same stream is rendered by:
 *   - the quick_plot process itself, onto a PDF/PS/SVG/PNG surface (hardcopy), and
 *   - the separate qp_cairo_viewer program, onto its on-screen window.
 * Having a single renderer guarantees that the screen and hardcopy look the same.
 *
 * Coordinates in the stream are inches from the lower left corner of the page,
 * matching the quick_plot "basic" routine conventions. Text and symbol heights
 * are in points. Line widths are in quick_plot units (1 unit = 1/200 inch, as in PGPLOT).
 *
 * Stream format: each op is
 *   uint8  opcode
 *   uint32 payload length in bytes
 *   payload (native-endian doubles / int32s)
 * Unknown opcodes are skipped using the length so the format can grow.
 */

#ifndef QP_CAIRO_RENDER_H
#define QP_CAIRO_RENDER_H

#include <stddef.h>
#include <stdint.h>
#include <cairo.h>

enum qpr_opcode {
  QPR_OP_PAGE         = 1,   /* page_w, page_h [inch]                          */
  QPR_OP_CLEAR        = 2,   /* Paint the page white.                          */
  QPR_OP_STATE        = 3,   /* Full qpr_state (used to resync after a clear). */
  QPR_OP_PUSH         = 4,
  QPR_OP_POP          = 5,
  QPR_OP_COLOR        = 6,   /* r, g, b in [0, 1]                              */
  QPR_OP_LINE_WIDTH   = 7,   /* width [qp units]                               */
  QPR_OP_LINE_PATTERN = 8,   /* int32 pattern                                  */
  QPR_OP_CLIP         = 9,   /* int32 on/off                                   */
  QPR_OP_GRAPH_POS    = 10,  /* x1, x2, y1, y2 [inch]                          */
  QPR_OP_CHAR_SIZE    = 11,  /* height [pt]                                    */
  QPR_OP_SYM_SIZE     = 12,  /* height [pt]                                    */
  QPR_OP_SYM_FILL     = 13,  /* int32 fill pattern                             */
  QPR_OP_TEXT_BG      = 14,  /* int32 opaque, r, g, b                          */
  QPR_OP_POLYLINE     = 15,  /* int32 n, x[n], y[n]                            */
  QPR_OP_RECT         = 16,  /* x1, x2, y1, y2, int32 fill pattern             */
  QPR_OP_SYMBOL       = 17,  /* x, y, int32 symbol type                        */
  QPR_OP_TEXT         = 18,  /* x, y, angle [deg], justify, int32 n, chars[n]  */
  QPR_OP_ARROW        = 19,  /* x1, y1, x2, y2, head_angle, head_barb, head_size, int32 filled */
  QPR_OP_TAG          = 20   /* int32 n, chars[n]: Label following ops (EG "r1.g2") for GUI hit-testing.
                                Ignored by the renderer. Reserved: not emitted yet. */
};

/* Transport framing between quick_plot and a viewer/GUI. Each message is:
 *   uint32 magic, uint32 version, uint32 type, uint32 reserved, uint64 nbytes, nbytes of op stream
 * RESET replaces the viewer's display list with the payload; APPEND adds to it.
 * Future GUI -> quick_plot events (mouse etc.) will use the same header on a return channel. */

#define QPR_MAGIC    0x4C445051u   /* "QPDL" */
#define QPR_VERSION  1u
enum { QPR_MSG_RESET = 1, QPR_MSG_APPEND = 2, QPR_MSG_CLOSE = 3 };

typedef struct {
  uint32_t magic, version, type, reserved;
  uint64_t nbytes;
} qpr_msg_header;

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
  int32_t sym_fill;
  int32_t text_bg_opaque;
  double text_bg_rgb[3];
} qpr_state;

#define QPR_MAX_DEPTH 64

typedef struct {
  cairo_t *cr;
  double scale;               /* device units per inch */
  double ox, oy;              /* device offset of the page's upper left corner */
  double page_w, page_h;      /* inch */
  double min_line_width;      /* device units. EG 1 pixel on screen so thin lines stay visible. */
  int draw;                   /* 0 => only track state (used for text measuring) */
  int depth;
  qpr_state st[QPR_MAX_DEPTH];
} qpr_renderer;

/* Initialize a renderer drawing onto cr. Scale and offset may be set later with qpr_set_view. */
void qpr_init (qpr_renderer *r, cairo_t *cr, double scale, double ox, double oy, int draw);
void qpr_set_view (qpr_renderer *r, double scale, double ox, double oy);

/* Execute all ops in buf. Returns 0 on success, -1 if the stream is truncated/corrupt. */
int qpr_exec (qpr_renderer *r, const uint8_t *buf, size_t len);

/* Length of text (with PGPLOT escape sequences) in inches at the current char size. */
double qpr_text_len (qpr_renderer *r, const char *text, int n);

/* Reads the PAGE op at the start of a stream (if any) to get the page size. Returns 0 if found. */
int qpr_peek_page (const uint8_t *buf, size_t len, double *page_w, double *page_h);

/* Default state applied at the start of a page. */
void qpr_default_state (qpr_state *s);

#endif
