/*
 * qp_cairo_render.c
 *
 * Drawing engine for the quick_plot Cairo backend. See qp_cairo_render.h.
 *
 * To look the same as the PLplot backend, the PLplot Cairo driver's rules are reproduced here:
 *   - Text is laid out by Pango from the same markup PLplot generates (font family, size in
 *     points = 0.25 * 72/25.4 * quick_plot height, super/subscript sizes and rises).
 *   - Plot symbols are the Unicode glyphs PLplot draws for each symbol number, at size
 *     72/152.4 * symbol height. Symbols with no Unicode glyph are drawn from PLplot's Hershey strokes.
 *   - Lines use bevel joins and butt caps, single segments and dashes round caps. Dashed lines use
 *     PLplot's dash patterns.
 *   - Filled areas are filled and then outlined with the current line width.
 * Inside text and symbols, one unit is 1/72 inch ("text points"): this is what one Cairo user unit
 * was for PLplot's 72 dpi hardcopy pages.
 */

#if defined (CESR_CAIRO)

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pango/pangocairo.h>
#include "qp_cairo_render.h"
#include "qp_cairo_hershey.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Font size [pt] for a quick_plot text height [pt]. */
#define QPR_TEXT_FONT_SCALE (0.25 * 72.0 / 25.4)

/* Font size [pt] for a quick_plot symbol height [pt]. */
#define QPR_SYM_FONT_SCALE (72.0 / 152.4)

/* Line width unit [inch]. */
#define QPR_LINE_WIDTH_UNIT 0.005

/* One PLplot mm (dash lengths, hatch spacing) in points. */
#define QPR_MM (72.0 / 25.4)

/* Arrow head length [inch] for head_size = 1. */
#define QPR_ARROW_HEAD_LEN 0.12

#define QPR_MARKUP_LEN 8192

/*-------------------------------------------------------------------------*/

static void default_state (qpr_state *s) {
  memset(s, 0, sizeof(*s));
  s->line_width = 0;   /* PLplot's default: drawn 1 point wide. */
  s->line_pattern = QPR_SOLID;
  s->char_size = 12;
  s->sym_size = 10;
}

void qpr_init (qpr_renderer *r, cairo_t *cr, double scale, double ox, double oy, int draw) {
  memset(r, 0, sizeof(*r));
  r->cr = cr;
  r->scale = r->scale_x = r->scale_y = scale;
  r->ox = ox;
  r->oy = oy;
  r->draw = draw;
  r->page_w = 8.5;
  r->page_h = 11;
  r->depth = 0;
  default_state(&r->st[0]);
}

void qpr_set_view (qpr_renderer *r, double scale, double ox, double oy) {
  r->scale = r->scale_x = r->scale_y = scale;
  r->ox = ox;
  r->oy = oy;
}

void qpr_set_xy_scale (qpr_renderer *r, double scale_x, double scale_y) {
  r->scale_x = scale_x;
  r->scale_y = scale_y;
}

/*-------------------------------------------------------------------------*/
/* Coordinate helpers: inches (origin lower left) -> device (origin upper left). */

static double dev_x (qpr_renderer *r, double x) {return r->ox + x * r->scale_x;}
static double dev_y (qpr_renderer *r, double y) {return r->oy + (r->page_h - y) * r->scale_y;}
static double pt_to_dev (qpr_renderer *r, double pt) {return pt * r->scale / 72.0;}
static qpr_state *cur (qpr_renderer *r) {return &r->st[r->depth];}

static void set_source (qpr_renderer *r) {
  qpr_state *s = cur(r);
  cairo_set_source_rgb(r->cr, s->rgb[0], s->rgb[1], s->rgb[2]);
}

/* Line width in device units. Zero width lines would be invisible in Cairo: use 1 text point. */

static double line_width_dev (qpr_renderer *r) {
  qpr_state *s = cur(r);
  if (s->line_width <= 0) return pt_to_dev(r, 1);
  return s->line_width * QPR_LINE_WIDTH_UNIT * r->scale;
}

/* Solid line style: what PLplot uses for polylines and outlines. */

static void solid_line_style (qpr_renderer *r) {
  set_source(r);
  cairo_set_line_width(r->cr, line_width_dev(r));
  cairo_set_line_join(r->cr, CAIRO_LINE_JOIN_BEVEL);
  cairo_set_line_cap(r->cr, CAIRO_LINE_CAP_BUTT);
  cairo_set_dash(r->cr, NULL, 0, 0);
}

/* Line style for the current line pattern. PLplot's patterns (in mm, mark then space) are used
 * for the quick_plot patterns with the same number. PLplot draws each dash as a separate line
 * with round caps. */

static void line_style (qpr_renderer *r) {
  static const double pat[5][4] = {
    {0, 0, 0, 0},           /* solid */
    {1.0, 1.0, 0, 0},       /* dashed */
    {2.0, 2.0, 0, 0},       /* dash_dot */
    {2.0, 1.0, 0, 0},       /* dotted */
    {2.5, 1.0, 1.0, 1.0}    /* dash_dot3 */
  };
  static const int npat[5] = {0, 2, 2, 2, 4};
  int p = cur(r)->line_pattern - 1, i;
  double dash[4];

  solid_line_style(r);
  if (p < 1 || p > 4) return;
  for (i = 0; i < npat[p]; i++) dash[i] = pt_to_dev(r, pat[p][i] * QPR_MM);
  cairo_set_dash(r->cr, dash, npat[p], 0);
  cairo_set_line_cap(r->cr, CAIRO_LINE_CAP_ROUND);
}

/* Begin a primitive: save cairo state and clip to the graph if clipping is on. */

static void begin_clipped (qpr_renderer *r) {
  qpr_state *s = cur(r);
  cairo_save(r->cr);
  cairo_new_path(r->cr);
  if (s->clip) {
    double x1 = dev_x(r, s->graph[0]), x2 = dev_x(r, s->graph[1]);
    double y1 = dev_y(r, s->graph[2]), y2 = dev_y(r, s->graph[3]);
    cairo_rectangle(r->cr, fmin(x1, x2), fmin(y1, y2), fabs(x2-x1), fabs(y2-y1));
    cairo_clip(r->cr);
  }
}

static void end_clipped (qpr_renderer *r) {
  cairo_restore(r->cr);
}

/*-------------------------------------------------------------------------*/
/* Text. */

static int hershey_to_unicode (int code) {
  int lo = 0, hi = (int)(sizeof(qpr_hershey_table) / sizeof(qpr_hershey_table[0])) - 1;
  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    if (qpr_hershey_table[mid].hershey == code) return qpr_hershey_table[mid].unicode;
    if (qpr_hershey_table[mid].hershey < code) lo = mid + 1; else hi = mid - 1;
  }
  return 0;
}

/* Greek letter for \gX, using PLplot's mnemonic table and Hershey code substitutions. */

static int greek_unicode (char c) {
  static const char mnemonic[] = "ABGDEZYHIKLMNCOPRSTUFXQWabgdezyhiklmncoprstufxqw";
  const char *p = c ? strchr(mnemonic, c) : NULL;
  int ig;
  if (!p) return 0;
  ig = (int)(p - mnemonic);
  if (ig >= 24) ig = ig + 100 - 24;
  ig = ig + 527;
  if (ig == 631) ig = 684;
  else if (ig == 634) ig = 685;
  else if (ig == 647) ig = 686;
  return hershey_to_unicode(ig);
}

typedef struct {
  char *s;
  size_t n;
} mbuf;

static void m_cat (mbuf *m, const char *t) {
  size_t k = strlen(t);
  if (m->n + k + 1 > QPR_MARKUP_LEN) return;
  memcpy(m->s + m->n, t, k + 1);
  m->n += k;
}

static void m_unicode (mbuf *m, int u) {
  char b[8];
  int k = 0;
  if (u == '&') {m_cat(m, "&#38;"); return;}
  if (u == '<') {m_cat(m, "&#60;"); return;}
  if (u == '>') {m_cat(m, "&#62;"); return;}
  if (u <= 0) return;
  if (u < 0x80) {b[k++] = (char)u;}
  else if (u < 0x800) {b[k++] = (char)(0xC0 | (u >> 6)); b[k++] = (char)(0x80 | (u & 0x3F));}
  else if (u < 0x10000) {b[k++] = (char)(0xE0 | (u >> 12)); b[k++] = (char)(0x80 | ((u >> 6) & 0x3F));
                         b[k++] = (char)(0x80 | (u & 0x3F));}
  else {b[k++] = (char)(0xF0 | (u >> 18)); b[k++] = (char)(0x80 | ((u >> 12) & 0x3F));
        b[k++] = (char)(0x80 | ((u >> 6) & 0x3F)); b[k++] = (char)(0x80 | (u & 0x3F));}
  b[k] = 0;
  m_cat(m, b);
}

/* Font: family 0 = sans, 1 = serif, 3 = script. */

typedef struct {
  int family;
  int italic;
} font_t;

/* PLplot's plP_script_scale. */

static void script_scale (int ifupper, int *level, double *old_scale, double *scale,
                          double *old_offset, double *offset) {
  if (*level == 0) {
    *old_scale = 1;
    *old_offset = 0;
  } else {
    *old_scale = *scale;
    *old_offset = *offset;
  }
  if ((*level >= 0 && ifupper) || (*level <= 0 && !ifupper)) {
    *scale = 0.75 * *old_scale;
    *offset = *old_offset + *old_scale;
  } else {
    *scale = *old_scale / 0.75;
    *offset = *old_offset - *scale;
  }
  if (ifupper) (*level)++; else (*level)--;
}

static void rise_span (mbuf *m, int ifsuperscript, double font_size, double mult, double rise) {
  char tag[100];
  double off = 1024.0 * 0.5 * font_size * (1.0 - mult);
  rise = 1024.0 * font_size * 0.8 * rise;
  if (ifsuperscript) snprintf(tag, sizeof tag, "<span rise=\"%d\" size=\"%d\">", (int)(rise + off), (int)(font_size * 1024 * mult));
  else               snprintf(tag, sizeof tag, "<span rise=\"%d\" size=\"%d\">", (int)-(rise - off), (int)(font_size * 1024 * mult));
  m_cat(m, tag);
}

static void open_span (mbuf *m, font_t f, double font_size, int up_down) {
  static const char *family[] = {"sans", "serif", "monospace", "sans,serif", "sans,serif"};
  char tag[200];
  int i, level = 0;
  double old_scale = 1, scale = 1, old_off = 0, off = 0;

  snprintf(tag, sizeof tag, "<span font_desc=\"%s\" size=\"%d\" style=\"%s\" weight=\"normal\">",
           family[f.family], (int)(font_size * 1024), f.italic ? "italic" : "normal");
  m_cat(m, tag);
  for (i = 0; i < up_down; i++) {
    script_scale(1, &level, &old_scale, &scale, &old_off, &off);
    rise_span(m, 1, font_size, scale, off);
  }
  for (i = 0; i > up_down; i--) {
    script_scale(0, &level, &old_scale, &scale, &old_off, &off);
    rise_span(m, 0, font_size, scale, off);
  }
}

static void close_span (mbuf *m, int up_down) {
  while (up_down > 0) {m_cat(m, "</span>"); up_down--;}
  while (up_down < 0) {m_cat(m, "</span>"); up_down++;}
  m_cat(m, "</span>");
}

/* Decode one UTF-8 character starting at t[i]. Bytes that are not valid UTF-8 are taken as Latin-1.
 * Returns the number of bytes used. */

static int utf8_char (const unsigned char *t, int i, int n, int *u) {
  unsigned char c = t[i];
  int k, need, v;
  if (c < 0x80) {*u = c; return 1;}
  if ((c & 0xE0) == 0xC0)      {need = 1; v = c & 0x1F;}
  else if ((c & 0xF0) == 0xE0) {need = 2; v = c & 0x0F;}
  else if ((c & 0xF8) == 0xF0) {need = 3; v = c & 0x07;}
  else {*u = c; return 1;}
  for (k = 1; k <= need; k++) {
    if (i + k >= n || (t[i+k] & 0xC0) != 0x80) {*u = c; return 1;}
    v = (v << 6) | (t[i+k] & 0x3F);
  }
  *u = v;
  return need + 1;
}

/* Convert quick_plot text (PGPLOT escapes) to Pango markup. Follows what plplot_interface.f90,
 * PLplot's text encoding and the Cairo driver's markup conversion do together. */

/* True if text (blanks allowed at either end) is a number, EG "-0.5", "+3", "1.2E-03". */

static int is_number (const char *t, int n) {
  int i = 0, nd = 0;
  while (n > 0 && t[n-1] == ' ') n--;
  while (i < n && t[i] == ' ') i++;
  if (i < n && (t[i] == '-' || t[i] == '+')) i++;
  for (; i < n && isdigit((unsigned char)t[i]); i++) nd++;
  if (i < n && t[i] == '.') for (i++; i < n && isdigit((unsigned char)t[i]); i++) nd++;
  if (nd == 0) return 0;
  if (i < n && strchr("eEdD", t[i])) {
    int ne = 0;
    i++;
    if (i < n && (t[i] == '-' || t[i] == '+')) i++;
    for (; i < n && isdigit((unsigned char)t[i]); i++) ne++;
    if (ne == 0) return 0;
  }
  return i == n;
}

static void text_to_markup (const char *text_in, int n_in, double font_size, mbuf *m) {
  font_t f = {0, 0};
  int up_down = 0, level = 0, i = 0, n = 0;
  int number = is_number(text_in, n_in);
  double old_scale = 1, scale = 1, old_off = 0, off = 0;
  char text[2048];

  /* As in plplot_interface.f90: \. -> \u.\d and \x -> x.
   * In a number (EG an axis tick label) a hyphen is drawn as a minus sign (U+2212). */
  while (n_in > 0 && text_in[n_in-1] == ' ') n_in--;
  for (i = 0; i < n_in && n < (int)sizeof(text) - 8; i++) {
    if (number && text_in[i] == '-') {
      memcpy(text + n, "\xe2\x88\x92", 3); n += 3;
    } else if (text_in[i] == '\\' && i + 1 < n_in && text_in[i+1] == '.') {
      memcpy(text + n, "\\u.\\d", 5); n += 5; i++;
    } else if (text_in[i] == '\\' && i + 1 < n_in && text_in[i+1] == 'x') {
      text[n++] = 'x'; i++;
    } else {
      text[n++] = text_in[i];
    }
  }

  open_span(m, f, font_size, 0);
  i = 0;

  while (i < n) {
    int u;

    if (text[i] == '\\' && i + 1 < n) {
      char e = text[i+1];
      i += 2;
      switch (e) {
      case 'u':
        if (up_down < 0) {m_cat(m, "</span>"); level++;}
        else {script_scale(1, &level, &old_scale, &scale, &old_off, &off); rise_span(m, 1, font_size, scale, off);}
        up_down++;
        break;
      case 'd':
        if (up_down > 0) {m_cat(m, "</span>"); level--;}
        else {script_scale(0, &level, &old_scale, &scale, &old_off, &off); rise_span(m, 0, font_size, scale, off);}
        up_down--;
        break;
      case 'f': case 'F': {
        font_t g = {0, 0};
        char k = (i < n) ? text[i] : 0;
        if (k == 'n')      g.family = 0;
        else if (k == 'r') g.family = 1;
        else if (k == 'i') {g.family = 1; g.italic = 1;}
        else if (k == 's') g.family = 3;
        else break;
        i++;
        close_span(m, up_down);
        f = g;
        open_span(m, f, font_size, up_down);
        break;
      }
      case 'g': case 'G':
        if (i < n) m_unicode(m, greek_unicode(text[i]));
        i++;
        break;
      case '(': {
        int code = 0;
        while (i < n && text[i] >= '0' && text[i] <= '9') code = 10 * code + (text[i++] - '0');
        if (i < n && text[i] == ')') i++;
        m_unicode(m, hershey_to_unicode(code));
        break;
      }
      case 'A':                 /* Angstrom */
        m_unicode(m, 0xC5);
        break;
      case '\\':
        m_cat(m, "\\");
        break;
      default:                  /* Unknown escape (EG \b): dropped, as with PLplot. */
        break;
      }
      continue;
    }

    i += utf8_char((const unsigned char *)text, i, n, &u);
    m_unicode(m, u);
  }

  close_span(m, up_down);
}

/* Create a Pango layout for markup drawn in text point units. */

static PangoLayout *make_layout (cairo_t *cr, const char *markup) {
  PangoLayout *layout = pango_cairo_create_layout(cr);
  pango_layout_set_markup(layout, markup, -1);
  return layout;
}

/* Draw text with its baseline starting at (x0, y0) [inch].
 * justify: 0 = left, 0.5 = center, 1 = right. Angle is in degrees. */

static void draw_text (qpr_renderer *r, double x0, double y0, double angle, double justify,
                       const char *text, int n) {
  qpr_state *s = cur(r);
  double font_size = QPR_TEXT_FONT_SCALE * s->char_size;
  char markup[QPR_MARKUP_LEN];
  mbuf m = {markup, 0};
  PangoLayout *layout;
  int w, h, baseline;

  markup[0] = 0;
  text_to_markup(text, n, font_size, &m);

  begin_clipped(r);
  cairo_translate(r->cr, dev_x(r, x0), dev_y(r, y0));
  cairo_rotate(r->cr, -angle * M_PI / 180.0);
  cairo_scale(r->cr, r->scale / 72.0, r->scale / 72.0);
  layout = make_layout(r->cr, markup);
  pango_layout_get_pixel_size(layout, &w, &h);
  baseline = pango_layout_get_baseline(layout);
  set_source(r);
  cairo_move_to(r->cr, -justify * w, -baseline / 1024.0);
  pango_cairo_show_layout(r->cr, layout);
  g_object_unref(layout);
  end_clipped(r);
}

double qpr_text_len (qpr_renderer *r, const char *text, int n) {
  int i, n_text;
  while (n > 0 && text[n-1] == ' ') n--;
  n_text = n;
  for (i = 0; i < n; i++) if (text[i] == '\\') n_text--;
  return 0.8 * n_text * cur(r)->char_size * 0.25 / 25.4;
}

/*-------------------------------------------------------------------------*/
/* Symbols. */

/* Unicode glyph for each quick_plot symbol number. These are the glyphs PLplot draws (Hershey
 * symbol -> Unicode). 0 => drawn from Hershey strokes. triangle_filled (13) is a filled triangle;
 * PLplot's font table mapped it to a square. */

static const int symbol_glyph[19] = {
  0x25A1, 0x22C5, 0x002B, 0x2217, 0x25CB, 0x00D7, 0x25A1, 0x25B3, 0x2295, 0x2299,
  0,      0x2662, 0x2729, 0x25B2, 0,      0x2721, 0x25A0, 0x2022, 0x22C6
};

/* Hershey strokes for the symbols with no Unicode glyph. Pen-up is marked by 64. y is up. */

static const signed char hershey_square_concave[] = {
  -9,9, -8,7, -7,3, -7,-3, -8,-7, -9,-9, 64,0,  9,9, 8,7, 7,3, 7,-3, 8,-7, 9,-9, 64,0,
  -9,9, -7,8, -3,7, 3,7, 7,8, 9,9, 64,0,  -9,-9, -7,-8, -3,-7, 3,-7, 7,-8, 9,-9, 64,64};

static const signed char hershey_red_cross[] = {
  -2,6, -2,2, -6,2, -6,-2, -2,-2, -2,-6, 2,-6, 2,-2, 6,-2, 6,2, 2,2, 2,6, -2,6, 64,64};

static void draw_hershey (qpr_renderer *r, double xd, double yd, double font_size, const signed char *v) {
  double u = pt_to_dev(r, 0.05 * font_size);
  int k = 0, pen_up = 1;
  solid_line_style(r);
  cairo_new_path(r->cr);
  for (;;) {
    int cx = v[k], cy = v[k+1];
    k += 2;
    if (cx == 64 && cy == 64) break;
    if (cx == 64) {pen_up = 1; continue;}
    if (pen_up) cairo_move_to(r->cr, xd + u * cx, yd - u * cy);
    else        cairo_line_to(r->cr, xd + u * cx, yd - u * cy);
    pen_up = 0;
  }
  cairo_stroke(r->cr);
}

static void draw_symbol (qpr_renderer *r, double xi, double yi, int type) {
  qpr_state *s = cur(r);
  double font_size = QPR_SYM_FONT_SCALE * s->sym_size;
  double xd = dev_x(r, xi), yd = dev_y(r, yi);
  char markup[300];
  mbuf m = {markup, 0};
  font_t f = {0, 0};
  PangoLayout *layout;
  int w, h, baseline;

  if (type < 0 || type > 18) return;

  begin_clipped(r);

  if (type == 10 || type == 14) {
    draw_hershey(r, xd, yd, font_size, type == 10 ? hershey_square_concave : hershey_red_cross);
    end_clipped(r);
    return;
  }

  markup[0] = 0;
  open_span(&m, f, font_size, 0);
  m_unicode(&m, symbol_glyph[type]);
  close_span(&m, 0);

  cairo_translate(r->cr, xd, yd);
  cairo_scale(r->cr, r->scale / 72.0, r->scale / 72.0);
  layout = make_layout(r->cr, markup);
  pango_layout_get_pixel_size(layout, &w, &h);
  baseline = pango_layout_get_baseline(layout);
  set_source(r);
  cairo_move_to(r->cr, -0.5 * w, 0.5 * font_size - baseline / 1024.0);
  pango_cairo_show_layout(r->cr, layout);
  g_object_unref(layout);
  end_clipped(r);
}

/*-------------------------------------------------------------------------*/
/* Rectangles. */

static void hatch (qpr_renderer *r, double x1, double y1, double x2, double y2, double angle_deg) {
  cairo_t *cr = r->cr;
  double sp = pt_to_dev(r, 2.0 * QPR_MM);
  double a = angle_deg * M_PI / 180, ca = cos(a), sa = sin(a);
  double cx = 0.5 * (x1 + x2), cy = 0.5 * (y1 + y2);
  double rad = 0.5 * hypot(x2 - x1, y2 - y1), t;
  /* Lines along direction (ca, -sa) (device y is down), spaced along the normal. */
  for (t = -rad; t <= rad; t += sp) {
    double px = cx + t * sa, py = cy + t * ca;
    cairo_move_to(cr, px - rad * ca, py + rad * sa);
    cairo_line_to(cr, px + rad * ca, py - rad * sa);
  }
}

static void draw_rect (qpr_renderer *r, double xa, double xb, double ya, double yb, int fill) {
  cairo_t *cr = r->cr;
  double x1 = fmin(dev_x(r, xa), dev_x(r, xb)), x2 = fmax(dev_x(r, xa), dev_x(r, xb));
  double y1 = fmin(dev_y(r, ya), dev_y(r, yb)), y2 = fmax(dev_y(r, ya), dev_y(r, yb));

  begin_clipped(r);
  solid_line_style(r);

  switch (fill) {
  case QPR_NO_FILL:
    cairo_rectangle(cr, x1, y1, x2-x1, y2-y1);
    cairo_stroke(cr);
    break;
  case QPR_HATCHED:
  case QPR_CROSS_HATCHED:
    cairo_rectangle(cr, x1, y1, x2-x1, y2-y1);
    cairo_clip(cr);
    cairo_new_path(cr);
    hatch(r, x1, y1, x2, y2, 45);
    if (fill == QPR_CROSS_HATCHED) hatch(r, x1, y1, x2, y2, -45);
    cairo_stroke(cr);
    break;
  default:
    cairo_rectangle(cr, x1, y1, x2-x1, y2-y1);
    cairo_fill_preserve(cr);
    cairo_stroke(cr);
  }
  end_clipped(r);
}

/*-------------------------------------------------------------------------*/

/* Points are read from the (possibly unaligned) stream payload. */

static void draw_polyline (qpr_renderer *r, int n, const uint8_t *px, const uint8_t *py) {
  int i;
  double x, y;
  if (n < 2) return;
  begin_clipped(r);
  line_style(r);
  /* PLplot draws a single segment as a line with round caps. */
  if (n == 2) cairo_set_line_cap(r->cr, CAIRO_LINE_CAP_ROUND);
  for (i = 0; i < n; i++) {
    memcpy(&x, px + i*sizeof(double), sizeof(double));
    memcpy(&y, py + i*sizeof(double), sizeof(double));
    if (i == 0) cairo_move_to(r->cr, dev_x(r, x), dev_y(r, y));
    else        cairo_line_to(r->cr, dev_x(r, x), dev_y(r, y));
  }
  cairo_stroke(r->cr);
  end_clipped(r);
}

static void draw_arrow (qpr_renderer *r, const double *a, int filled) {
  cairo_t *cr = r->cr;
  double x1 = dev_x(r, a[0]), y1 = dev_y(r, a[1]), x2 = dev_x(r, a[2]), y2 = dev_y(r, a[3]);
  double half = 0.5 * a[4] * M_PI / 180, barb = a[5];
  double len = a[6] * QPR_ARROW_HEAD_LEN * r->scale;
  double dx = x2 - x1, dy = y2 - y1, d = sqrt(dx*dx + dy*dy);
  double ux, uy, nx, ny, w;

  if (d == 0) return;
  ux = dx / d; uy = dy / d; nx = -uy; ny = ux;
  w = len * tan(half);

  begin_clipped(r);
  line_style(r);
  cairo_move_to(cr, x1, y1);
  cairo_line_to(cr, x2 - ux * len * (1 - barb), y2 - uy * len * (1 - barb));
  cairo_stroke(cr);

  solid_line_style(r);
  cairo_move_to(cr, x2, y2);
  cairo_line_to(cr, x2 - ux*len + nx*w, y2 - uy*len + ny*w);
  cairo_line_to(cr, x2 - ux*len*(1-barb), y2 - uy*len*(1-barb));
  cairo_line_to(cr, x2 - ux*len - nx*w, y2 - uy*len - ny*w);
  cairo_close_path(cr);
  if (filled) cairo_fill(cr); else cairo_stroke(cr);
  end_clipped(r);
}

/*-------------------------------------------------------------------------*/
/* Stream decoding. */

typedef struct {
  const uint8_t *p;
  size_t n;
  int bad;
} reader;

static double rd_d (reader *rd) {
  double v = 0;
  if (rd->n < sizeof(double)) {rd->bad = 1; return 0;}
  memcpy(&v, rd->p, sizeof(double));
  rd->p += sizeof(double); rd->n -= sizeof(double);
  return v;
}

static int32_t rd_i (reader *rd) {
  int32_t v = 0;
  if (rd->n < sizeof(int32_t)) {rd->bad = 1; return 0;}
  memcpy(&v, rd->p, sizeof(int32_t));
  rd->p += sizeof(int32_t); rd->n -= sizeof(int32_t);
  return v;
}

int qpr_exec (qpr_renderer *r, const uint8_t *buf, size_t len) {
  size_t off = 0;

  while (off < len) {
    uint8_t op;
    uint32_t plen;
    reader rd;
    qpr_state *s = cur(r);

    if (off + 5 > len) return -1;
    op = buf[off];
    memcpy(&plen, buf + off + 1, 4);
    if (off + 5 + plen > len) return -1;
    rd.p = buf + off + 5; rd.n = plen; rd.bad = 0;
    off += 5 + plen;

    switch (op) {
    case QPR_OP_PAGE:
      r->page_w = rd_d(&rd);
      r->page_h = rd_d(&rd);
      break;
    case QPR_OP_CLEAR:
      if (r->draw) {
        cairo_save(r->cr);
        cairo_set_source_rgb(r->cr, 1, 1, 1);
        cairo_rectangle(r->cr, r->ox, r->oy, r->page_w * r->scale_x, r->page_h * r->scale_y);
        cairo_fill(r->cr);
        cairo_restore(r->cr);
      }
      break;
    case QPR_OP_STATE:
      if (plen >= sizeof(qpr_state)) memcpy(s, rd.p, sizeof(qpr_state));
      break;
    case QPR_OP_PUSH:
      if (r->depth < QPR_MAX_DEPTH - 1) {r->st[r->depth+1] = r->st[r->depth]; r->depth++;}
      break;
    case QPR_OP_POP:
      if (r->depth > 0) r->depth--;
      break;
    case QPR_OP_COLOR:
      s->rgb[0] = rd_d(&rd); s->rgb[1] = rd_d(&rd); s->rgb[2] = rd_d(&rd);
      break;
    case QPR_OP_LINE_WIDTH:   s->line_width = rd_d(&rd); break;
    case QPR_OP_LINE_PATTERN: s->line_pattern = rd_i(&rd); break;
    case QPR_OP_CLIP:         s->clip = rd_i(&rd); break;
    case QPR_OP_GRAPH_POS:
      s->graph[0] = rd_d(&rd); s->graph[1] = rd_d(&rd); s->graph[2] = rd_d(&rd); s->graph[3] = rd_d(&rd);
      break;
    case QPR_OP_CHAR_SIZE:    s->char_size = rd_d(&rd); break;
    case QPR_OP_SYM_SIZE:     s->sym_size = rd_d(&rd); break;
    case QPR_OP_POLYLINE: {
      int32_t n = rd_i(&rd);
      if (n < 0 || rd.n < 2 * (size_t)n * sizeof(double)) return -1;
      if (r->draw) draw_polyline(r, n, rd.p, rd.p + n * sizeof(double));
      break;
    }
    case QPR_OP_RECT: {
      double x1 = rd_d(&rd), x2 = rd_d(&rd), y1 = rd_d(&rd), y2 = rd_d(&rd);
      int32_t fill = rd_i(&rd);
      if (r->draw && !rd.bad) draw_rect(r, x1, x2, y1, y2, fill);
      break;
    }
    case QPR_OP_SYMBOL: {
      double x = rd_d(&rd), y = rd_d(&rd);
      int32_t type = rd_i(&rd);
      if (r->draw && !rd.bad) draw_symbol(r, x, y, type);
      break;
    }
    case QPR_OP_TEXT: {
      double x = rd_d(&rd), y = rd_d(&rd), ang = rd_d(&rd), just = rd_d(&rd);
      int32_t n = rd_i(&rd);
      if (rd.bad || n < 0 || (size_t)n > rd.n) return -1;
      if (r->draw) draw_text(r, x, y, ang, just, (const char *)rd.p, n);
      break;
    }
    case QPR_OP_ARROW: {
      double a[7];
      int i;
      int32_t filled;
      for (i = 0; i < 7; i++) a[i] = rd_d(&rd);
      filled = rd_i(&rd);
      if (r->draw && !rd.bad) draw_arrow(r, a, filled);
      break;
    }
    default:
      break;   /* Unknown op: skip. */
    }
    if (rd.bad) return -1;
  }
  return 0;
}

#endif
