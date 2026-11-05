/*
 * qp_cairo_render.c
 *
 * Shared drawing engine for the quick_plot Cairo backend. See qp_cairo_render.h.
 * Compiled into sim_utils (when CESR_CAIRO is defined) and into qp_cairo_viewer.
 */

#if defined (CESR_CAIRO) || defined (QP_CAIRO_VIEWER)

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "qp_cairo_render.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Inches spacing of hatch lines. */
#define QPR_HATCH_SPACING 0.06

/* Arrow head length [inch] for head_size = 1. */
#define QPR_ARROW_HEAD_LEN 0.12

/* Symbol half-width as a fraction of the symbol height. */
#define QPR_SYM_RADIUS_FACTOR 0.30

/* Super/subscript size reduction and baseline shift (fraction of current height). */
#define QPR_SCRIPT_SCALE 0.7
#define QPR_SCRIPT_RISE  0.5

/* Font em size as a fraction of the quick_plot text height. Matches the PLPLOT backend's text size. */
#define QPR_FONT_SCALE 0.82

/* Line width unit [inch]. Matches the PLPLOT backend's line weights. */
#define QPR_LINE_WIDTH_UNIT (1.0/400.0)

/*-------------------------------------------------------------------------*/

void qpr_default_state (qpr_state *s) {
  memset(s, 0, sizeof(*s));
  s->line_width = 1;
  s->line_pattern = QPR_SOLID;
  s->char_size = 12;
  s->sym_size = 10;
  s->sym_fill = QPR_SOLID_FILL;
  s->text_bg_opaque = 0;
  s->text_bg_rgb[0] = s->text_bg_rgb[1] = s->text_bg_rgb[2] = 1;
}

void qpr_init (qpr_renderer *r, cairo_t *cr, double scale, double ox, double oy, int draw) {
  memset(r, 0, sizeof(*r));
  r->cr = cr;
  r->scale = scale;
  r->ox = ox;
  r->oy = oy;
  r->draw = draw;
  r->page_w = 8.5;
  r->page_h = 11;
  r->depth = 0;
  qpr_default_state(&r->st[0]);
}

void qpr_set_view (qpr_renderer *r, double scale, double ox, double oy) {
  r->scale = scale;
  r->ox = ox;
  r->oy = oy;
}

/*-------------------------------------------------------------------------*/
/* Coordinate helpers: inches (origin lower left) -> device (origin upper left). */

static double dev_x (qpr_renderer *r, double x) {return r->ox + x * r->scale;}
static double dev_y (qpr_renderer *r, double y) {return r->oy + (r->page_h - y) * r->scale;}
static double pt_to_dev (qpr_renderer *r, double pt) {return pt * r->scale / 72.0;}
static qpr_state *cur (qpr_renderer *r) {return &r->st[r->depth];}

static void set_source (qpr_renderer *r) {
  qpr_state *s = cur(r);
  cairo_set_source_rgb(r->cr, s->rgb[0], s->rgb[1], s->rgb[2]);
}

static void apply_line_style (qpr_renderer *r) {
  qpr_state *s = cur(r);
  double lw = s->line_width * r->scale * QPR_LINE_WIDTH_UNIT;
  double u;
  double dash[8];
  int nd = 0;

  if (lw < r->min_line_width) lw = r->min_line_width;
  set_source(r);
  cairo_set_line_width(r->cr, lw);
  cairo_set_line_join(r->cr, CAIRO_LINE_JOIN_ROUND);
  cairo_set_line_cap(r->cr, CAIRO_LINE_CAP_BUTT);

  /* Dash lengths scale with line width a little so thick dashed lines still look dashed. */
  u = pt_to_dev(r, 1.0) * (1 + 0.15 * (s->line_width > 1 ? s->line_width - 1 : 0));

  switch (s->line_pattern) {
  case QPR_DASHED:    dash[0] = 6*u; dash[1] = 4*u; nd = 2; break;
  case QPR_DASH_DOT:  dash[0] = 6*u; dash[1] = 3*u; dash[2] = 1*u; dash[3] = 3*u; nd = 4; break;
  case QPR_DOTTED:    dash[0] = 1*u; dash[1] = 3*u; nd = 2; break;
  case QPR_DASH_DOT3: dash[0] = 6*u; dash[1] = 3*u; dash[2] = 1*u; dash[3] = 3*u;
                      dash[4] = 1*u; dash[5] = 3*u; dash[6] = 1*u; dash[7] = 3*u; nd = 8; break;
  default: nd = 0;
  }
  cairo_set_dash(r->cr, dash, nd, 0);
}

/* Begin a primitive: save cairo state and clip to the graph if clipping is on. */

static void begin_clipped (qpr_renderer *r) {
  qpr_state *s = cur(r);
  cairo_save(r->cr);
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
/* Text. Implements the PGPLOT escape sequences that quick_plot users rely on:
 *   \u \d      super/subscript (nestable)     \b       backspace one character
 *   \fn \fr \fi \fs   normal/roman/italic/script font
 *   \gX        Greek letter X (PGPLOT letter mapping)
 *   \A  Angstrom   \.  centered dot   \x  multiplication sign   \\  backslash
 *   \(nnnn)    Hershey character: a few common ones are mapped, others dropped.
 * Bytes >= 128 are taken to be Latin-1 and converted to UTF-8 so cairo never sees invalid UTF-8.
 */

static const char *greek_utf8 (char c) {
  /* PGPLOT mapping: A B G D E Z Y H I K L M N C O P R S T U F X Q W */
  static const char latin[] = "ABGDEZYHIKLMNCOPRSTUFXQW";
  static const char *upper[] = {"Α","Β","Γ","Δ","Ε","Ζ","Η","Θ","Ι","Κ","Λ","Μ","Ν","Ξ","Ο","Π","Ρ","Σ","Τ","Υ","Φ","Χ","Ψ","Ω"};
  static const char *lower[] = {"α","β","γ","δ","ε","ζ","η","θ","ι","κ","λ","μ","ν","ξ","ο","π","ρ","σ","τ","υ","φ","χ","ψ","ω"};
  const char *p = strchr(latin, toupper((unsigned char)c));
  if (!p || c == 0) return NULL;
  return isupper((unsigned char)c) ? upper[p - latin] : lower[p - latin];
}

static const char *hershey_utf8 (int code) {
  switch (code) {
  case 2245: return "∝";
  case 2246: return "∞";
  case 2248: return "≈";
  case 2239: return "≠";
  case 2242: return "≤";
  case 2243: return "≥";
  case 2233: return "±";
  case 2235: return "×";
  case 2237: return "÷";
  case 2266: return "√";
  case 2267: return "∫";
  case 2268: return "∮";
  case 2270: return "∂";
  case 2271: return "∇";
  case 2236: return "·";
  case 2278: return "†";
  case 2721: return "°";
  default:   return NULL;
  }
}

static void symbol_geometry (cairo_t *cr, double x, double y, double rad, int type, int fill);

typedef struct {
  qpr_renderer *r;
  int draw;
  double h;          /* base height, device units */
  double x;          /* current pen position along the baseline */
  double last_adv;
  int level;
  char font;
} text_ctx;

static void select_font (text_ctx *t) {
  cairo_font_slant_t slant = CAIRO_FONT_SLANT_NORMAL;
  const char *family = "sans-serif";
  switch (t->font) {
  case 'r': family = "serif"; break;
  case 'i': family = "serif"; slant = CAIRO_FONT_SLANT_ITALIC; break;
  case 's': family = "serif"; slant = CAIRO_FONT_SLANT_ITALIC; break;
  default: break;
  }
  cairo_select_font_face(t->r->cr, family, slant, CAIRO_FONT_WEIGHT_NORMAL);
}

static void emit_utf8 (text_ctx *t, const char *s) {
  cairo_text_extents_t ext;
  double size = t->h * pow(QPR_SCRIPT_SCALE, abs(t->level));
  double rise = 0;
  int i, lev = t->level;

  /* Each level shifts the baseline by a fraction of the height of the level it is attached to. */
  for (i = 0; i < abs(lev); i++) rise += QPR_SCRIPT_RISE * t->h * pow(QPR_SCRIPT_SCALE, i);
  if (lev < 0) rise = -rise;

  select_font(t);
  cairo_set_font_size(t->r->cr, QPR_FONT_SCALE * size);
  cairo_text_extents(t->r->cr, s, &ext);
  if (t->draw) {
    cairo_move_to(t->r->cr, t->x, -rise);
    cairo_show_text(t->r->cr, s);
  }
  t->x += ext.x_advance;
  t->last_adv = ext.x_advance;
}

static void emit_marker (text_ctx *t, int type) {
  double size = t->h * pow(QPR_SCRIPT_SCALE, abs(t->level));
  double rad = 0.3 * size;
  if (t->draw) {
    cairo_save(t->r->cr);
    cairo_set_dash(t->r->cr, NULL, 0, 0);
    symbol_geometry(t->r->cr, t->x + 1.3*rad, -0.35*size, rad, type, 1);
    cairo_restore(t->r->cr);
  }
  t->x += 2.6 * rad;
  t->last_adv = 2.6 * rad;
}

/* Walk the text, measuring (and drawing if t->draw). Returns the total advance in device units. */

static double layout_text (text_ctx *t, const char *text, int n) {
  int i = 0;
  char buf[8];

  while (i < n) {
    unsigned char c = (unsigned char)text[i];

    if (c == '\\' && i+1 < n) {
      char e = text[i+1];
      const char *g;
      switch (e) {
      case 'u': case 'U': t->level++; i += 2; continue;
      case 'd': case 'D': t->level--; i += 2; continue;
      case 'b': case 'B': t->x -= t->last_adv; i += 2; continue;
      case 'f': case 'F':
        if (i+2 < n) t->font = (char)tolower((unsigned char)text[i+2]);
        i += 3; continue;
      case 'g': case 'G':
        if (i+2 < n && (g = greek_utf8(text[i+2])) != NULL) emit_utf8(t, g);
        i += 3; continue;
      case 'A': emit_utf8(t, "Å"); i += 2; continue;
      case '.': emit_utf8(t, "·"); i += 2; continue;
      case 'x': emit_utf8(t, "×"); i += 2; continue;
      case '\\': emit_utf8(t, "\\"); i += 2; continue;
      case '(': {
        int j = i + 2, code = 0;
        while (j < n && isdigit((unsigned char)text[j])) code = 10*code + (text[j++] - '0');
        if (j < n && text[j] == ')') j++;
        if ((g = hershey_utf8(code)) != NULL) emit_utf8(t, g);
        i = j; continue;
      }
      case 'm': case 'M': {      /* \mn or \m(n): Graph marker n drawn inline. */
        int j = i + 2, code = 0, paren = 0;
        if (j < n && text[j] == '(') {paren = 1; j++;}
        while (j < n && isdigit((unsigned char)text[j])) code = 10*code + (text[j++] - '0');
        if (paren && j < n && text[j] == ')') j++;
        emit_marker(t, code);
        i = j; continue;
      }
      default: break;           /* Not an escape: draw the backslash literally. */
      }
    }

    if (c < 128) {
      buf[0] = (char)c; buf[1] = 0;
    } else {                    /* Latin-1 -> UTF-8 */
      buf[0] = (char)(0xC0 | (c >> 6));
      buf[1] = (char)(0x80 | (c & 0x3F));
      buf[2] = 0;
    }
    emit_utf8(t, buf);
    i++;
  }
  return t->x;
}

static double text_width_dev (qpr_renderer *r, const char *text, int n) {
  text_ctx t = {r, 0, pt_to_dev(r, cur(r)->char_size), 0, 0, 0, 'n'};
  return layout_text(&t, text, n);
}

double qpr_text_len (qpr_renderer *r, const char *text, int n) {
  return text_width_dev(r, text, n) / r->scale;
}

static void draw_text (qpr_renderer *r, double x0, double y0, double angle, double justify,
                       const char *text, int n) {
  qpr_state *s = cur(r);
  double w = text_width_dev(r, text, n);
  double h = pt_to_dev(r, s->char_size);
  text_ctx t = {r, 1, h, 0, 0, 0, 'n'};

  begin_clipped(r);
  cairo_translate(r->cr, dev_x(r, x0), dev_y(r, y0));
  cairo_rotate(r->cr, -angle * M_PI / 180.0);
  cairo_translate(r->cr, -justify * w, 0);

  if (s->text_bg_opaque) {
    cairo_set_source_rgb(r->cr, s->text_bg_rgb[0], s->text_bg_rgb[1], s->text_bg_rgb[2]);
    cairo_rectangle(r->cr, 0, -0.8*h, w, 1.05*h);
    cairo_fill(r->cr);
  }

  set_source(r);
  layout_text(&t, text, n);
  end_clipped(r);
}

/*-------------------------------------------------------------------------*/
/* Symbols. Drawn as geometry so all sizes are consistent. Types follow quick_plot. */

static void star_path (cairo_t *cr, double x, double y, double rad, int npt, double inner) {
  int i;
  for (i = 0; i < 2*npt; i++) {
    double a = -M_PI/2 + i * M_PI / npt;
    double rr = (i % 2 == 0) ? rad : rad * inner;
    if (i == 0) cairo_move_to(cr, x + rr*cos(a), y + rr*sin(a));
    else        cairo_line_to(cr, x + rr*cos(a), y + rr*sin(a));
  }
  cairo_close_path(cr);
}

static void triangle_path (cairo_t *cr, double x, double y, double rad, int up) {
  double s = up ? 1 : -1;
  cairo_move_to(cr, x, y - s*rad);
  cairo_line_to(cr, x + 0.866*rad, y + s*0.5*rad);
  cairo_line_to(cr, x - 0.866*rad, y + s*0.5*rad);
  cairo_close_path(cr);
}

/* Draw symbol type centered at (x, y) [device] using the current source and line width. */

static void symbol_geometry (cairo_t *cr, double x, double y, double rad, int type, int fill) {
  cairo_new_path(cr);

  switch (type) {
  case 0:  /* square */
    cairo_rectangle(cr, x-rad, y-rad, 2*rad, 2*rad); cairo_stroke(cr); break;
  case 1:  /* dot */
    cairo_arc(cr, x, y, 0.35*rad, 0, 2*M_PI); cairo_fill(cr); break;
  case 2:  /* plus */
    cairo_move_to(cr, x-rad, y); cairo_line_to(cr, x+rad, y);
    cairo_move_to(cr, x, y-rad); cairo_line_to(cr, x, y+rad); cairo_stroke(cr); break;
  case 3:  /* times (asterisk) */
    { int i; for (i = 0; i < 3; i++) {
        double a = i * M_PI / 3 + M_PI/2;
        cairo_move_to(cr, x - rad*cos(a), y - rad*sin(a));
        cairo_line_to(cr, x + rad*cos(a), y + rad*sin(a));
      } cairo_stroke(cr); } break;
  case 4:  /* circle */
    cairo_arc(cr, x, y, rad, 0, 2*M_PI); cairo_stroke(cr); break;
  case 5:  /* x */
    cairo_move_to(cr, x-rad, y-rad); cairo_line_to(cr, x+rad, y+rad);
    cairo_move_to(cr, x-rad, y+rad); cairo_line_to(cr, x+rad, y-rad); cairo_stroke(cr); break;
  case 7:  /* triangle */
    triangle_path(cr, x, y, 1.15*rad, 1); cairo_stroke(cr); break;
  case 8:  /* circle_plus */
    cairo_arc(cr, x, y, rad, 0, 2*M_PI);
    cairo_move_to(cr, x-rad, y); cairo_line_to(cr, x+rad, y);
    cairo_move_to(cr, x, y-rad); cairo_line_to(cr, x, y+rad); cairo_stroke(cr); break;
  case 9:  /* circle_dot */
    cairo_arc(cr, x, y, rad, 0, 2*M_PI); cairo_stroke(cr);
    cairo_arc(cr, x, y, 0.25*rad, 0, 2*M_PI); cairo_fill(cr); break;
  case 10: /* square_concave */
    cairo_move_to(cr, x-rad, y-rad);
    cairo_curve_to(cr, x-0.3*rad, y, x-0.3*rad, y, x-rad, y+rad);
    cairo_curve_to(cr, x, y+0.3*rad, x, y+0.3*rad, x+rad, y+rad);
    cairo_curve_to(cr, x+0.3*rad, y, x+0.3*rad, y, x+rad, y-rad);
    cairo_curve_to(cr, x, y-0.3*rad, x, y-0.3*rad, x-rad, y-rad);
    cairo_close_path(cr); cairo_stroke(cr); break;
  case 11: /* diamond */
    cairo_move_to(cr, x, y-1.2*rad); cairo_line_to(cr, x+rad, y);
    cairo_line_to(cr, x, y+1.2*rad); cairo_line_to(cr, x-rad, y);
    cairo_close_path(cr); cairo_stroke(cr); break;
  case 12: /* star5 */
    star_path(cr, x, y, 1.2*rad, 5, 0.4); cairo_stroke(cr); break;
  case 13: /* triangle_filled */
    triangle_path(cr, x, y, 1.15*rad, 1); if (fill) cairo_fill(cr); else cairo_stroke(cr); break;
  case 14: /* red_cross */
    cairo_move_to(cr, x-rad, y); cairo_line_to(cr, x+rad, y);
    cairo_move_to(cr, x, y-rad); cairo_line_to(cr, x, y+rad);
    cairo_set_line_width(cr, 2.5 * cairo_get_line_width(cr)); cairo_stroke(cr); break;
  case 15: /* star_of_david */
    triangle_path(cr, x, y, 1.15*rad, 1); triangle_path(cr, x, y, 1.15*rad, 0); cairo_stroke(cr); break;
  case 16: /* square_filled */
    cairo_rectangle(cr, x-rad, y-rad, 2*rad, 2*rad); if (fill) cairo_fill(cr); else cairo_stroke(cr); break;
  case 17: /* circle_filled */
    cairo_arc(cr, x, y, rad, 0, 2*M_PI); if (fill) cairo_fill(cr); else cairo_stroke(cr); break;
  case 18: /* star5_filled */
    star_path(cr, x, y, 1.2*rad, 5, 0.4); if (fill) cairo_fill(cr); else cairo_stroke(cr); break;
  default:
    if (type >= 32 && type < 127) {   /* A character used as a marker. */
      char str[2] = {(char)type, 0};
      cairo_text_extents_t ext;
      cairo_select_font_face(cr, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
      cairo_set_font_size(cr, 2.5*rad);
      cairo_text_extents(cr, str, &ext);
      cairo_move_to(cr, x - ext.x_bearing - ext.width/2, y - ext.y_bearing - ext.height/2);
      cairo_show_text(cr, str);
    } else if (type >= 0) {           /* Hershey symbol numbers etc: fall back to a circle. */
      cairo_arc(cr, x, y, rad, 0, 2*M_PI); cairo_stroke(cr);
    }                                 /* Negative: do_not_draw. */
  }
}

static void draw_symbol (qpr_renderer *r, double xi, double yi, int type) {
  qpr_state *s = cur(r);
  if (type < 0) return;
  begin_clipped(r);
  apply_line_style(r);
  cairo_set_dash(r->cr, NULL, 0, 0);
  symbol_geometry(r->cr, dev_x(r, xi), dev_y(r, yi), QPR_SYM_RADIUS_FACTOR * pt_to_dev(r, s->sym_size),
                  type, s->sym_fill != QPR_NO_FILL);
  end_clipped(r);
}

/*-------------------------------------------------------------------------*/

static void hatch (qpr_renderer *r, double x1, double y1, double x2, double y2, int dir) {
  cairo_t *cr = r->cr;
  double sp = QPR_HATCH_SPACING * r->scale;
  double w = x2 - x1, h = y2 - y1, c;
  for (c = -h; c < w; c += sp) {
    if (dir > 0) {cairo_move_to(cr, x1 + c, y2); cairo_line_to(cr, x1 + c + h, y1);}
    else         {cairo_move_to(cr, x1 + c, y1); cairo_line_to(cr, x1 + c + h, y2);}
  }
}

static void draw_rect (qpr_renderer *r, double xa, double xb, double ya, double yb, int fill) {
  cairo_t *cr = r->cr;
  double x1 = fmin(dev_x(r, xa), dev_x(r, xb)), x2 = fmax(dev_x(r, xa), dev_x(r, xb));
  double y1 = fmin(dev_y(r, ya), dev_y(r, yb)), y2 = fmax(dev_y(r, ya), dev_y(r, yb));

  /* Like PGPLOT, rectangles are not clipped to the graph. */
  cairo_save(cr);
  apply_line_style(r);
  cairo_new_path(cr);

  switch (fill) {
  case QPR_NO_FILL:
    cairo_rectangle(cr, x1, y1, x2-x1, y2-y1);
    cairo_stroke(cr);
    break;
  case QPR_HATCHED:
  case QPR_CROSS_HATCHED:
    cairo_rectangle(cr, x1, y1, x2-x1, y2-y1);
    cairo_clip(cr);
    cairo_set_dash(cr, NULL, 0, 0);
    hatch(r, x1, y1, x2, y2, 1);
    if (fill == QPR_CROSS_HATCHED) hatch(r, x1, y1, x2, y2, -1);
    cairo_stroke(cr);
    break;
  default:
    cairo_set_antialias(cr, CAIRO_ANTIALIAS_NONE);   /* Abutting fills must not show seams. */
    cairo_rectangle(cr, x1, y1, x2-x1, y2-y1);
    cairo_fill(cr);
  }
  cairo_restore(cr);
}

/* Points are read from the (possibly unaligned) stream payload. */

static void draw_polyline (qpr_renderer *r, int n, const uint8_t *px, const uint8_t *py) {
  int i;
  double x, y;
  if (n < 2) return;
  begin_clipped(r);
  apply_line_style(r);
  cairo_new_path(r->cr);
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
  apply_line_style(r);
  cairo_new_path(cr);
  cairo_move_to(cr, x1, y1);
  cairo_line_to(cr, x2 - ux * len * (1 - barb), y2 - uy * len * (1 - barb));
  cairo_stroke(cr);

  cairo_set_dash(cr, NULL, 0, 0);
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

int qpr_peek_page (const uint8_t *buf, size_t len, double *page_w, double *page_h) {
  size_t off = 0;
  while (off + 5 <= len) {
    uint32_t plen;
    memcpy(&plen, buf + off + 1, 4);
    if (off + 5 + plen > len) return -1;
    if (buf[off] == QPR_OP_PAGE && plen >= 2*sizeof(double)) {
      memcpy(page_w, buf + off + 5, sizeof(double));
      memcpy(page_h, buf + off + 5 + sizeof(double), sizeof(double));
      return 0;
    }
    off += 5 + plen;
  }
  return -1;
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
        cairo_rectangle(r->cr, r->ox, r->oy, r->page_w * r->scale, r->page_h * r->scale);
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
    case QPR_OP_SYM_FILL:     s->sym_fill = rd_i(&rd); break;
    case QPR_OP_TEXT_BG:
      s->text_bg_opaque = rd_i(&rd);
      s->text_bg_rgb[0] = rd_d(&rd); s->text_bg_rgb[1] = rd_d(&rd); s->text_bg_rgb[2] = rd_d(&rd);
      break;
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
