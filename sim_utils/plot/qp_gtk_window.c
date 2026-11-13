/*
 * qp_gtk_window.c
 *
 * GTK 3 + Cairo plot window for quick_plot. See qp_gtk_window.h.
 */

#if defined (CESR_GTK)

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <poll.h>
#include <math.h>
#include <dlfcn.h>
#include <gtk/gtk.h>
#include <readline/readline.h>
#include "qp_gtk_window.h"

/* In get_tty_char_c.c: called while waiting for a key in single character mode. */
extern void (*acc_tty_idle_hook) (void);

static struct {
  int initialized;
  int no_display;
  GtkWidget *win;
  GtkWidget *area;
  cairo_surface_t *surf;      /* Off-screen picture, device pixels, device scale set. */
  int surf_w, surf_h;         /* Picture size in logical pixels. */
  int surf_scale;
  int want_w, want_h;         /* Window size to use when (re)opening. */
  char title[256];
  qpg_resize_fn resize_fn;
  void *resize_data;
  guint resize_idle;          /* Pending resize redraw (coalesces resize events). */
  cairo_t *page_cr;           /* Context given to PLplot for the open page. */
  rl_getc_func_t *old_getc;   /* readline getc function in place before ours. */
  int hook_installed;
  GdkPixbuf *icon;
} g;

/*-------------------------------------------------------------------------*/

static int scale_factor (void) {
  return g.area ? gtk_widget_get_scale_factor(g.area) : 1;
}

/* Make sure the off-screen surface has the given size at the current screen scale. */

static void ensure_surface (int w, int h) {
  int s = scale_factor();
  if (g.page_cr) return;   /* PLplot is drawing into the current surface: keep it. */
  if (w < 1) w = 1;
  if (h < 1) h = 1;
  if (g.surf && g.surf_w == w && g.surf_h == h && g.surf_scale == s) return;
  if (g.surf) cairo_surface_destroy(g.surf);
  g.surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, w * s, h * s);
  cairo_surface_set_device_scale(g.surf, s, s);
  g.surf_w = w;
  g.surf_h = h;
  g.surf_scale = s;
}

/*-------------------------------------------------------------------------*/
/* Signal handlers. */

static gboolean on_draw (GtkWidget *w, cairo_t *cr, gpointer data) {
  int aw = gtk_widget_get_allocated_width(w), ah = gtk_widget_get_allocated_height(w);

  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_paint(cr);
  if (!g.surf) return TRUE;

  if (aw == g.surf_w && ah == g.surf_h) {
    cairo_set_source_surface(cr, g.surf, 0, 0);
  } else {
    /* Window and picture sizes differ (resize pending or no resize callback): scale to fit. */
    double sx = (double)aw / g.surf_w, sy = (double)ah / g.surf_h;
    double sc = sx < sy ? sx : sy;
    cairo_translate(cr, 0.5 * (aw - sc * g.surf_w), 0.5 * (ah - sc * g.surf_h));
    cairo_scale(cr, sc, sc);
    cairo_set_source_surface(cr, g.surf, 0, 0);
    cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_GOOD);
  }
  cairo_paint(cr);
  return TRUE;
}

static gboolean do_resize (gpointer data) {
  g.resize_idle = 0;
  if (g.area && g.resize_fn) {
    int aw = gtk_widget_get_allocated_width(g.area), ah = gtk_widget_get_allocated_height(g.area);
    if (aw != g.surf_w || ah != g.surf_h || scale_factor() != g.surf_scale) g.resize_fn(aw, ah, g.resize_data);
  }
  return G_SOURCE_REMOVE;
}

static void on_size_allocate (GtkWidget *w, GdkRectangle *a, gpointer data) {
  g.want_w = a->width;
  g.want_h = a->height;
  /* Redraw once things settle rather than for every intermediate size. */
  if (g.resize_fn && !g.resize_idle) g.resize_idle = g_idle_add(do_resize, NULL);
}

static void on_scale_changed (GObject *obj, GParamSpec *p, gpointer data) {
  /* Window moved between a Retina and a non-Retina screen. */
  if (g.resize_fn && !g.resize_idle) g.resize_idle = g_idle_add(do_resize, NULL);
}

static void on_destroy (GtkWidget *w, gpointer data) {
  /* The user closed the window. It is reopened on the next qpg_begin_draw. */
  g.win = NULL;
  g.area = NULL;
}

/*-------------------------------------------------------------------------*/

int qpg_init (void) {
  if (!g.initialized) {
    g.initialized = 1;
    g.no_display = !gtk_init_check(NULL, NULL);
    g.want_w = 600;
    g.want_h = 400;
    strcpy(g.title, "quick_plot");
  }
  return g.no_display ? -1 : 0;
}

/*-------------------------------------------------------------------------*/
/* Window icon: a small plot drawn with Cairo, so no image file needs to be installed. */

static cairo_surface_t *make_icon_surface (int n) {
  cairo_surface_t *surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, n, n);
  cairo_t *cr = cairo_create(surf);
  double u = n / 64.0, x;
  int i;

  cairo_scale(cr, u, u);

  /* Page: white rounded square with a gray border. */
  cairo_new_sub_path(cr);
  cairo_arc(cr, 52, 12, 8, -M_PI/2, 0);
  cairo_arc(cr, 52, 52, 8, 0, M_PI/2);
  cairo_arc(cr, 12, 52, 8, M_PI/2, M_PI);
  cairo_arc(cr, 12, 12, 8, M_PI, 3*M_PI/2);
  cairo_close_path(cr);
  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_fill_preserve(cr);
  cairo_set_source_rgb(cr, 0.45, 0.45, 0.45);
  cairo_set_line_width(cr, 2);
  cairo_stroke(cr);

  /* Axes box. */
  cairo_set_source_rgb(cr, 0, 0, 0);
  cairo_set_line_width(cr, 1.5);
  cairo_rectangle(cr, 12, 12, 40, 40);
  cairo_stroke(cr);

  /* Two curves, as in a Tao beta function plot. */
  cairo_set_line_width(cr, 2.5);
  cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  for (i = 0; i < 2; i++) {
    for (x = 0; x <= 40.01; x += 0.5) {
      double y = i == 0 ? 38 - 14 * pow(sin(M_PI * x / 20), 2) : 22 + 14 * pow(sin(M_PI * (x + 10) / 20), 2);
      if (x == 0) cairo_move_to(cr, 12 + x, y);
      else        cairo_line_to(cr, 12 + x, y);
    }
    if (i == 0) cairo_set_source_rgb(cr, 0, 0, 1);
    else        cairo_set_source_rgb(cr, 1, 0.65, 0);
    cairo_stroke(cr);
  }

  cairo_destroy(cr);
  cairo_surface_flush(surf);
  return surf;
}

#if defined (GDK_WINDOWING_QUARTZ)

/* macOS: GTK does not set the Dock icon (gtk_window_set_icon is ignored), so set it with
 * [NSApp setApplicationIconImage:]. The Objective-C runtime is looked up with dlsym (it is
 * loaded by GTK's Quartz backend) so no extra link flags are needed. */

typedef void *qpg_id;
typedef void *qpg_sel;

static cairo_status_t png_to_bytes (void *closure, const unsigned char *data, unsigned int len) {
  g_byte_array_append((GByteArray *)closure, data, len);
  return CAIRO_STATUS_SUCCESS;
}

static void set_macos_dock_icon (void) {
  qpg_id (*get_class) (const char *) = (qpg_id (*)(const char *))dlsym(RTLD_DEFAULT, "objc_getClass");
  qpg_sel (*sel) (const char *) = (qpg_sel (*)(const char *))dlsym(RTLD_DEFAULT, "sel_registerName");
  void *msg = dlsym(RTLD_DEFAULT, "objc_msgSend");
  qpg_id data, img, app;
  GByteArray *png;
  cairo_surface_t *surf;

  if (!get_class || !sel || !msg) return;
  if (!get_class("NSApplication") || !get_class("NSImage") || !get_class("NSData")) return;

  surf = make_icon_surface(256);
  png = g_byte_array_new();
  cairo_surface_write_to_png_stream(surf, png_to_bytes, png);
  cairo_surface_destroy(surf);

  data = ((qpg_id (*)(qpg_id, qpg_sel))msg)(get_class("NSData"), sel("alloc"));
  data = ((qpg_id (*)(qpg_id, qpg_sel, const void *, unsigned long))msg)(data, sel("initWithBytes:length:"),
                                                                          png->data, (unsigned long)png->len);
  img = ((qpg_id (*)(qpg_id, qpg_sel))msg)(get_class("NSImage"), sel("alloc"));
  if (data && img) img = ((qpg_id (*)(qpg_id, qpg_sel, qpg_id))msg)(img, sel("initWithData:"), data);
  app = ((qpg_id (*)(qpg_id, qpg_sel))msg)(get_class("NSApplication"), sel("sharedApplication"));
  if (app && img) ((void (*)(qpg_id, qpg_sel, qpg_id))msg)(app, sel("setApplicationIconImage:"), img);
  g_byte_array_unref(png);
}

#endif

static void set_icon (void) {
  if (!g.icon) {
    cairo_surface_t *surf = make_icon_surface(64);
    g.icon = gdk_pixbuf_get_from_surface(surf, 0, 0, 64, 64);
    cairo_surface_destroy(surf);
#if defined (GDK_WINDOWING_QUARTZ)
    if (strcmp(G_OBJECT_TYPE_NAME(gdk_display_get_default()), "GdkQuartzDisplay") == 0) set_macos_dock_icon();
#endif
  }
  if (g.icon) gtk_window_set_icon(GTK_WINDOW(g.win), g.icon);
}

/*-------------------------------------------------------------------------*/

static void create_window (void) {
  g.win = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(g.win), g.title);
  set_icon();
  gtk_window_set_default_size(GTK_WINDOW(g.win), g.want_w, g.want_h);
  gtk_window_set_focus_on_map(GTK_WINDOW(g.win), FALSE);   /* Leave keyboard focus in the terminal. */

  g.area = gtk_drawing_area_new();
  gtk_container_add(GTK_CONTAINER(g.win), g.area);
  g_signal_connect(g.area, "draw", G_CALLBACK(on_draw), NULL);
  g_signal_connect(g.area, "size-allocate", G_CALLBACK(on_size_allocate), NULL);
  g_signal_connect(g.area, "notify::scale-factor", G_CALLBACK(on_scale_changed), NULL);
  g_signal_connect(g.win, "destroy", G_CALLBACK(on_destroy), NULL);

  gtk_widget_show_all(g.win);
  qpg_pump();
}

int qpg_open (int width, int height, const char *title) {
  if (qpg_init() != 0) return -1;
  g.want_w = width;
  g.want_h = height;
  if (title) {
    strncpy(g.title, title, sizeof(g.title) - 1);
    g.title[sizeof(g.title) - 1] = 0;
  }
  if (!g.win) create_window();
  ensure_surface(width, height);
  return 0;
}

cairo_t *qpg_begin_draw (int width, int height) {
  if (qpg_init() != 0) return NULL;
  if (!g.win) create_window();
  ensure_surface(width, height);
  return cairo_create(g.surf);
}

void qpg_end_draw (cairo_t *cr) {
  if (!cr) return;
  cairo_destroy(cr);
  cairo_surface_flush(g.surf);
  if (g.area) gtk_widget_queue_draw(g.area);
}

void qpg_pump (void) {
  if (!g.initialized || g.no_display) return;
  while (gtk_events_pending()) gtk_main_iteration_do(FALSE);
}

void qpg_set_resize_callback (qpg_resize_fn fn, void *data) {
  g.resize_fn = fn;
  g.resize_data = data;
}

void qpg_close (void) {
  if (g.page_cr) {cairo_destroy(g.page_cr); g.page_cr = NULL;}
  if (g.win) gtk_widget_destroy(g.win);
  if (g.surf) cairo_surface_destroy(g.surf);
  g.surf = NULL;
  g.surf_w = g.surf_h = 0;
  qpg_pump();
}

int qpg_display_info (double *width_mm, double *height_mm, double *x_res, double *y_res,
                      int *width_px, int *height_px, int *scale) {
  GdkDisplay *disp;
  GdkMonitor *mon;
  GdkRectangle geom;
  int wmm, hmm;

  *width_mm = *height_mm = *x_res = *y_res = 0;
  *width_px = *height_px = 0;
  *scale = 1;
  if (qpg_init() != 0) return -1;

  disp = gdk_display_get_default();
  mon = gdk_display_get_primary_monitor(disp);
  if (!mon) mon = gdk_display_get_monitor(disp, 0);
  if (!mon) return -1;

  gdk_monitor_get_geometry(mon, &geom);
  *width_px = geom.width;
  *height_px = geom.height;
  wmm = gdk_monitor_get_width_mm(mon);
  hmm = gdk_monitor_get_height_mm(mon);
  *scale = gdk_monitor_get_scale_factor(mon);
  if (wmm <= 0 || hmm <= 0) return -1;

  *width_mm = wmm;
  *height_mm = hmm;
  *x_res = (double)geom.width / wmm;
  *y_res = (double)geom.height / hmm;
  return 0;
}

int qpg_write_png (const char *file) {
  if (!g.surf) return -1;
  return cairo_surface_write_to_png(g.surf, file) == CAIRO_STATUS_SUCCESS ? 0 : -1;
}

/*-------------------------------------------------------------------------*/
/* Page interface for PLplot. */

/* readline getc function: process window events every 20 ms until a key is available, then
 * read it with readline's own getc function.
 *
 * This is used instead of readline's rl_event_hook: on macOS with the conda-forge readline the
 * rl_event_hook input path was seen to stop reading keys and spin at 100% CPU (issue #2276).
 * Here readline reads every key the normal way.
 *
 * Signals that readline has caught (Ctrl-C, terminal resize) are handled while waiting so
 * that Ctrl-C works even if no key is pressed afterwards. */

static int getc_with_events (FILE *stream) {
  struct pollfd pfd;
  pfd.fd = fileno(stream);
  pfd.events = POLLIN;
  pfd.revents = 0;

  for (;;) {
    int n;
#if defined (RL_READLINE_VERSION) && RL_READLINE_VERSION >= 0x0800
    if (rl_pending_signal()) rl_check_signals();
#endif
    n = poll(&pfd, 1, 20);
    if (n > 0 || (n < 0 && errno != EINTR)) break;   /* Key available (or error): let readline read it. */
    qpg_pump();
  }

  return g.old_getc ? g.old_getc(stream) : rl_getc(stream);
}

static void tty_idle (void) {qpg_pump();}

static void install_hook (void) {
  if (!g.hook_installed && !g.no_display) {
    g.old_getc = rl_getc_function;
    rl_getc_function = getc_with_events;
    acc_tty_idle_hook = tty_idle;
    g.hook_installed = 1;
  }
}

static void remove_hook (void) {
  if (g.hook_installed) {
    rl_getc_function = g.old_getc;
    acc_tty_idle_hook = NULL;
    g.hook_installed = 0;
  }
}

cairo_t *qpg_page_open (int width, int height, const char *title) {
  if (g.page_cr) qpg_page_close();

  if (qpg_init() == 0) {
    qpg_open(width, height, title);
  } else {
    fprintf(stderr, "quick_plot: No display available: plots will not be shown on screen.\n");
    ensure_surface(width, height);
  }
  g.page_cr = cairo_create(g.surf);
  install_hook();
  return g.page_cr;
}

void qpg_page_present (void) {
  if (!g.page_cr || g.no_display) return;
  cairo_surface_flush(g.surf);
  if (!g.win) create_window();   /* The user closed the window: show it again. */
  if (g.area) gtk_widget_queue_draw(g.area);
  qpg_pump();
  /* Testing aid: save what the window shows. */
  if (getenv("QP_GTK_SNAPSHOT")) qpg_write_png(getenv("QP_GTK_SNAPSHOT"));
}

void qpg_page_close (void) {
  remove_hook();
  qpg_close();
}

/*-------------------------------------------------------------------------*/
/* Screen interface for the quick_plot Cairo backend (qp_cairo.c). */

int qpg_screen_open (int width, int height, const char *title) {
  if (qpg_init() != 0) return -1;
  qpg_open(width, height, title);
  install_hook();
  return 0;
}

void qpg_screen_close (void) {
  remove_hook();
  qpg_close();
}

#if defined (CESR_PLPLOT)

/* From plplot.h (whose headers are not installed with the Bmad PLplot build). */
extern void pl_cmd (int32_t op, void *ptr);
#define QPG_PLESC_DEVINIT 26

int qpg_plplot_open (int width, int height) {
  cairo_t *cr = qpg_page_open(width, height, "quick_plot");
  if (!cr) return -1;
  pl_cmd(QPG_PLESC_DEVINIT, cr);
  return 0;
}

#endif

#endif
