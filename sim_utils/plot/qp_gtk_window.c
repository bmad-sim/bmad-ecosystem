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
#include <gtk/gtk.h>
#include <readline/readline.h>
#include "qp_gtk_window.h"

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
  rl_hook_func_t *old_hook;   /* readline event hook in place before ours. */
  int hook_installed;
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

static void create_window (void) {
  g.win = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(g.win), g.title);
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
  qpg_pump();
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

/* Called by readline about 50 times a second while it waits for a key. */
static int readline_hook (void) {
  qpg_pump();
  if (g.old_hook) return g.old_hook();
  return 0;
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

  if (!g.hook_installed && !g.no_display) {
    g.old_hook = rl_event_hook;
    rl_event_hook = readline_hook;
    rl_set_keyboard_input_timeout(20000);   /* microseconds */
    g.hook_installed = 1;
  }
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
  if (g.hook_installed) {
    rl_event_hook = g.old_hook;
    g.hook_installed = 0;
  }
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
