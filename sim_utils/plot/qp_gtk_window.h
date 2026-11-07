/*
 * qp_gtk_window.h
 *
 * Plot window for quick_plot built on GTK 3 and Cairo. Plain C, so it builds with gcc.
 * Compiled when CESR_GTK is defined (ACC_PLOT_DISPLAY_TYPE = "GTK").
 *
 * The plot is drawn by PLplot's extcairo driver into an off-screen Cairo image surface kept at
 * the screen's full pixel resolution (sharp on HiDPI/Retina displays). The window shows that
 * surface, so it can be repainted (uncovered, moved, resized) without redrawing the plot.
 *
 * GTK only does work when its events are processed. qpg_page_open installs a readline event
 * hook (rl_event_hook) so events are processed while the program (EG Tao) waits for input,
 * and qpg_page_present processes them after each drawing update.
 *
 * Units: GTK logical pixels ("points" on macOS).
 */

#ifndef QP_GTK_WINDOW_H
#define QP_GTK_WINDOW_H

#include <cairo.h>

/*---- Page interface used by plplot_interface.f90 (via qpg_plplot_open) ----*/

/* Open the plot window and return a Cairo context for drawing into it. If no display is
 * available (EG ssh without X, CI) the context draws off-screen so the program carries on. */
cairo_t *qpg_page_open (int width, int height, const char *title);

/* Show what has been drawn so far and process pending window events.
 * If the user closed the window it is reopened. */
void qpg_page_present (void);

/* Close the window and free the page. */
void qpg_page_close (void);

/* Open the window and attach it to the current PLplot stream (which must use the extcairo
 * device and have been initialized with plinit). Returns 0 on success. */
int qpg_plplot_open (int width, int height);

/*---- Lower level interface (also used by the standalone test program) ----*/

typedef void (*qpg_resize_fn) (int width, int height, void *data);

/* Initialize GTK. Returns 0 if a display is available, -1 if not. */
int qpg_init (void);

int qpg_open (int width, int height, const char *title);
cairo_t *qpg_begin_draw (int width, int height);
void qpg_end_draw (cairo_t *cr);

/* Process pending window events without blocking. */
void qpg_pump (void);

/* Redraw at the new size when the window is resized. Without a callback the last picture
 * is scaled to fit. */
void qpg_set_resize_callback (qpg_resize_fn fn, void *data);

void qpg_close (void);

/* Primary monitor size [mm], resolution [logical pixels / mm] and HiDPI scale factor.
 * Returns 0 on success. */
int qpg_display_info (double *width_mm, double *height_mm, double *x_res, double *y_res, int *scale);

/* Write the current off-screen picture (at device resolution) to a PNG file. */
int qpg_write_png (const char *file);

#endif
