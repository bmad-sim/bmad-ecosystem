/*
 * qp_cairo_viewer
 *
 * On-screen window for the quick_plot Cairo backend.
 *
 * Started by quick_plot (see qp_cairo.c). Reads display-list messages (qpr_msg_header + op stream)
 * on stdin and draws them with the shared renderer (qp_cairo_render.c) so the screen matches
 * PDF/PS/SVG/PNG output exactly. The plot is scaled to fit the window, keeping its aspect ratio.
 *
 * A reader thread does the blocking reads; the main thread only handles window events, so the
 * window stays responsive (resizing, exposing) while the program that started it waits for input.
 *
 * Usage: qp_cairo_viewer [--width <points>] [--height <points>] [--title <title>]
 *
 * Window system access is through SDL2. Only this file depends on SDL.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <SDL.h>
#include <cairo.h>
#include "qp_cairo_render.h"

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

typedef struct {
  SDL_mutex *lock;
  uint8_t *list;          /* Current display list. */
  size_t len, cap;
  Uint32 ev_update;       /* SDL user event type: display list changed. */
} shared_t;

static shared_t sh;

/*-------------------------------------------------------------------------*/

static int read_full (void *p, size_t n) {
  return fread(p, 1, n, stdin) == n ? 0 : -1;
}

static void push_event (Uint32 type) {
  SDL_Event e;
  SDL_zero(e);
  e.type = type;
  SDL_PushEvent(&e);
}

static int reader_thread (void *arg) {
  qpr_msg_header h;
  uint8_t *tmp = NULL;
  size_t tmpcap = 0;
  (void)arg;

  for (;;) {
    if (read_full(&h, sizeof h) != 0) break;
    if (h.magic != QPR_MAGIC || h.version != QPR_VERSION) {
      fprintf(stderr, "qp_cairo_viewer: Bad message header. Exiting.\n");
      break;
    }
    if (h.type == QPR_MSG_CLOSE) break;

    if (h.nbytes > tmpcap) {
      tmpcap = (size_t)h.nbytes;
      tmp = (uint8_t *)realloc(tmp, tmpcap);
      if (!tmp) break;
    }
    if (h.nbytes > 0 && read_full(tmp, (size_t)h.nbytes) != 0) break;

    SDL_LockMutex(sh.lock);
    if (h.type == QPR_MSG_RESET) sh.len = 0;
    if (sh.len + h.nbytes > sh.cap) {
      size_t nc = sh.cap ? sh.cap : 65536;
      while (nc < sh.len + h.nbytes) nc *= 2;
      sh.list = (uint8_t *)realloc(sh.list, nc);
      sh.cap = nc;
    }
    memcpy(sh.list + sh.len, tmp, (size_t)h.nbytes);
    sh.len += (size_t)h.nbytes;
    SDL_UnlockMutex(sh.lock);

    push_event(sh.ev_update);
  }

  free(tmp);
  push_event(SDL_QUIT);   /* Input closed: the plotting program has finished with this window. */
  return 0;
}

/*-------------------------------------------------------------------------*/

typedef struct {
  SDL_Window *win;
  SDL_Renderer *ren;
  SDL_Texture *tex;
  int tw, th;
  cairo_surface_t *surf;
} view_t;

static void redraw (view_t *v) {
  int pw, ph, ww, wh;
  double page_w = 8.5, page_h = 11, scale, ox, oy;
  cairo_t *cr;
  qpr_renderer r;

  SDL_GetRendererOutputSize(v->ren, &pw, &ph);
  SDL_GetWindowSize(v->win, &ww, &wh);
  if (pw <= 0 || ph <= 0) return;

  if (!v->tex || pw != v->tw || ph != v->th) {
    if (v->tex) SDL_DestroyTexture(v->tex);
    if (v->surf) cairo_surface_destroy(v->surf);
    /* Cairo ARGB32 is native-endian 32 bit ARGB which is what SDL_PIXELFORMAT_ARGB8888 means. */
    v->tex = SDL_CreateTexture(v->ren, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, pw, ph);
    v->surf = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, pw, ph);
    v->tw = pw; v->th = ph;
  }

  cr = cairo_create(v->surf);
  cairo_set_source_rgb(cr, 0.85, 0.85, 0.85);   /* Area outside the page. */
  cairo_paint(cr);

  SDL_LockMutex(sh.lock);
  qpr_peek_page(sh.list, sh.len, &page_w, &page_h);
  scale = pw / page_w;
  if (ph / page_h < scale) scale = ph / page_h;
  ox = 0.5 * (pw - scale * page_w);
  oy = 0.5 * (ph - scale * page_h);

  cairo_set_source_rgb(cr, 1, 1, 1);
  cairo_rectangle(cr, ox, oy, scale * page_w, scale * page_h);
  cairo_fill(cr);

  qpr_init(&r, cr, scale, ox, oy, 1);
  r.min_line_width = (ww > 0) ? (double)pw / ww : 1;   /* One screen pixel (HiDPI aware). */
  qpr_exec(&r, sh.list, sh.len);
  SDL_UnlockMutex(sh.lock);

  cairo_destroy(cr);
  cairo_surface_flush(v->surf);

  /* Debug/testing aid: save what the window shows. */
  if (getenv("QP_CAIRO_VIEWER_SNAPSHOT")) cairo_surface_write_to_png(v->surf, getenv("QP_CAIRO_VIEWER_SNAPSHOT"));

  SDL_UpdateTexture(v->tex, NULL, cairo_image_surface_get_data(v->surf), cairo_image_surface_get_stride(v->surf));
  SDL_RenderClear(v->ren);
  SDL_RenderCopy(v->ren, v->tex, NULL, NULL);
  SDL_RenderPresent(v->ren);
}

/*-------------------------------------------------------------------------*/

int main (int argc, char **argv) {
  int width = 600, height = 800, i;
  const char *title = "quick_plot";
  view_t v;
  SDL_Thread *th;

  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--width") && i+1 < argc) width = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--height") && i+1 < argc) height = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--title") && i+1 < argc) title = argv[++i];
  }
  if (width < 50) width = 50;
  if (height < 50) height = 50;

#ifdef _WIN32
  _setmode(_fileno(stdin), _O_BINARY);
#else
  /* The viewer shares the terminal's process group with the program that started it, so a
   * Ctrl-C meant to interrupt that program (EG a Tao optimization) would also reach the viewer. */
  signal(SIGINT, SIG_IGN);
#endif

  /* SDL's own SIGINT/SIGTERM handlers would undo the above and keep SIGTERM from ending the viewer. */
  SDL_SetHint(SDL_HINT_NO_SIGNAL_HANDLERS, "1");

#ifdef SDL_HINT_WINDOW_NO_ACTIVATION_WHEN_SHOWN
  /* Do not take keyboard focus away from the terminal the user is typing commands into. */
  SDL_SetHint(SDL_HINT_WINDOW_NO_ACTIVATION_WHEN_SHOWN, "1");
#endif
  SDL_SetHint(SDL_HINT_VIDEO_ALLOW_SCREENSAVER, "1");

  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    fprintf(stderr, "qp_cairo_viewer: SDL_Init failed: %s\n", SDL_GetError());
    return 1;
  }

  memset(&v, 0, sizeof v);
  v.win = SDL_CreateWindow(title, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, width, height,
                           SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
  if (!v.win) {
    fprintf(stderr, "qp_cairo_viewer: Cannot create window: %s\n", SDL_GetError());
    return 1;
  }
  v.ren = SDL_CreateRenderer(v.win, -1, SDL_RENDERER_ACCELERATED);
  if (!v.ren) v.ren = SDL_CreateRenderer(v.win, -1, SDL_RENDERER_SOFTWARE);
  if (!v.ren) {
    fprintf(stderr, "qp_cairo_viewer: Cannot create renderer: %s\n", SDL_GetError());
    return 1;
  }

  sh.lock = SDL_CreateMutex();
  sh.ev_update = SDL_RegisterEvents(1);
  th = SDL_CreateThread(reader_thread, "qp_reader", NULL);
  SDL_DetachThread(th);

  redraw(&v);

  for (;;) {
    SDL_Event e;
    int need_redraw = 0, quit = 0;

    if (!SDL_WaitEvent(&e)) break;
    do {   /* Coalesce queued events into a single redraw. */
      if (e.type == SDL_QUIT) quit = 1;
      else if (e.type == sh.ev_update) need_redraw = 1;
      else if (e.type == SDL_WINDOWEVENT &&
               (e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED || e.window.event == SDL_WINDOWEVENT_EXPOSED))
        need_redraw = 1;
    } while (SDL_PollEvent(&e));

    if (quit) break;
    if (need_redraw) redraw(&v);
  }

  /* Do not wait for the reader thread: it may be blocked reading stdin. */
  if (v.surf) cairo_surface_destroy(v.surf);
  SDL_Quit();
  _Exit(0);
}
