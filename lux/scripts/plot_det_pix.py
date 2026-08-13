#!/usr/bin/env python3

"""
Interactive viewer for the Lux det_pix_out_file and its ".x", ".y" and
".histogram" companions.

  plot_det_pix.py spherical.det_pix

Once the window is up, everything is changeable without rerunning the script:
pick any column in the file from the list on the left, switch between linear
and log color scales, cycle colormaps, and click on the image to cut horizontal
and vertical lineouts through that point. Zoom and pan with the usual matplotlib
toolbar.

Run with -h for the command line arguments, and press "?" in the window for the
key bindings.

To work with the data yourself, in a script or in an interactive session, use
the reader module directly:

  from det_pix import read
  dp = read('spherical.det_pix')
  x, y = dp.curve('intensity')           # For a .x or .y file.
  grid = dp.image('e_ave', scale = 1e3)  # For the det_pix file itself.
"""

import os
import sys
import glob
import argparse
import warnings

import numpy as np
import matplotlib
from matplotlib import ticker
from matplotlib.colors import Normalize, LogNorm
from matplotlib.widgets import RadioButtons, CheckButtons, Button

import det_pix

plt = None      # matplotlib.pyplot, imported by main() after the backend is set.

#------------------------------------------------------------------------------
# Colormaps. Sequential (single ramp, dark to light) for magnitudes, diverging
# (two hues about a neutral middle) for quantities that change sign, and cyclic
# for phases. First in each list is the default.

CMAPS = {
  'sequential': ['inferno', 'viridis', 'magma', 'cividis', 'Greys_r', 'gnuplot2'],
  'diverging':  ['RdBu_r', 'coolwarm', 'PuOr_r'],
  'cyclic':     ['twilight', 'twilight_shifted', 'hsv'],
}

EMPTY_COLOR = '#9a9a9a'    # Pixels with no data. Neither end of any colormap.
LINE_COLOR  = '#1a6fdf'
CUT_COLOR   = '#e8590c'

UNITS = {1: 'm', 1e2: 'cm', 1e3: 'mm', 1e6: r'$\mu$m', 1e9: 'nm'}

# Figure layout, in fractions of the figure. The image gets IMG_RECT; the
# lineout panels and the colorbar are hung off whatever box it ends up with.

IMG_RECT = [0.28, 0.10, 0.42, 0.60]
PAD, TOP_H, RIGHT_W, CBAR_W, CBAR_PAD = 0.012, 0.15, 0.10, 0.017, 0.045

#------------------------------------------------------------------------------

def color_kind(name, z):
  """Which family of colormap suits column <name> with values <z>."""
  if name.startswith('phase'): return 'cyclic'
  z = z[np.isfinite(z)]
  if len(z) and z.min() < 0 < z.max(): return 'diverging'
  return 'sequential'

#------------------------------------------------------------------------------

def make_norm(z, kind, log):
  """Color scale limits for values <z>."""
  z = z[np.isfinite(z)]
  if len(z) == 0: return Normalize(0, 1)

  if kind == 'cyclic':
    return Normalize(-np.pi, np.pi)

  if kind == 'diverging':
    # Symmetric about zero so that the neutral middle color means zero. The
    # percentile keeps a single wild pixel from flattening everything else.
    a = np.percentile(np.abs(z), 99.5)
    if a == 0: a = max(np.abs(z).max(), 1e-30)
    return Normalize(-a, a)

  if log:
    pos = z[z > 0]
    if len(pos): return LogNorm(max(pos.min(), z.max()/1e6), z.max())

  lo, hi = z.min(), z.max()
  return Normalize(lo, hi if hi > lo else lo + 1)

#------------------------------------------------------------------------------

def step_curve(centers, values, d):
  """Polyline drawing <values> as a bar-like step curve over the pixel bins."""
  edges = np.concatenate([centers - d/2, [centers[-1] + d/2]])
  return np.repeat(edges, 2)[1:-1], np.repeat(values, 2)

#------------------------------------------------------------------------------

class Viewer:
  """The det_pix window: image, lineouts and the controls for both."""

  HELP = '''
Keys:
  n / p    Next / previous quantity        l    Linear <-> log color scale
  m        Lineout mode: sum <-> cut       c    Next colormap
  a        Equal aspect on / off           ?    This message
  Click on the image to move the cut. Toolbar keys (s = save, o = zoom,
  h = unzoom, q = quit) work as usual.
'''

  def __init__(self, dp, args):
    self.dp     = dp
    self.args   = args
    self.name   = dp.resolve(args.plot) if args.plot else self.first_quantity()
    self.names  = dp.plottable()
    self.log    = args.log
    self.equal  = True
    self.cut    = None            # (ix, iy) array indices, or None for sums.
    self.cmap_i = 0
    self.unit   = UNITS.get(args.scale, '%g m' % (1/args.scale))
    self.unit_text = self.unit.replace(r'$\mu$', 'u')     # No math for the toolbar.
    self.build()

  def first_quantity(self):
    for name in ('intensity', 'intens_x'):
      if name in self.dp.columns: return name
    return self.dp.plottable()[0]

  #----------------------------------------------------------------------------

  def build(self):
    self.fig = plt.figure(figsize = (11.5, 7.5))
    self.fig.canvas.manager.set_window_title(os.path.basename(self.dp.file_name))

    self.ax_img   = self.fig.add_axes(IMG_RECT)
    self.ax_top   = self.fig.add_axes([0, 0, 1, 1], sharex = self.ax_img)
    self.ax_right = self.fig.add_axes([0, 0, 1, 1], sharey = self.ax_img)
    self.ax_cbar  = self.fig.add_axes([0, 0, 1, 1])
    self.sync_panels()

    for ax in (self.ax_top, self.ax_right):
      ax.grid(True, lw = 0.4, color = '#cccccc', alpha = 0.7)
      ax.tick_params(labelsize = 7)
      for side in ('top', 'right'): ax.spines[side].set_visible(False)
    plt.setp(self.ax_top.get_xticklabels(), visible = False)
    plt.setp(self.ax_right.get_yticklabels(), visible = False)
    self.ax_top.yaxis.get_offset_text().set_fontsize(6)
    self.ax_right.xaxis.get_offset_text().set_fontsize(6)

    self.ax_img.set_xlabel('x (%s)' % self.unit)
    self.ax_img.set_ylabel('y (%s)' % self.unit)
    self.ax_img.tick_params(labelsize = 8)

    grid = self.grid()
    self.im = self.ax_img.imshow(grid.z, origin = 'lower', extent = grid.extent,
                                 interpolation = 'nearest')
    self.cbar = self.fig.colorbar(self.im, cax = self.ax_cbar)
    self.ax_cbar.tick_params(labelsize = 7)

    self.line_top,   = self.ax_top.plot([], [], lw = 1.1, color = LINE_COLOR)
    self.line_right, = self.ax_right.plot([], [], lw = 1.1, color = LINE_COLOR)
    self.vline = self.ax_img.axvline(np.nan, lw = 0.8, color = CUT_COLOR, alpha = 0.9)
    self.hline = self.ax_img.axhline(np.nan, lw = 0.8, color = CUT_COLOR, alpha = 0.9)

    self.cut_label = self.ax_top.text(0.5, 1.22, '', va = 'bottom', ha = 'center', fontsize = 8,
                                      color = '#444444', transform = self.ax_top.transAxes)
    self.ax_img.format_coord = self.format_coord   # Cursor readout in the toolbar.
    self.fig.suptitle('%s\n%s' % (self.dp.title(), self.run_summary()),
                      fontsize = 8.5, color = '#444444', linespacing = 1.6)
    self.build_controls()
    self.connect()
    self.refresh(reset_limits = True)

  #----------------------------------------------------------------------------

  def sync_panels(self, event = None):
    """Keep the lineouts and the colorbar glued to the image.

    Needed because matplotlib shrinks the image box, not the axes rectangle,
    when the aspect ratio is locked -- and it shrinks it by a different amount
    every time the plot is zoomed.
    """
    p = self.ax_img.get_position()
    wanted = ((self.ax_top,   [p.x0, p.y1 + PAD, p.width, TOP_H]),
              (self.ax_right, [p.x1 + PAD, p.y0, RIGHT_W, p.height]),
              (self.ax_cbar,  [p.x1 + 2*PAD + RIGHT_W + CBAR_PAD, p.y0, CBAR_W, p.height]))

    moved = False
    for ax, box in wanted:
      now = ax.get_position()
      if max(abs(a - b) for a, b in zip([now.x0, now.y0, now.width, now.height], box)) > 1e-4:
        ax.set_position(box)
        moved = True
    if moved and event is not None: self.fig.canvas.draw_idle()

  def build_controls(self):
    n = len(self.names)
    ax = self.fig.add_axes([0.012, max(0.28, 0.90 - 0.031*n), 0.20, min(0.62, 0.031*n)])
    ax.set_title('Quantity', fontsize = 8.5, loc = 'left')
    ax.set_facecolor('#f7f7f7')
    self.radio = RadioButtons(ax, self.names, active = self.names.index(self.name))
    for label in self.radio.labels: label.set_fontsize(7)
    self.radio.on_clicked(self.on_quantity)

    ax = self.fig.add_axes([0.012, 0.155, 0.20, 0.10])
    ax.set_facecolor('#f7f7f7')
    self.check = CheckButtons(ax, ['log scale', 'equal aspect', 'cut at cursor'],
                              [self.log, True, False])
    for label in self.check.labels: label.set_fontsize(7.5)
    self.check.on_clicked(self.on_check)

    self.button = Button(self.fig.add_axes([0.012, 0.10, 0.20, 0.04]), 'colormap',
                         color = '#eeeeee')
    self.button.label.set_fontsize(8)
    self.button.on_clicked(lambda event: self.cycle_cmap())

    self.fig.text(0.012, 0.055, 'n/p quantity   l log   m mode\nc colormap   a aspect   ? help',
                  fontsize = 7, family = 'monospace', color = '#777777')

  def connect(self):
    # Free up the keys used below from the default matplotlib bindings.
    for param, key in (('keymap.yscale', 'l'), ('keymap.pan', 'p'), ('keymap.back', 'c')):
      plt.rcParams[param] = [k for k in plt.rcParams[param] if k != key]
    self.fig.canvas.mpl_connect('draw_event', self.sync_panels)
    self.fig.canvas.mpl_connect('key_press_event', self.on_key)
    self.fig.canvas.mpl_connect('button_press_event', self.on_click)

  #----------------------------------------------------------------------------

  def grid(self, name = None):
    return self.dp.image(name or self.name, margin = self.args.margin, scale = self.args.scale)

  def refresh(self, reset_limits = False):
    """Redraw everything for the current quantity and settings."""
    g = self.grid()
    self.z, self.gx, self.gy = g.z, g.x, g.y
    self.dx, self.dy = g.dx, g.dy
    self.ix0, self.iy0 = g.ix0, g.iy0

    kind = color_kind(self.name, g.z)
    name = self.args.cmap or CMAPS[kind][self.cmap_i % len(CMAPS[kind])]
    cmap = plt.get_cmap(name).with_extremes(bad = EMPTY_COLOR)   # Needs matplotlib 3.4.

    log = self.log and kind == 'sequential'
    self.im.set_data(g.z)
    self.im.set_cmap(cmap)
    self.im.set_norm(make_norm(g.z, kind, log))
    self.im.set_extent(g.extent)
    self.cbar.update_normal(self.im)
    self.ax_cbar.set_ylabel('%s%s' % (self.name, '  (log)' if log else ''), fontsize = 8)
    self.ax_img.set_aspect(self.args.aspect if self.equal else 'auto')

    if reset_limits:
      self.ax_img.set_xlim(g.extent[0], g.extent[1])
      self.ax_img.set_ylim(g.extent[2], g.extent[3])

    self.refresh_cuts()
    self.fig.canvas.draw_idle()

  def refresh_cuts(self):
    """Redraw the two lineout panels."""
    z = self.z
    if self.cut is None:
      how = np.nansum if self.name in det_pix._EXTENSIVE else np.nanmean
      with warnings.catch_warnings():
        warnings.simplefilter('ignore')          # Blank rows: the mean of nothing.
        vx, vy = how(z, axis = 0), how(z, axis = 1)
      what = 'sum' if how is np.nansum else 'mean'
      title = '%s of %s over all pixels' % (what, self.name)
      self.vline.set_xdata([np.nan]);  self.hline.set_ydata([np.nan])
    else:
      ix, iy = self.cut
      vx, vy = z[iy,:], z[:,ix]
      title = '%s cut at x = %.4g, y = %.4g %s' % (self.name, self.gx[ix], self.gy[iy], self.unit)
      self.vline.set_xdata([self.gx[ix]]);  self.hline.set_ydata([self.gy[iy]])

    vx = np.where(np.isfinite(vx), vx, np.nan)
    vy = np.where(np.isfinite(vy), vy, np.nan)
    self.line_top.set_data(*step_curve(self.gx, vx, self.dx))
    self.line_right.set_data(*reversed(step_curve(self.gy, vy, self.dy)))
    self.cut_label.set_text(title)

    good = vx[np.isfinite(vx)]
    log = self.log and len(good) > 0 and good.min() >= 0
    self.ax_top.set_yscale('log' if log else 'linear')
    self.ax_right.set_xscale('log' if log else 'linear')

    # The lineout panels are narrow, so hold down the number of tick labels.
    for axis in (self.ax_top.yaxis, self.ax_right.xaxis):
      if log:
        axis.set_major_locator(ticker.LogLocator(numticks = 4))
        axis.set_minor_locator(ticker.NullLocator())
      else:
        axis.set_major_locator(ticker.MaxNLocator(3))
        fmt = ticker.ScalarFormatter()
        fmt.set_powerlimits((-3, 4))       # Long numbers become a common factor.
        axis.set_major_formatter(fmt)

    for ax, scaley in ((self.ax_top, True), (self.ax_right, False)):
      ax.relim()
      ax.autoscale_view(scalex = not scaley, scaley = scaley)
    self.fig.canvas.draw_idle()

  def run_summary(self):
    """Second title line: the header numbers worth having in front of you."""
    p, bits = self.dp.params, []
    for label, key, is_coord in (('at detector', 'n_at_detec', False),
                                 ('intensity', 'intensity_norm', False),
                                 ('x center', 'x_center_det', True),
                                 ('y center', 'y_center_det', True),
                                 ('x rms', 'x_rms_det', True), ('y rms', 'y_rms_det', True)):
      if key not in p: continue
      if is_coord: bits.append('%s = %.4g %s' % (label, p[key] * self.args.scale, self.unit))
      else: bits.append('%s = %.5g' % (label, p[key]))
    return '   |   '.join(bits)

  #----------------------------------------------------------------------------
  # Callbacks.

  def on_quantity(self, name):
    self.name = name
    self.refresh()

  def on_check(self, label):
    self.log, self.equal, cut_mode = self.check.get_status()
    if not cut_mode: self.cut = None
    elif self.cut is None: self.cut = (len(self.gx)//2, len(self.gy)//2)
    self.refresh()

  def pixel_at(self, x, y):
    """Array indices (ix, iy) of the pixel at plot position (x, y), or None."""
    if x is None or y is None: return None
    ix = int(round((x - self.gx[0]) / self.dx))
    iy = int(round((y - self.gy[0]) / self.dy))
    if 0 <= ix < len(self.gx) and 0 <= iy < len(self.gy): return (ix, iy)
    return None

  def on_click(self, event):
    if event.inaxes is not self.ax_img or event.button != 1: return
    if self.fig.canvas.toolbar and self.fig.canvas.toolbar.mode: return   # Zooming.
    ij = self.pixel_at(event.xdata, event.ydata)
    if ij is None: return
    self.cut = ij
    if not self.check.get_status()[2]: self.check.set_active(2)           # Calls on_check.
    else: self.refresh_cuts()

  def format_coord(self, x, y):
    """What the toolbar shows as the cursor moves over the image."""
    ij = self.pixel_at(x, y)
    if ij is None: return 'x = %.4g   y = %.4g %s' % (x, y, self.unit_text)
    ix, iy = ij
    return 'x = %.4g   y = %.4g %s      pixel (%d, %d)      %s = %.6g' % (
           self.gx[ix], self.gy[iy], self.unit_text, ix + self.ix0, iy + self.iy0,
           self.name, self.z[iy,ix])

  def on_key(self, event):
    if event.key in ('n', 'p'):
      i = (self.names.index(self.name) + (1 if event.key == 'n' else -1)) % len(self.names)
      self.radio.set_active(i)                     # Calls on_quantity.
    elif event.key == 'l':   self.check.set_active(0)
    elif event.key == 'a':   self.check.set_active(1)
    elif event.key == 'm':   self.check.set_active(2)
    elif event.key == 'c':   self.cycle_cmap()
    elif event.key == '?':   print(self.HELP)

  def cycle_cmap(self):
    self.args.cmap = None
    self.cmap_i += 1
    self.refresh()

#------------------------------------------------------------------------------

def plot_1d(dp, args):
  """Viewer for the .x, .y and .histogram files: one curve, pick the column."""

  names = dp.plottable()
  name = dp.resolve(args.plot) if args.plot else ('intensity' if 'intensity' in names else names[0])
  unit = UNITS.get(args.scale, '%g m' % (1/args.scale))

  fig = plt.figure(figsize = (10, 6))
  fig.canvas.manager.set_window_title(os.path.basename(dp.file_name))
  fig.suptitle(dp.title(), fontsize = 9, color = '#444444')
  ax = fig.add_axes([0.31, 0.10, 0.64, 0.80])
  ax.grid(True, lw = 0.4, color = '#cccccc', alpha = 0.7)
  for side in ('top', 'right'): ax.spines[side].set_visible(False)
  ax.set_xlabel(('%s (%s)' % (dp.kind, unit)) if dp.kind in ('x', 'y') else dp.columns[1])

  r, v = dp.curve(name, scale = args.scale)
  curve, = ax.plot(r, v, lw = 1.2, color = LINE_COLOR)
  ax.set_ylabel(name)
  if args.log: ax.set_yscale('log')

  ax_radio = fig.add_axes([0.02, max(0.16, 0.90 - 0.031*len(names)), 0.22, min(0.74, 0.031*len(names))])
  ax_radio.set_title('Quantity', fontsize = 8.5, loc = 'left')
  ax_radio.set_facecolor('#f7f7f7')
  radio = RadioButtons(ax_radio, names, active = names.index(name))
  for label in radio.labels: label.set_fontsize(7)

  ax_check = fig.add_axes([0.02, 0.09, 0.22, 0.045])
  ax_check.set_facecolor('#f7f7f7')
  check = CheckButtons(ax_check, ['log scale'], [args.log])
  for label in check.labels: label.set_fontsize(7.5)

  def on_pick(new_name):
    curve.set_ydata(dp.curve(new_name, scale = args.scale)[1])
    ax.set_ylabel(new_name)
    ax.relim();  ax.autoscale_view()
    fig.canvas.draw_idle()

  def on_log(label):
    ax.set_yscale('log' if check.get_status()[0] else 'linear')
    fig.canvas.draw_idle()

  def on_key(event):
    if event.key in ('n', 'p'):
      i = (names.index(ax.get_ylabel()) + (1 if event.key == 'n' else -1)) % len(names)
      radio.set_active(i)
    elif event.key == 'l': check.set_active(0)

  plt.rcParams['keymap.yscale'] = [k for k in plt.rcParams['keymap.yscale'] if k != 'l']
  radio.on_clicked(on_pick)
  check.on_clicked(on_log)
  fig.canvas.mpl_connect('key_press_event', on_key)
  fig._widgets = (radio, check)        # Keep the widgets alive.
  return fig

#------------------------------------------------------------------------------

def find_data_file():
  """Default input file: det.pix, or the only det_pix file in this directory."""
  if os.path.exists('det.pix'): return 'det.pix'
  hits = [f for f in glob.glob('*det_pix*') + glob.glob('*.pix')
          if not f.endswith(('.x', '.y', '.histogram'))]
  if len(hits) == 1:
    print('Using data file: ' + hits[0])
    return hits[0]
  sys.exit('No data file given and cannot guess one. Candidates here: %s' %
           (', '.join(hits) if hits else 'none'))

#------------------------------------------------------------------------------

def parse_args():
  p = argparse.ArgumentParser(description = 'Interactive plot of a Lux det_pix file.',
                              formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  p.add_argument('data_file', nargs = '?', help = 'det_pix file. Default: det.pix')
  p.add_argument('-plot', '--plot', metavar = 'COL', default = None,
                 help = 'Column to plot. Name, abbreviation, or one of the old '
                        'x (Intens_x), y (Intens_y), i (Intensity), e (E_ave). '
                        'Default: intensity')
  p.add_argument('-scale', '--scale', type = float, default = 1e3,
                 help = 'Coordinate scale. 1e3 -> mm')
  p.add_argument('-aspect', '--aspect', type = float, default = 1.0,
                 help = 'Aspect ratio of the image')
  p.add_argument('-margin', '--margin', type = int, default = 4,
                 help = 'Blank pixels around the active area')
  p.add_argument('-cmap', '--cmap', default = None,
                 help = 'Colormap. Default: chosen to suit the quantity')
  p.add_argument('-log', '--log', action = 'store_true', help = 'Start with a log color scale')
  p.add_argument('-save', '--save', metavar = 'FILE', default = None, help = 'Write the plot to FILE')
  p.add_argument('-noshow', '--noshow', action = 'store_true', help = 'Do not open a window')
  p.add_argument('-dpi', '--dpi', type = int, default = 150, help = 'Resolution of the saved plot')
  p.add_argument('-list', '--list', action = 'store_true',
                 help = 'List the file contents and exit')
  return p.parse_args()

#------------------------------------------------------------------------------

def main():
  global plt

  args = parse_args()
  file_name = args.data_file or find_data_file()

  try:
    dp = det_pix.read(file_name)
  except (OSError, det_pix.DetPixError) as err:
    sys.exit('%s' % err)

  if args.list:
    print('File:    %s  (%s data)' % (dp.file_name, dp.kind))
    print('Rows:    %d' % len(dp.data))
    print('\nHeader parameters:')
    for k, v in dp.params.items(): print('  %-35s %s' % (k, v))
    print('\nPlottable columns:')
    for c in dp.plottable():
      v = dp.col(c)
      print('  %-20s min = %11.4g  max = %11.4g' % (c, v.min(), v.max()))
    return

  if args.noshow: matplotlib.use('Agg')
  import matplotlib.pyplot as pyplot
  plt = pyplot

  try:
    view = Viewer(dp, args) if dp.is_2d else plot_1d(dp, args)
  except det_pix.DetPixError as err:
    sys.exit('%s' % err)

  fig = view.fig if isinstance(view, Viewer) else view
  if args.save:
    fig.savefig(args.save, dpi = args.dpi)
    print('Wrote: ' + args.save)
  if not args.noshow:
    print(Viewer.HELP)
    plt.show()

#------------------------------------------------------------------------------

if __name__ == '__main__':
  main()
