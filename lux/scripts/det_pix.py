#!/usr/bin/env python3

"""
Reader for the data files written by Lux: det_pix_out_file and the associated
".x", ".y" and ".histogram" files. All four share the same layout:

  # <param> = <value>          Header parameters.
  ...
  #---------------------       Separator.
  # ... banner ...             Zero or more banner lines.
  #   ix    iy   x_pix  ...    Column names (last comment line).
  <numbers>                    Data rows.

Nothing here plots anything -- see plot_det_pix.py for that. This module is
meant to be usable on its own:

  >>> from det_pix import read
  >>> dp = read('spherical.det_pix')
  >>> dp.params['n_track_tot']
  1000000
  >>> dp.columns
  ['ix', 'iy', 'x_pix', ..., 'init_ang_y_rms']
  >>> dp.col('intensity').sum()
  2.409e+08
  >>> g = dp.image('intensity', scale = 1e3)   # 2D grid in mm
  >>> plt.imshow(g.z, extent = g.extent, origin = 'lower')

Column names are taken from the file itself and lowercased. Names to the right
of the "|" in the header (the "Init" block) are prefixed with "init_" so that,
for example, the initial and final Ang_x_ave columns do not collide.
"""

import os
import re
import numpy as np

# Short names accepted anywhere a column name is wanted. The single letters are
# the -plot arguments used by older versions of plot_det_pix.py.

ALIASES = {
  'x': 'intens_x',
  'y': 'intens_y',
  'i': 'intensity',
  'e': 'e_ave',
  'n': 'n_photon',
  'intens': 'intensity',
  'energy': 'e_ave',
}

# Columns that add up over pixels (as opposed to averages). Empty pixels are
# zero for these and undefined (NaN) for everything else.

_EXTENSIVE = ('intensity', 'intens_x', 'intens_y', 'n_photon')

_SEPARATOR = re.compile(r'^\s*#?\s*-{3,}\s*$')
_PARAM     = re.compile(r'^\s*#?\s*([A-Za-z_]\w*)\s*=\s*(.*)$')

#------------------------------------------------------------------------------

class DetPixError(Exception):
  pass

#------------------------------------------------------------------------------

class Grid:
  """A column mapped onto the detector pixel grid.

  z          -- 2D array indexed [iy, ix], ready for imshow(origin = 'lower').
  x, y       -- 1D arrays of pixel center coordinates.
  extent     -- (x_min, x_max, y_min, y_max) of the pixel *edges*, for imshow.
  dx, dy     -- Pixel size.
  ix0, iy0   -- Pixel index of z[0,0]. So z[iy,ix] is pixel (ix+ix0, iy+iy0).
  name       -- Name of the plotted column.
  """

  def __init__(self, z, x, y, dx, dy, ix0, iy0, name):
    self.z, self.x, self.y = z, x, y
    self.dx, self.dy = dx, dy
    self.ix0, self.iy0 = ix0, iy0
    self.name = name
    self.extent = (x[0] - dx/2, x[-1] + dx/2, y[0] - dy/2, y[-1] + dy/2)

  def index_at(self, x, y):
    """Array indices (iy, ix) of the pixel containing point (x, y), or None."""
    ix = int(round((x - self.x[0]) / self.dx))
    iy = int(round((y - self.y[0]) / self.dy))
    if 0 <= ix < len(self.x) and 0 <= iy < len(self.y): return (iy, ix)
    return None

#------------------------------------------------------------------------------

class DetPix:
  """Contents of one Lux data file. Use det_pix.read() to make one."""

  def __init__(self, file_name, params, columns, data):
    self.file_name = file_name
    self.params    = params
    self.columns   = columns
    self.data      = data

    if 'ix' in columns and 'iy' in columns:
      self.kind = 'pix'         # The det_pix file itself: 2D pixel data.
      self.coords = ['ix', 'iy', 'x_pix', 'y_pix']
    elif 'ix' in columns:
      self.kind = 'x'           # det_pix.x: projection onto the x axis.
      self.coords = ['ix', 'x_pix']
    elif 'iy' in columns:
      self.kind = 'y'           # det_pix.y: projection onto the y axis.
      self.coords = ['iy', 'y_pix']
    else:
      self.kind = 'histogram'   # det_pix.histogram. Column 1 is the bin value.
      self.coords = columns[:2]

  #----------------------------------------------------------------------------

  @property
  def is_2d(self):
    return self.kind == 'pix'

  def resolve(self, name):
    """Full column name from an alias, exact name or unambiguous abbreviation."""
    key = name.strip().lower()
    key = ALIASES.get(key, key)
    if key in self.columns: return key

    hit = [c for c in self.columns if c.startswith(key)]
    if len(hit) == 1: return hit[0]
    if len(hit) > 1:
      raise DetPixError('Ambiguous column name "%s". Matches: %s' % (name, ', '.join(hit)))
    raise DetPixError('No such column: "%s". Known columns: %s' % (name, ', '.join(self.columns)))

  def col(self, name):
    """1D array of column <name>."""
    return self.data[:, self.columns.index(self.resolve(name))]

  def plottable(self):
    """Columns worth plotting -- that is, everything but the coordinates."""
    return [c for c in self.columns if c not in self.coords]

  #----------------------------------------------------------------------------

  def pixel_size(self, axis):
    """Pixel size along 'x' or 'y', from the header if it is there."""
    p = self.params.get('d%s_pixel' % axis)
    if p: return float(p)

    # Fall back on the data: <coord> = index * pixel_size + offset.
    ix, r = self.col('i' + axis), self.col('%s_pix' % axis)
    if len(ix) > 1 and ix.max() > ix.min():
      return float(np.polyfit(ix, r, 1)[0])
    raise DetPixError('Cannot determine the pixel size along %s.' % axis)

  def _axis_coords(self, axis, i_min, i_max):
    """Pixel center coordinates for pixel indices i_min ... i_max."""
    d = self.pixel_size(axis)
    ix, r = self.col('i' + axis), self.col('%s_pix' % axis)
    offset = float(np.mean(r - ix * d)) if len(ix) else 0.0   # Detector r0.
    return d, np.arange(i_min, i_max+1) * d + offset

  def active_range(self, axis):
    """Pixel index range holding the data: (min, max)."""
    lo, hi = self.params.get('n%s_active_min' % axis), self.params.get('n%s_active_max' % axis)
    idx = self.col('i' + axis)
    if len(idx) == 0:
      if lo is None: raise DetPixError('%s: no data.' % self.file_name)
      return int(lo), int(hi)
    lo = int(idx.min()) if lo is None else min(int(lo), int(idx.min()))
    hi = int(idx.max()) if hi is None else max(int(hi), int(idx.max()))
    return lo, hi

  #----------------------------------------------------------------------------

  def image(self, name, margin = 0, scale = 1.0, fill = None):
    """Column <name> put onto the 2D pixel grid. Returns a Grid.

    margin -- Number of blank pixels to pad the active area with.
    scale  -- Multiplies the coordinates. 1e3 gives millimeters.
    fill   -- Value for pixels not in the file. Defaults to 0 for intensity-like
              columns and NaN for the rest (an average over no photons is not a
              number, and NaN keeps it out of the color scale).
    """
    if not self.is_2d:
      raise DetPixError('%s is not a 2D pixel file.' % self.file_name)

    name = self.resolve(name)
    if fill is None: fill = 0.0 if name in _EXTENSIVE else np.nan

    ix_min, ix_max = self.active_range('x')
    iy_min, iy_max = self.active_range('y')
    ix_min -= margin;  ix_max += margin
    iy_min -= margin;  iy_max += margin

    z = np.full((iy_max+1-iy_min, ix_max+1-ix_min), float(fill))
    ix = np.rint(self.col('ix')).astype(int) - ix_min
    iy = np.rint(self.col('iy')).astype(int) - iy_min
    z[iy, ix] = self.col(name)

    dx, x = self._axis_coords('x', ix_min, ix_max)
    dy, y = self._axis_coords('y', iy_min, iy_max)
    return Grid(z, scale*x, scale*y, scale*dx, scale*dy, ix_min, iy_min, name)

  def curve(self, name, scale = 1.0):
    """(abscissa, value) arrays for a 1D (.x, .y or .histogram) file.

    <scale> applies to detector positions only. The abscissa of a histogram is
    whatever lux_param%histogram_variable was, so it is left alone.
    """
    if self.is_2d: raise DetPixError('%s is 2D pixel data. Use image().' % self.file_name)
    if self.kind == 'histogram': return self.data[:, 1], self.col(name)
    return scale * self.col('%s_pix' % self.kind), self.col(name)

  #----------------------------------------------------------------------------

  def title(self):
    """One line description of where the data came from."""
    bits = [os.path.basename(self.file_name)]
    if 'lattice_file' in self.params: bits.append(str(self.params['lattice_file']))
    n = self.params.get('n_track_tot')
    if n: bits.append('%g photons tracked' % n)
    return '   |   '.join(bits)

#------------------------------------------------------------------------------

def _parse_value(text):
  """Value of a header parameter: string, bool, int or float."""
  text = text.strip()
  if text[:1] in ('"', "'"):
    end = text.find(text[0], 1)
    return text[1:end] if end > 0 else text[1:]

  text = text.split('#')[0].strip()     # Lux puts a comment after some values.
  if text.upper() in ('T', '.TRUE.'):  return True
  if text.upper() in ('F', '.FALSE.'): return False
  for cast in (int, float):
    try: return cast(text)
    except ValueError: pass
  return text

#------------------------------------------------------------------------------

def _parse_columns(line, n_col):
  """Column names from the header line. Names after a "|" get an init_ prefix."""
  groups = line.lstrip('#').split('|')
  names = [c.lower() for c in groups[0].split()]
  for g in groups[1:]:
    names += ['init_' + c.lower() for c in g.split()]

  if len(names) != n_col:      # Do not guess if the header does not add up.
    names = ['col%d' % i for i in range(n_col)]
  return names

#------------------------------------------------------------------------------

def read(file_name):
  """Read a Lux det_pix file (or its .x, .y, .histogram companion)."""

  params, comments, rows = {}, [], []
  in_header = True

  with open(file_name) as f:
    for line in f:
      line = line.rstrip('\n')
      if not line.strip(): continue

      if line.lstrip().startswith('#') or in_header:
        if _SEPARATOR.match(line):
          in_header = False
          comments = []
          continue
        match = _PARAM.match(line) if in_header else None
        if match:
          params[match.group(1).lower()] = _parse_value(match.group(2))
        else:
          comments.append(line)
        continue

      rows.append(line.split())

  if not rows: raise DetPixError('No data rows found in: ' + file_name)

  n_col = max(len(r) for r in rows)
  if any(len(r) != n_col for r in rows):
    raise DetPixError('Rows with differing numbers of columns in: ' + file_name)

  try:
    data = np.array(rows, dtype = float)
  except ValueError as err:      # Fortran writes "****" for numbers that overflow.
    raise DetPixError('Cannot read the numbers in %s: %s' % (file_name, err))

  columns = _parse_columns(comments[-1] if comments else '', n_col)
  return DetPix(file_name, params, columns, data)

#------------------------------------------------------------------------------

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 2:
    sys.exit('Usage: det_pix.py <file>    # Print a summary of the file.')

  dp = read(sys.argv[1])
  print('File:    %s  (%s)' % (dp.file_name, dp.kind))
  print('Rows:    %d' % len(dp.data))
  print('\nHeader parameters:')
  for k, v in dp.params.items(): print('  %-35s %s' % (k, v))
  print('\nColumns:')
  for c in dp.columns:
    v = dp.col(c)
    print('  %-20s min = %12.5g   max = %12.5g   sum = %12.5g' % (c, v.min(), v.max(), v.sum()))
