#!/usr/bin/env python3

"""
Combine Lux det_pix files into one, for when a simulation has to be broken up
into multiple runs.

  combine_det_pix_files.py run1.det_pix run2.det_pix ... {-out <file>}

The result is written in the same format Lux uses, so it can be fed to
plot_det_pix.py or read with det_pix.py like any other det_pix file.

What "combining" means here: the intensity Lux writes is normalized by the
number of photons tracked (see lux_param%intensity_normalization_coef), so the
runs are not simply added up. Each file is first put back on an unnormalized
footing, the runs are summed, and the total is renormalized by the total number
of photons tracked over all of the runs. Averages (E_ave, Ang_x_ave, ...) are
combined weighted by unnormalized intensity, which is how Lux forms them in the
first place, and the rms values are combined through the second moments.

Only pixels that are in the input files can be in the output. Lux leaves out
pixels below lux_param%intensity_min_det_pixel_cutoff, so set that low enough in
the runs being combined.
"""

import os
import sys
import argparse

import numpy as np

import det_pix

# The (average, rms) column pairs. Both members of a pair are combined together
# because the new rms depends on the old averages.

AVE_RMS = [('e_ave', 'e_rms'),
           ('ang_x_ave', 'ang_x_rms'),
           ('ang_y_ave', 'ang_y_rms'),
           ('init_x_ave', 'init_x_rms'),
           ('init_ang_x_ave', 'init_ang_x_rms'),
           ('init_y_ave', 'init_y_rms'),
           ('init_ang_y_ave', 'init_ang_y_rms')]

# Runs that disagree on any of these are not measuring the same thing.

MUST_MATCH = ['intensity_normalization_coef', 'normalization_includes_pixel_area',
              'dx_pixel', 'dy_pixel']

# Header parameters this script knows about. Anything else in a header comes
# from lux_param%bmad_parameters_out and is passed through.

KNOWN = MUST_MATCH + ['master_parameter_file', 'lattice_file', 'normalization',
                      'intensity_x_unnorm', 'intensity_x_norm', 'intensity_y_unnorm',
                      'intensity_y_norm', 'intensity_unnorm', 'intensity_norm',
                      'n_track_tot', 'n_at_detec', 'n_hit_pixel',
                      'nx_active_min', 'nx_active_max', 'ny_active_min', 'ny_active_max',
                      'x_center_det', 'y_center_det', 'x_rms_det', 'y_rms_det']

#------------------------------------------------------------------------------

def moments(ave, rms, weight):
  """Combine averages and rms values of several samples into one (ave, rms).

  Arrays may have an extra leading axis, in which case the combination is over
  the first axis and the result is an array.
  """
  ave, rms, weight = np.asarray(ave), np.asarray(rms), np.asarray(weight)
  total = weight.sum(axis = 0)
  bad = (total == 0)
  total = np.where(bad, 1.0, total)                 # Avoid 0/0 for empty pixels.

  new_ave = (ave * weight).sum(axis = 0) / total
  mean_sq = ((rms**2 + ave**2) * weight).sum(axis = 0) / total
  new_rms = np.sqrt(np.maximum(0.0, mean_sq - new_ave**2))
  return np.where(bad, 0.0, new_ave), np.where(bad, 0.0, new_rms)

#------------------------------------------------------------------------------

def check_input(dps):
  """Complain about anything that makes the files uncombinable."""

  for dp in dps:
    if not dp.is_2d:
      sys.exit('Not a det_pix file (the .x, .y and .histogram files cannot be '
               'combined): ' + dp.file_name)
    for key in MUST_MATCH + ['normalization', 'n_track_tot']:
      if key not in dp.params:
        sys.exit('Header parameter "%s" is missing from: %s' % (key, dp.file_name))

  first = dps[0]
  for dp in dps[1:]:
    for key in MUST_MATCH:
      if dp.params[key] != first.params[key]:
        sys.exit('Files disagree on %s: %s in %s but %s in %s' %
                 (key, first.params[key], first.file_name, dp.params[key], dp.file_name))

  if any(np.any(dp.col('phase_x')) or np.any(dp.col('phase_y')) for dp in dps):
    print('WARNING: the input has non-zero phases, so it is from a coherent\n'
          '         simulation. The runs are being added incoherently and the\n'
          '         phases are dropped from the output.', file = sys.stderr)

#------------------------------------------------------------------------------

def combine(dps):
  """Merge the pixel data of several files. Returns (pixels, params)."""

  # Every pixel that appears in any of the files, in the order Lux writes them.

  keys = sorted(set().union(*[set(zip(dp.col('ix').astype(int).tolist(),
                                      dp.col('iy').astype(int).tolist())) for dp in dps]))
  where = {key: n for n, key in enumerate(keys)}
  n_pix, n_file = len(keys), len(dps)

  intens_x = np.zeros(n_pix)
  intens_y = np.zeros(n_pix)
  intensity = np.zeros(n_pix)
  n_photon = np.zeros(n_pix, dtype = np.int64)
  x_pix, y_pix = np.zeros(n_pix), np.zeros(n_pix)

  # Averages are combined all at once at the end, so keep the per-file values.

  ave = np.zeros((n_file, len(AVE_RMS), n_pix))
  rms = np.zeros((n_file, len(AVE_RMS), n_pix))
  weight = np.zeros((n_file, 1, n_pix))
  filled = np.zeros(n_pix, dtype = bool)

  dx = dps[0].params['dx_pixel']
  dy = dps[0].params['dy_pixel']

  for n, dp in enumerate(dps):
    norm = dp.params['normalization']       # Undo it: the runs are different sizes.
    at = np.array([where[k] for k in zip(dp.col('ix').astype(int).tolist(),
                                         dp.col('iy').astype(int).tolist())])
    w = dp.col('intensity') / norm

    # A pixel index has to mean the same place on the detector in every file.
    off = filled[at] & ((np.abs(x_pix[at] - dp.col('x_pix')) > 0.01 * dx) |
                        (np.abs(y_pix[at] - dp.col('y_pix')) > 0.01 * dy))
    if off.any():
      sys.exit('The detector is not in the same place in every file. Pixel %s is at '
               '(%g, %g) in %s but elsewhere in another file.' %
               (keys[at[off][0]], dp.col('x_pix')[off][0], dp.col('y_pix')[off][0], dp.file_name))

    np.add.at(intens_x,  at, dp.col('intens_x') / norm)
    np.add.at(intens_y,  at, dp.col('intens_y') / norm)
    np.add.at(intensity, at, w)
    np.add.at(n_photon,  at, dp.col('n_photon').astype(np.int64))
    x_pix[at], y_pix[at] = dp.col('x_pix'), dp.col('y_pix')
    filled[at] = True
    weight[n, 0, at] = w
    for m, (a, r) in enumerate(AVE_RMS):
      ave[n, m, at], rms[n, m, at] = dp.col(a), dp.col(r)

  ave, rms = moments(ave, rms, weight)

  #--------------------------------------------------------------------
  # Header parameters. The detector-wide numbers are combined from the header
  # values of the individual runs, weighted by their total intensity.

  first = dps[0]
  n_track = sum(dp.params['n_track_tot'] for dp in dps)
  norm_tot = first.params['normalization'] * first.params['n_track_tot'] / n_track
  w_file = [dp.params.get('intensity_unnorm', 0.0) for dp in dps]

  params = {}
  for key in ('master_parameter_file', 'lattice_file'):
    names = [str(dp.params.get(key, '?')) for dp in dps]
    params[key] = names[0] if len(set(names)) == 1 else '[' + ', '.join(names) + ']'

  for key in MUST_MATCH: params[key] = first.params[key]
  params['normalization'] = norm_tot

  for key in ('intensity_x_unnorm', 'intensity_y_unnorm', 'intensity_unnorm'):
    params[key] = sum(dp.params.get(key, 0.0) for dp in dps)
    params[key.replace('unnorm', 'norm')] = params[key] * norm_tot

  for key in ('n_track_tot', 'n_at_detec', 'n_hit_pixel'):
    params[key] = sum(dp.params.get(key, 0) for dp in dps)

  for key in ('nx_active_min', 'ny_active_min'):
    params[key] = min(dp.params[key] for dp in dps if key in dp.params)
  for key in ('nx_active_max', 'ny_active_max'):
    params[key] = max(dp.params[key] for dp in dps if key in dp.params)

  for axis in 'xy':
    center, spread = moments([dp.params.get('%s_center_det' % axis, 0.0) for dp in dps],
                             [dp.params.get('%s_rms_det' % axis, 0.0) for dp in dps], w_file)
    params['%s_center_det' % axis] = float(center)
    params['%s_rms_det' % axis] = float(spread)

  # Whatever lux_param%bmad_parameters_out put in the header.

  for key in first.params:
    if key in KNOWN: continue
    values = [dp.params.get(key) for dp in dps]
    params[key] = values[0] if len(set(values)) == 1 else 'CONFLICTING VALUES'

  pixels = dict(keys = keys, x_pix = x_pix, y_pix = y_pix, intens_x = intens_x,
                intens_y = intens_y, intensity = intensity, n_photon = n_photon,
                ave = ave, rms = rms)
  return pixels, params

#------------------------------------------------------------------------------

def write(file_name, pixels, params, extra_keys):
  """Write the combined data in the format Lux itself uses."""

  norm = params['normalization']
  dx, dy = params['dx_pixel'], params['dy_pixel']

  with open(file_name, 'w') as f:
    f.write('# master_parameter_file             = "%s"\n' % params['master_parameter_file'])
    f.write('# lattice_file                      = "%s"\n' % params['lattice_file'])
    f.write('# intensity_normalization_coef      = %12.4E\n' % params['intensity_normalization_coef'])
    f.write('# normalization_includes_pixel_area = %2s\n' %
            ('T' if params['normalization_includes_pixel_area'] else 'F'))
    f.write('# normalization       = %14.6E\n' % norm)
    for key in ('intensity_x_unnorm', 'intensity_x_norm', 'intensity_y_unnorm',
                'intensity_y_norm', 'intensity_unnorm', 'intensity_norm'):
      f.write('# %-19s = %16.5E\n' % (key, params[key]))
    for key in ('n_track_tot', 'n_at_detec', 'n_hit_pixel'):
      f.write('# %-19s = %d\n' % (key, params[key]))
    f.write('# dx_pixel            = %10.6f\n' % dx)
    f.write('# dy_pixel            = %10.6f\n' % dy)
    for key in ('nx_active_min', 'nx_active_max', 'ny_active_min', 'ny_active_max'):
      f.write('# %-19s = %8d\n' % (key, params[key]))

    for key in extra_keys:                      # From lux_param%bmad_parameters_out.
      value = params[key]
      f.write('%s = %16.8E\n' % (key, value) if isinstance(value, float)
              else '%s = "%s"\n' % (key, value))

    for what in ('center', 'rms'):
      for axis, d in (('x', dx), ('y', dy)):
        key = '%s_%s_det' % (axis, what)
        f.write('# %-19s =%16.5E  # %.4g pixels\n' % (key, params[key], params[key] / d))

    f.write('#-----------------------------------------------------\n')
    f.write('%s|%s\n' % ('#'.ljust(165), 'Init'.rjust(46)))
    f.write('#   ix    iy      x_pix      y_pix   Intens_x    Phase_x   Intens_y    Phase_y'
            '  Intensity    N_photon     E_ave     E_rms  Ang_x_ave  Ang_x_rms  Ang_y_ave'
            '  Ang_y_rms|     X_ave      X_rms  Ang_x_ave  Ang_x_rms      Y_ave      Y_rms'
            '  Ang_y_ave  Ang_y_rms\n')

    for n, (ix, iy) in enumerate(pixels['keys']):
      row = [pixels['x_pix'][n], pixels['y_pix'][n],
             pixels['intens_x'][n] * norm, 0.0, pixels['intens_y'][n] * norm, 0.0,
             pixels['intensity'][n] * norm]
      f.write('%6d%6d%s%12d%s%s\n' % (ix, iy,
              ''.join(('%11.3E' % v) for v in row), pixels['n_photon'][n],
              ''.join(('%10.3f' % pixels[w][0,n]) for w in ('ave', 'rms')),
              ''.join(('%11.3E' % pixels[w][m,n]) for m in range(1, len(AVE_RMS))
                                                  for w in ('ave', 'rms'))))

#------------------------------------------------------------------------------

def main():
  p = argparse.ArgumentParser(description = 'Combine Lux det_pix files into one.')
  p.add_argument('files', nargs = '+', metavar = 'FILE', help = 'det_pix files to combine')
  p.add_argument('-out', '--out', metavar = 'FILE', default = 'combined.det_pix',
                 help = 'Output file. Default: combined.det_pix')
  args = p.parse_args()

  try:
    dps = [det_pix.read(name) for name in args.files]
  except (OSError, det_pix.DetPixError) as err:
    sys.exit('%s' % err)

  if os.path.abspath(args.out) in [os.path.abspath(f) for f in args.files]:
    sys.exit('The output file is also one of the input files: ' + args.out)

  check_input(dps)
  pixels, params = combine(dps)
  extra_keys = [key for key in params if key not in KNOWN]
  write(args.out, pixels, params, extra_keys)

  print('Combined %d file%s: %d photons tracked, %d pixels with data.' %
        (len(dps), '' if len(dps) == 1 else 's', params['n_track_tot'], len(pixels['keys'])))
  print('Wrote: ' + args.out)

#------------------------------------------------------------------------------

if __name__ == '__main__':
  main()
