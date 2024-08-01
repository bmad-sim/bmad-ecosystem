#!/usr/bin/env python3

# Script to combine det_pix output files from Lux into one det_pix file.
# This is useful when a simulation has to be broken up into multiple runs.

import sys
import math

T = True    # For converting from det_pix file header T/F notation
F = False

class pix_class:
  def __init__(self):
    self.x_pix          = 0
    self.y_pix          = 0
    self.intens_x       = 0
    self.phase_x        = 0
    self.intens_y       = 0
    self.phase_y        = 0
    self.intensity      = 0
    self.n_photon       = 0
    self.e_ave          = 0
    self.e_rms          = 0
    self.ang_x_ave      = 0
    self.ang_x_rms      = 0
    self.ang_y_ave      = 0
    self.ang_y_rms      = 0
    self.x_ave_init     = 0
    self.x_rms_init     = 0
    self.ang_x_ave_init = 0
    self.ang_x_rms_init = 0
    self.y_ave_init     = 0
    self.y_rms_init     = 0
    self.ang_y_ave_init = 0
    self.ang_y_rms_init = 0
  
#----------------------------------------------

def ave(ave_old, intens_old, ave_now, intens_now):
  ave_now = float(ave_now)
  intens_now = float(intens_now)
  return (ave_old * intens_old + ave_now * intens_now) / (intens_old + intens_now)

#----------------------------------------------

def rms(rms_old, ave_old, intens_old, rms_now, ave_now, intens_now):
  rms_now = float(rms_now)
  ave_now = float(ave_now)
  intens_now = float(intens_now)

  rms2_old = rms_old**2 + ave_old**2
  rms2_now = rms_now**2 + ave_now**2
  rms2_new = (rms2_old * intens_old + rms2_now * intens_now) / (intens_old + intens_now)
  ave_new = (ave_old * intens_old + ave_now * intens_now) / (intens_old + intens_now)
  return math.sqrt(max(0.0, rms2_new - ave_new**2))

#----------------------------------------------

pix_table = {}

for arg in sys.argv[1:]:
  det_file = open(arg, 'r')

  for n_header in range(1, 1000):
    line = det_file.readline()
    if line[0:3] == '#--': break
    exec (line)              # This defines photons_tracked, etc. Look at a det_pix file to understand this.
    if arg == sys.argv[1]:   # If first time.
      z = line.split('=')
      exec (z[0].strip() + '_tot =' +  z[1])  # This initalizes photons_tracked_tot, etc.

  if arg != sys.argv[1]:
    master_parameter_file_tot  = master_parameter_file_tot + ', ' + master_parameter_file
    lattice_file_tot           = lattice_file_tot + ', ' + lattice_file
    photons_tracked_tot        = photons_tracked_tot + photons_tracked
    nx_active_min_tot = min(nx_active_min_tot, nx_active_min)
    nx_active_max_tot = max(nx_active_max_tot, nx_active_max)
    ny_active_min_tot = min(ny_active_min_tot, ny_active_min)
    ny_active_max_tot = max(ny_active_max_tot, ny_active_max)
    intensity_x_unnorm_tot = intensity_x_unnorm_tot + intensity_x_unnorm
    intensity_y_unnorm_tot = intensity_y_unnorm_tot + intensity_y_unnorm

    intensity_unnorm_old   = intensity_unnorm_tot
    intensity_unnorm_tot   = intensity_unnorm_tot + intensity_unnorm

    normilization_tot    = normalization * photons_tracked / photons_tracked_tot
    intensity_x_norm_tot = intensity_x_unnorm_tot * normalization_tot
    intensity_y_norm_tot = intensity_y_norm_tot + intensity_y_norm * normalization_tot
    intensity_norm_tot   = intensity_norm_tot + intensity_norm * normalization_tot
    x_rms_det            = rms(x_rms_det_tot, x_center_det_tot, intensity_unnorm_old, x_rms_det, x_center_det, intensity_unnorm)
    y_rms_det            = rms(y_rms_det_tot, y_center_det_tot, intensity_unnorm_old, y_rms_det, y_center_det, intensity_unnorm)
    x_center_det_tot     = ave(x_center_det_tot, intensity_unnorm_old, x_center_det, intensity_unnorm)
    y_center_det_tot     = ave(y_center_det_tot, intensity_unnorm_old, y_center_det, intensity_unnorm)

  #----------------------------

  while line:
    line = det_file.readline()
    if not line: break
    if line[0] == '#': continue
    col = line.rstrip().split()
    pos = (col[0], col[1])
    if not pos in pix_table: pix_table[pos] = pix_class()
    pix = pix_table[pos]
    pix.x_pix          = col[2]
    pix.y_pix          = col[3]
    pix.intens_x       = pix.intens_x + float(col[4])
    pix.phase_x        = pix.phase_x
    pix.intens_y       = pix.intens_y + float(col[6])
    pix.phase_y        = pix.phase_y

    intensity_old = pix.intensity
    pix.intensity      = pix.intensity + float(col[8])

    pix.n_photon       = pix.n_photon + int(col[9])
    pix.e_rms          = rms(pix.e_rms, pix.e_ave, intensity_old, col[11], col[10], col[8])
    pix.e_ave          = ave (pix.e_ave, intensity_old, col[10], col[8])

    pix.ang_x_rms      = rms(pix.ang_x_rms, pix.ang_x_ave, intensity_old, col[13], col[12], col[8])
    pix.ang_x_ave      = ave(pix.ang_x_ave, intensity_old, col[12], col[8])

    pix.ang_y_rms      = rms(pix.ang_x_rms, pix.ang_x_ave, intensity_old, col[15], col[14], col[8])
    pix.ang_y_ave      = ave(pix.ang_x_ave, intensity_old, col[14], col[8])

    pix.x_rms_init     = rms(pix.x_rms_init, pix.x_ave_init, intensity_old, col[17], col[16], col[8])
    pix.x_ave_init     = ave(pix.x_ave_init, intensity_old, col[16], col[8])
    pix.ang_x_rms_init = rms(pix.ang_x_rms_init, pix.ang_x_ave_init, intensity_old, col[19], col[18], col[8])
    pix.ang_x_ave_init = ave(pix.ang_x_ave_init, intensity_old, col[18], col[8])

    pix.y_rms_init     = rms(pix.y_rms_init, pix.y_ave_init, intensity_old, col[21], col[20], col[8])
    pix.y_ave_init     = ave(pix.y_ave_init, intensity_old, col[20], col[8])
    pix.ang_y_rms_init = rms(pix.ang_y_rms_init, pix.ang_y_ave_init, intensity_old, col[23], col[22], col[8])
    pix.ang_y_ave_init = ave(pix.ang_y_ave_init, intensity_old, col[22], col[8])

  det_file.close()

#----------------------------

def to_str (str, fmt):
  fmt = '{:' + fmt + '}'
  if 'd' in fmt:
    return fmt.format(int(str))
  else:
    return fmt.format(float(str))

#----------------------------

# key is [ix_pix, iy_pix] 
def sort_order(key):
  return 100000*int(key[0]) + int(key[1])

#----------------------------

out_file = open('dp.dat', 'w')

out_file.write('master_parameter_file             = "[' + master_parameter_file_tot + ']"\n')
out_file.write('lattice_file                      = "[' + lattice_file_tot + ']"\n')
out_file.write('intensity_normalization_coef      = ' + to_str(intensity_normalization_coef_tot, '.6g') + '\n')
out_file.write('normalization_includes_pixel_area = ' + to_str(normalization_includes_pixel_area_tot, '.6g') + '\n')
out_file.write('normalization       = ' + to_str(normalization_tot, '.6g') + '\n')
out_file.write('intensity_x_unnorm  = ' + to_str(intensity_x_unnorm_tot, '.6g') + '\n')
out_file.write('intensity_x_norm    = ' + to_str(intensity_x_norm_tot, '.6g') + '\n')
out_file.write('intensity_y_unnorm  = ' + to_str(intensity_y_unnorm_tot, '.6g') + '\n')
out_file.write('intensity_y_norm    = ' + to_str(intensity_y_norm_tot, '.6g') + '\n')
out_file.write('intensity_unnorm    = ' + to_str(intensity_unnorm_tot, '.6g') + '\n')
out_file.write('intensity_norm      = ' + to_str(intensity_norm_tot, '.6g') + '\n')
out_file.write('dx_pixel            = ' + to_str(dx_pixel_tot, '.6g') + '\n')
out_file.write('dy_pixel            = ' + to_str(dy_pixel_tot, '.6g') + '\n')
out_file.write('nx_active_min       = ' + to_str(nx_active_min_tot, '.6g') + '\n')
out_file.write('nx_active_max       = ' + to_str(nx_active_max_tot, '.6g') + '\n')
out_file.write('ny_active_min       = ' + to_str(ny_active_min_tot, '.6g') + '\n')
out_file.write('ny_active_max       = ' + to_str(ny_active_max_tot, '.6g') + '\n')
out_file.write('x_center_det        = ' + to_str(x_center_det_tot, '.6g') + '\n')
out_file.write('y_center_det        = ' + to_str(y_center_det_tot, '.6g') + '\n')
out_file.write('x_rms_det           = ' + to_str(x_rms_det_tot, '.6g') + '\n')
out_file.write('y_rms_det           = ' + to_str(y_rms_det_tot, '.6g') + '\n')

out_file.write('''
#-----------------------------------------------------
#                                                                                                                                                                    |                                          Init
#   ix    iy      x_pix      y_pix   Intens_x    Phase_x   Intens_y    Phase_y  Intensity    N_photon     E_ave     E_rms  Ang_x_ave  Ang_x_rms  Ang_y_ave  Ang_y_rms|     X_ave      X_rms  Ang_x_ave  Ang_x_rms      Y_ave      Y_rms  Ang_y_ave  Ang_y_rms
''')

for k in sorted(pix_table.keys(), key = sort_order):
  p = pix_table[k]
  out_file.write(to_str(k[0], '6d') +  
                 to_str(k[1], '6d') +  
                 to_str(p.x_pix, '11.3e') + 
                 to_str(p.y_pix, '11.3e') + 
                 to_str(p.intens_x, '11.3e') + 
                 to_str(p.phase_x, '11.3e') + 
                 to_str(p.intens_y, '11.3e') + 
                 to_str(p.phase_y, '11.3e') + 
                 to_str(p.intensity, '11.3e') + 
                 to_str(p.n_photon, '12d') + 
                 to_str(p.e_ave, '10.3f') + 
                 to_str(p.e_rms, '10.3f') + 
                 to_str(p.ang_x_ave, '11.3e') + 
                 to_str(p.ang_x_rms, '11.3e') + 
                 to_str(p.ang_y_ave, '11.3e') + 
                 to_str(p.ang_y_rms, '11.3e') + 
                 to_str(p.x_ave_init, '11.3e') + 
                 to_str(p.x_rms_init, '11.3e') + 
                 to_str(p.ang_x_ave_init, '11.3e') + 
                 to_str(p.ang_x_rms_init, '11.3e') + 
                 to_str(p.y_ave_init, '11.3e') + 
                 to_str(p.y_rms_init, '11.3e') + 
                 to_str(p.ang_y_ave_init, '11.3e') + 
                 to_str(p.ang_y_rms_init, '11.3e') + '\n')

print ('Wrote: ' + out_file.name)
