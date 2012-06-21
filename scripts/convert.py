#!/usr/bin/python

import sys
import shutil
import os
import re

g = re.compile(r'(\s*)(call)( element_at_s.*?,.*?,) *(\w*)(.*)', re.IGNORECASE)

print 'Number of files: ' + str(len(sys.argv[1:]))

for arg in sys.argv[1:]:
  in_file = open(arg)
  out_file = open(r'temp.temp', 'w')
  found = False

  for line in in_file.readlines(): 

    if '%do_rad_int_calc_data' in line:
      found = True
      line = line.replace('%do_rad_int_calc_data', '%calc%rad_int_for_data')

    if '%do_rad_int_calc_plotting' in line:
      found = True
      line = line.replace('%do_rad_int_calc_plotting', '%calc%rad_int_for_plotting')

    if '%do_chrom_calc' in line:
      found = True
      line = line.replace('%do_chrom_calc', '%calc%chrom')

    if '%lattice_recalc' in line:
      found = True
      line = line.replace('%lattice_recalc', '%calc%lattice')
  
    if '%mat6_recalc_on' in line:
      found = True
      line = line.replace('%mat6_recalc_on', '%calc%mat6')

    if '%track_recalc_on' in line:
      found = True
      line = line.replace('%track_recalc_on', '%calc%track')

    if '%beam_saved_at' in line:
      found = True
      line = line.replace('%beam_saved_at', '%beam%saved_at')

    if '%current_beam' in line:
      found = True
      line = line.replace('%current_beam', '%beam%current')

    if '%beam_info' in line:
      found = True
      line = line.replace('%beam_info', '%beam')

    out_file.write(line)

  out_file.close()

  if found: 
    shutil.copy('temp.temp', arg)
    print 'File modified: ' + arg
  else:
    os.remove('temp.temp')
