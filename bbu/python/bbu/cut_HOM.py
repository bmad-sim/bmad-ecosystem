#!/usr/bin/env python

from bbu import bbu_main, find_threshold
import os, sys
import glob

# From command line:  python3 cut_HOM.py 	# HOMs to cut to		'directory of original dat files'		'output cut directory'

def main( argv ):
  max_nHOMs = int(argv[0])
  print('MAX HOMS : ', max_nHOMs)
  dat_dir = argv[1]
  cut_dat_dir = argv[2]

  List = []
  Listcut = []
  files_in_dir = os.listdir(dat_dir)
#  files_in_dir = glob.glob( dat_dir.rstrip('/')+'/*dat' )
  for file_in_dir in files_in_dir:
    if (not file_in_dir.startswith('vhoms')): continue 
    if (not file_in_dir.endswith('.dat')): continue
    del List[:]
    del Listcut[:]
    f = open(os.path.join(dat_dir,file_in_dir), 'r')
    #Skip the first lines of file
    next(f)
    next(f)
    next(f)
    next(f)
    for line in f:
      s = line.split()
      if len(s) > 1:
        lr_name  = s[0]
        HOM_freq = float(s[2])
        RoQ = float(s[3])
        Q = float(s[4])
        mode = int(s[5])
        pol_angle = float(s[6])
        rating = 1/( RoQ * Q * HOM_freq)
        List.append([lr_name, HOM_freq, RoQ, Q, mode, pol_angle, rating])
    
    n = len(List) - max_nHOMs
    if (n < 0): n  = 0

    if (len(List)>0):
      Listcut = (sorted(List, key = lambda x: x[6]))
      if (len(List) > 1):
        del Listcut[-n:]
    
    if not os.path.exists(cut_dat_dir):
      os.makedirs(cut_dat_dir)
    
    # Format new data files for bmad
    f_dat = os.path.join(cut_dat_dir,'cut_'+file_in_dir)
    f = open(f_dat,'w')
    f.write('\t\t\tFrequency       R/Q                     Q       mode    Polarization_Angle\n  \
                        (Hz)            Ohm/m^(2n)                              (Radians/2pi)\n&long_range_modes\n')

    for l in range(len(Listcut)):
      name = '\tlr('+str(l+1)+') ='
      f.write('{0:10s} {1:20s} {2:20s} {3:20s} {4:10s} {5:10s}'.format(name,str(Listcut[l][1]),str(Listcut[l][2]),str(Listcut[l][3]),str(Listcut[l][4]),str(Listcut[l][5]) )+'\n')
    f.write('/')
    f.close()

if __name__ == "__main__":
  main(sys.argv[1:])
