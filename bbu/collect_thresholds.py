#!/usr/bin/env python

import os, sys

def main(argv):

  total = 0
  didnot = 0
  d = sys.argv[1]  # Directory of bbu threshold files 
  files = os.listdir(d)
  for file in files:
    if (not file.startswith('bbu_threshold')): continue
    f = open ( os.path.join(d, file),'r')
    contents = f.readlines()
    f.close
    with open(os.path.join(d,'bbu_combined_thresholds.txt'),'a') as myfile:
      for line in contents:
        total = total + 1
        if ( line.strip() == 'DID NOT CONVERGE'):
          didnot = didnot + 1
        else:
          myfile.write(line)
      myfile.close  

  print ("Of ", total, " threshold calculations, ", didnot, " did not converge")

if __name__ == "__main__":
  main(sys.argv[1:])

