#!/usr/bin/env python

# Script to convert phase data from Cesrv output "show phase" to Tao format.
# Input to this script should be the output of "show phase" including the two header lines but 
# excluding anything anything after the last phase data line.

# NOTE: This script is not needed for the wave analysis tutorial and is only included for historical documentation.
# NOTE: This script makes assumptions about the arrangement of BPMs in CESR that may change in the future.

import sys

in_file = open(sys.argv[1], encoding='utf-8')
out_file = open(sys.argv[1]+'.tao', 'w', encoding='utf-8')

phase_a = []
phase_b = []

for line in in_file.readlines(): 
  word = line.split()
  # Skip the header lines and the "beginning" data line.
  if word[0] == '|' or word[0] == 'BEGINNI': continue
  phase_a.append(word[3])
  phase_b.append(word[6])

# Insert zeros for det_12w2 (#18) and det_12w3 (#19) since cesrv does not return any data for them.

phase_a.insert(17, '0.0')
phase_a.insert(17, '0.0')
phase_b.insert(17, '0.0')
phase_b.insert(17, '0.0')

# And print the data in chunks of 20.

n = len(phase_a)
for i in range(0, n, 20):
  print ('set data phase.a[' + str(i+1) + ':' + str(min(i+20, n)) + ']|meas = [' + ', '.join(phase_a[i:i+20]) + ']')

for i in range(0, n, 20):
  print ('set data phase.b[' + str(i+1) + ':' + str(min(i+20, n)) + ']|meas = [' + ', '.join(phase_b[i:i+20]) + ']')
