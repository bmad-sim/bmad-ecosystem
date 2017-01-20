#!/usr/bin/env python

# Script for processing the Mathematica output for exact bend multipoles into Fortran form.
# The output is a temporary that must be cut and pasted as appropriate.

import sys
import shutil
import os
import re

# Coefficients for evaluating F_j

in_file = open('exact_coefs', 'r')
out_file = open('coefs.out', 'w')

#for n_ord in range(2, 23):
for n_ord in range(1, 23):
  # limits
  nums = in_file.readline().strip().split(',')
  out_file.write('  exact_bend_coef_struct(' + str(n_ord) + ', ' + nums[1] + ', ' + nums[2] + ', &\n')

  # F(r) non-log coefs
  coefs = in_file.readline().strip().split(',')
  if coefs == ['']: coefs = []    # For n_ord = 1
  out_file.write('    ' + str(len(coefs)-1) + ', [')
  coefs += ['0.0'] * (12 - len(coefs))

  for ic in range(len(coefs)-1):
    if coefs[ic][-1] == '.': coefs[ic] = coefs[ic] + '0'
    if coefs[ic] == '0': coefs[ic] = '0.0'
    out_file.write(coefs[ic] + '_rp, ')
  out_file.write(coefs[len(coefs)-1] + '_rp], &\n')

  # F(r) log coefs
  coefs = in_file.readline().strip().split(',')
  if coefs == ['']: coefs = ['1.0']  # For n_ord = 1
  out_file.write('    ' + str(len(coefs)-1) + ', [')
  coefs += ['0.0'] * (11 - len(coefs))

  for ic in range(len(coefs)-1):
    if coefs[ic][-1] == '.': coefs[ic] = coefs[ic] + '0'
    if coefs[ic] == '0': coefs[ic] = '0.0'
    out_file.write(coefs[ic] + '_rp, ')
  out_file.write(coefs[len(coefs)-1] + '_rp], &\n')

  # Pade numerator coefs
  coefs = in_file.readline().strip().split(',')[n_ord:]
  out_file.write('    ' + str(len(coefs)-1) + ', [')
  coefs += ['0.0'] * (8 - len(coefs))

  for ic in range(len(coefs)-1):
    if coefs[ic][-1] == '.': coefs[ic] = coefs[ic] + '0'
    if coefs[ic] == '0': coefs[ic] = '0.0'
    out_file.write(coefs[ic] + '_rp, ')
  out_file.write(coefs[len(coefs)-1] + '_rp], &\n')

  # Pade denominator coefs
  coefs = in_file.readline().strip().split(',')
  out_file.write('    ' + str(len(coefs)-1) + ', [')
  coefs += ['0.0'] * (9 - len(coefs))

  for ic in range(len(coefs)-1):
    if coefs[ic][-1] == '.': coefs[ic] = coefs[ic] + '0'
    if coefs[ic] == '0': coefs[ic] = '0.0'
    out_file.write(coefs[ic] + '_rp, ')
  out_file.write(coefs[len(coefs)-1] + '_rp], &\n')

  # Pade 2nd derivative numerator coefs
  mx = max(n_ord-2, 0)
  coefs = in_file.readline().strip().split(',')[mx:]
  out_file.write('    ' + str(len(coefs)-1) + ', [')
  coefs += ['0.0'] * (8 - len(coefs))

  for ic in range(len(coefs)-1):
    if coefs[ic][-1] == '.': coefs[ic] = coefs[ic] + '0'
    if coefs[ic] == '0': coefs[ic] = '0.0'
    out_file.write(coefs[ic] + '_rp, ')
  out_file.write(coefs[len(coefs)-1] + '_rp], &\n')

  # Pade 2nd derivative denominator coefs
  coefs = in_file.readline().strip().split(',')
  out_file.write('    ' + str(len(coefs)-1) + ', [')
  coefs += ['0.0'] * (9 - len(coefs))

  for ic in range(len(coefs)-1):
    if coefs[ic][-1] == '.': coefs[ic] = coefs[ic] + '0'
    if coefs[ic] == '0': coefs[ic] = '0.0'
    out_file.write(coefs[ic] + '_rp, ')
  out_file.write(coefs[len(coefs)-1] + '_rp]), &\n')

  # Ignore blank line

  in_file.readline() 

in_file.close()

#-------------------------------------------------------
# Coefficients for converting between horizontally pure and vertically pure.

in_file = open('pure_conversion', 'r')

out_file.write('''
! Horizontal to Vertical Imaginary Potential

type pure_bend_multipole_struct
  real(rp) convert(0:n_pole_maxx)
end type

''')

for iz in range(4):

  if iz == 0:
    out_file.write('type (pure_bend_multipole_struct) :: v_to_h_imag(0:n_pole_maxx) = [ &\n')
  elif iz == 1:
    out_file.write('type (pure_bend_multipole_struct) :: h_to_v_imag(0:n_pole_maxx) = [ &\n')
  elif iz == 2:
    out_file.write('type (pure_bend_multipole_struct) :: v_to_h_real(0:n_pole_maxx) = [ &\n')
  elif iz == 3:
    out_file.write('type (pure_bend_multipole_struct) :: h_to_v_real(0:n_pole_maxx) = [ &\n')

  in_file.readline()   # Ignore header string

  for i in range(22):
    out_file.write('     pure_bend_multipole_struct([')
    coefs = in_file.readline().strip().split()
    for ic in range(len(coefs)):
      if coefs[ic] == '0': coefs[ic] = '0.0'
      if coefs[ic] == '1.': coefs[ic] = '1.0'
      if coefs[ic] == '-1.': coefs[ic] = '-1.0'
      out_file.write(coefs[ic] + '_rp')
      if ic == 11:
        out_file.write(', &\n        ')
      elif i == 21 and ic == 21:
        out_file.write(']) &\n]\n\n')
      elif ic == 21:
        out_file.write(']), &\n')
      else:
        out_file.write(', ')


in_file.close()

out_file.close()
