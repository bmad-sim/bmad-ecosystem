#!/usr/bin/env python

# Script for processing the Mathematica output into Fortran.

import sys
import shutil
import os
import re

in_file = open('exact_coefs')
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
  out_file.write(coefs[len(coefs)-1] + '_rp]), &\n')

  #

  in_file.readline()

in_file.close()
out_file.close()


