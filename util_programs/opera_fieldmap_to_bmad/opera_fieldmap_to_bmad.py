#! /usr/bin/python
####PARSE OPERA FIELD MAP ####
####TABLE TO BMAD FIELD MAP FORMAT####

# Parse opera field map table to bmad field map format.
# Developed by: Henry Lovelace III

import os, sys, stat, math
import glob
import numpy as np
import time
import matplotlib.pyplot as plt
#########################################################################

def main( argv ):
  opera_file = sys.argv[1]
  o_f = open(opera_file,'r')
  bmad_parse = os.path.join(os.getcwd(),'bmad_parse_'+opera_file)
  b_p=open(bmad_parse,'w')

#########################################################################
#Gauss (1e-4)or Tesla (1)
  units = 1e-4

#ITERATOR AND POINTS
  i=8
  x=0
  y=0
  z=0
###Set increment for points###
  increment = .01
##############################
  nodes=o_f.readlines()
  tot_nodes=nodes[0].split()
  all_nodes=int(tot_nodes[0])*int(tot_nodes[1])*int(tot_nodes[2])

###Set increment for points for file###
  if 'M' in nodes[1]:increment_x = 1.0
  if 'CM' in nodes[1]:increment_x = .01
  if 'LENGU' in nodes[1]:increment_x = .01
  if 'MM' in nodes[1]:increment_x = .001
  
  if 'M' in nodes[2]:increment_y = 1.0
  if 'CM' in nodes[2]:increment_y = .01
  if 'LENGU' in nodes[2]:increment_y = .01
  if 'MM' in nodes[2]:increment_y = .001

  if 'M' in nodes[3]:increment_z = 1.0
  if 'CM' in nodes[3]:increment_z = .01
  if 'LENGU' in nodes[3]:increment_z = .01
  if 'MM' in nodes[3]:increment_z = .001
 
# Hard coded limit to array in Python is 536,870,912 on 32 bit system
  b_p.write('{ geometry = xyz, \n')
  b_p.write('  field_type = magnetic, \n')
  b_p.write('  field_scale = 1.0, \n')
  b_p.write('  ele_anchor_pt = center, \n') # double check your field map, you may want to change this.

#(x,y,z) = dr * (ix,iy,iz) + r0 + r_anchor

  initial_value=nodes[8].split()
  initial_x=float(initial_value[0])*increment_x#*1e-2         
  initial_y=float(initial_value[1])*increment_y#*1e-2
  initial_z=float(initial_value[2])*increment_z#*1e-2
  b_p.write('r0=('+str(initial_x)+', '+str(initial_y)+', '+str(initial_z)+'),\n')
#Stupid iteration I will find something smarter 
  test_x = initial_x
  test_y = initial_y
  test_z = initial_z
  
  dr_x = 0
  dr_y = 0
  dr_z = 0

  for x in range (int(tot_nodes[0])):
   for y in range (int(tot_nodes[1])):
    for z in range (int(tot_nodes[2])):       
     find_dr=nodes[i].split()
     if i==int(all_nodes+7):
      break
     else:
      dr_x = np.abs(float(find_dr[0])-test_x)
      if dr_x != 0:
       dr_x_set = dr_x
      test_x = float(find_dr[0])
      i=i+1
  i = 8

  for x in range (int(tot_nodes[0])):
   for y in range (int(tot_nodes[1])):
    for z in range (int(tot_nodes[2])):       
     find_dr=nodes[i].split()
     if i==int(all_nodes+7):
      break
     else:
      dr_y = np.abs(float(find_dr[1])-test_y)
      if dr_y != 0:
       dr_y_set = dr_y
      test_y = float(find_dr[1])
      i=i+1
  i = 8

  for x in range (int(tot_nodes[0])):
   for y in range (int(tot_nodes[1])):
    for z in range (int(tot_nodes[2])):       
     find_dr=nodes[i].split()
     if i==int(all_nodes+7):
      break
     else:
      dr_z = np.abs(float(find_dr[2])-float(test_z))
      if dr_z != 0:
       dr_z_set = dr_z
      test_z = float(find_dr[2])
      i=i+1

  i = 8
  dr_x_set = dr_x_set*increment_x
  dr_y_set = dr_y_set*increment_y
  dr_z_set = dr_z_set*increment_z
  b_p.write('dr=('+str(dr_x_set)+', '+str(dr_y_set)+', '+str(dr_z_set)+'), \n')
#end of stupid iteration
################################################################################
  for x in range (int(tot_nodes[0])):
   for y in range (int(tot_nodes[1])):
    for z in range (int(tot_nodes[2])):       
     B_fields=nodes[i].split()
     magnet_fields_Bx = float(B_fields[3])*units
     magnet_fields_By = float(B_fields[4])*units
     magnet_fields_Bz = float(B_fields[5])*units
     
     if i==int(all_nodes+7):
      b_p.write('pt( '+str(x)+', '+str(y)+', '+str(z)+') = ('+str(magnet_fields_Bx)+', '+str(magnet_fields_By)+', '+str(magnet_fields_Bz)+')}\n')  
      break         
     else:
      b_p.write('pt( '+str(x)+', '+str(y)+', '+str(z)+') = ('+str(magnet_fields_Bx)+', '+str(magnet_fields_By)+', '+str(magnet_fields_Bz)+'),\n')
      i=i+1


def append_zero(i):
 if i<10:
   return '0'+str(i)
 else:
   return str(i)


def append_two_zero(i):
  if i<10:
    return '00'+str(i)
  elif i < 100:
    return '0'+str(i)
  else:
    return str(i)



if __name__ == "__main__":
  main(sys.argv[1:])


