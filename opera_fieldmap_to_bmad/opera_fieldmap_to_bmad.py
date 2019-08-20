#! /usr/bin/python
 
# Parse opera field map table to bmad field map format.
# Developed by: Henry Lovelace III

import os, sys, stat, math
import glob
import numpy as np
import time

#########################################################################

def main( argv ):
  #target_dir =  os.getcwd()  
  opera_file = sys.argv[1]
  #opera_file = os.getcwd(opera_file)
  o_f = open(opera_file,'r')
  bmad_parse = os.path.join(os.getcwd(),'bmad_parse_'+opera_file)
  b_p=open(bmad_parse,'w')

#########################################################################
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
  if 'CM' in nodes[1]:increment_x = .01
  if 'MM' in nodes[1]:increment_x = .001
  
  if 'CM' in nodes[2]:increment_y = .01
  if 'MM' in nodes[2]:increment_y = .001

  if 'CM' in nodes[3]:increment_z = .01
  if 'MM' in nodes[3]:increment_z = .001
  
  #print all_nodes
  #time.sleep(5)
# Hard coded limit to array in Python is 536,870,912 on 32 bit system
  b_p.write('{ master_parameter = b_field,\n')
  b_p.write('  geometry = xyz, \n')
  b_p.write('  field_type = mixed, \n')

  initial_value=nodes[8].split()
  initial_x=initial_value[0]         
  initial_y=initial_value[1]
  initial_z=initial_value[2]
  b_p.write('r0=('+str(initial_x)+'*1E-2,  '+str(initial_y)+'*1E-2,  '+str(initial_z)+'*1E-2),\n')
  b_p.write('dr=('+str(increment_x)+', '+str(increment_y)+',  '+str(increment_z)+'), \n')
  for x in range (int(tot_nodes[2])):
   for y in range (int(tot_nodes[1])):
    for z in range (int(tot_nodes[0])):       
     B_fields=nodes[i].split()
     magnet_fields_Bx = B_fields[3]
     magnet_fields_By = B_fields[4]
     magnet_fields_Bz = B_fields[5]
     
     if i==int(all_nodes+7):
      b_p.write('pt( '+str(x)+', '+str(y)+', '+str(z)+') = ( 0, 0, 0, '+str(magnet_fields_Bx)+'*1E-4,  '+str(magnet_fields_By)+'*1E-4,  '+str(magnet_fields_Bz)+'*1E-4)}\n')  
      break         
     else:
      b_p.write('pt( '+str(x)+', '+str(y)+', '+str(z)+') = ( 0, 0, 0, '+str(magnet_fields_Bx)+'*1E-4,  '+str(magnet_fields_By)+'*1E-4,  '+str(magnet_fields_Bz)+'*1E-4),\n')
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

