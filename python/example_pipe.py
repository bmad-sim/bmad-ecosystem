#!/usr/bin/env python3

# export PYTHONPATH=$PYTHONPATH:/home/cem52/nfs/linux_lib/tao/python

from pytao.core import pipe

# Pipe init
Tao=pipe.tao_io('example')
           
#Send command to Tao and print stdout
output=Tao.cmd_in('sho ele 0')
print(output)

# Send special pyton command and view its data
Tao.cmd('python help')

print('Data from Tao command: "python help"')
print(Tao.data())

