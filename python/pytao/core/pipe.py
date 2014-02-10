"""
Class for piping in commands to Tao from python and for grabbing Tao output

Example:
  from pytao.core import pipe
  pipe = pipe.tao_io('-lat my_lat.bmad')   # Init with standard command-line switches
  
  output = pipe.cmd('use var *')           # Issue a command, nothing returned
  
  output = pipe.cmd_in('show top10')       # Issue a command, capturing stdout
  print(output)                            # Get the stdout of command
  
  pipe.cmd('python help')                  # Issue a Tao python command. 
                                           # Data is stored internally as a list of strings
  data=pipe.data()                         # Retrieve this data
  print(data)
"""

from ctypes import CDLL, c_char_p
from pytao.core import stdout_redirect

class tao_io:

  #-----------------------------------------------------------------
  #
  
  def __init__(self, initargs, taolib='/home/cem52/nfs/linux_lib/production/lib/libtao.so'):
    self.initargs = initargs
    if (self.initargs == 'example'):
      self.initargs =  '-init /home/cem52/nfs/linux_lib/tao/example/tao.init -noplot'
    try: 
      self.pipe = CDLL(taolib)
      self.pipe.tao_c_set_init_args(self.initargs.encode('utf-8'))
      self.pipe.tao_c_scratch_line.restype = c_char_p
      self.is_open = True
    except:
      print('Could not load shared library: '+taolib)
      self.is_open = False

  #-----------------------------------------------------------------
  # tao_io.cmd_in method sends a command to Tao and returns the output string

  def cmd_in (self, cmd_str):
    self.cmd_str = cmd_str
    if not self.is_open: 
      print ('Not connected to Tao...')
      return None
    stdout_redirect.set()
    self.pipe.tao_c_command(cmd_str.encode('utf-8'))
    stdout_redirect.unset()
    return stdout_redirect.read()

  #-----------------------------------------------------------------
  # tao_io.cmd('tao command') sends a command to tao without returning the output. 

  def cmd (self, cmd_str):
    if not self.is_open:
      print ('Not connected to Tao...')
      return None
    self.pipe.tao_c_command(cmd_str.encode('utf-8'))

  #-----------------------------------------------------------------
  # tao_io.data() method gets the results from a Tao>python command from the scratch array
    
  def data(self):
    n_lines = self.pipe.tao_c_scratch_n_lines()
    if (n_lines>0):
      return [self.pipe.tao_c_scratch_line(n+1).decode('utf-8') for n in range(n_lines)]
    else:
      return None
