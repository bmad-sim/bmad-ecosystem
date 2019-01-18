import os
import ctypes
import numpy as np

#--------------------------------------

class Tao:
  """
Class to run and interact with Tao. Requires libtao shared object. 

Setup:

import os
import sys
TAO_PYTHON_DIR=os.environ['ACC_ROOT_DIR'] + '/tao/python'
sys.path.insert(0, TAO_PYTHON_DIR)
import pytao
tao = pytao.Tao()
tao.init("command line args here...")
  """

  #---------------------------------------------

  def __init__(self, init='', so_lib = ''):
    if so_lib == '':
      BASE_DIR=os.environ['ACC_ROOT_DIR'] + '/production/lib/'
      if os.path.isfile(BASE_DIR + 'libtao.so'):
        self.so_lib = ctypes.CDLL(BASE_DIR + 'libtao.so')
      elif os.path.isfile(BASE_DIR + 'libtao.dylib'):
        self.so_lib = ctypes.CDLL(BASE_DIR + 'libtao.dylib')
      elif os.path.isfile(BASE_DIR + 'libtao.dll'):
        self.so_lib = ctypes.CDLL(BASE_DIR + 'libtao.dll')
      else:
        raise ValueError ('Shared object libtao library not found in: ' + BASE_DIR)
    else:
      self.so_lib = ctypes.CDLL(so_lib)

    self.so_lib.tao_c_out_io_buffer_get_line.restype = ctypes.c_char_p
    self.so_lib.tao_c_out_io_buffer_reset.restype = None
    
    # Attributes
    self.initialized = False
 
    try:
        self.register_cell_magic()
    except:
        pass
        #print('unable to register cell magic')
        
    if init:
        # Call init
        cmd = '-init '+init
        self.init(cmd)

  #---------------------------------------------
  # Used by init and cmd routines

  def get_output(self):
    n_lines = self.so_lib.tao_c_out_io_buffer_num_lines()
    lines = [self.so_lib.tao_c_out_io_buffer_get_line(i).decode('utf-8') for i in range(1, n_lines+1)]
    self.so_lib.tao_c_out_io_buffer_reset()
    return lines
 
  #---------------------------------------------
  # Init Tao

  def init(self, cmd):
    if not self.initialized:
        self.so_lib.tao_c_init_tao(cmd.encode('utf-8'))
        self.initialized = True
        return self.get_output()
    else:
        # Reinit
        self.cmd('reinit tao')
        self.initialized = True
        return self.get_output()
        
 
  #---------------------------------------------
  # Send a command to Tao and return the output

  def cmd(self, cmd):
    self.so_lib.tao_c_command(cmd.encode('utf-8'))
    return self.get_output()

  #---------------------------------------------
  # Get real array output. 
  # Only python commands that load the real array buffer can be used with this method.

  def cmd_real (self, cmd):
    self.so_lib.tao_c_command(cmd.encode('utf-8'))
    n = self.so_lib.tao_c_real_array_size()
    self.so_lib.tao_c_get_real_array.restype = ctypes.POINTER(ctypes.c_double * n)
    
    # Old way:
    #array = []
    #for re in self.so_lib.tao_c_get_real_array().contents: array.append(re)
    #return array

    #NumPy way:
    # This is a pointer to the scratch space. 
    array = np.ctypeslib.as_array(
        (ctypes.c_double * n).from_address(ctypes.addressof(self.so_lib.tao_c_get_real_array().contents)))
    # Return a copy
    return np.copy(array)
    
    
  #----------
  # Get integer array output. 
  # Only python commands that load the integer array buffer can be used with this method.

  def cmd_integer (self, cmd):
    self.so_lib.tao_c_command(cmd.encode('utf-8'))
    n = self.so_lib.tao_c_integer_array_size()
    self.so_lib.tao_c_get_integer_array.restype = ctypes.POINTER(ctypes.c_int * n)
    array = []
    for inte in self.so_lib.tao_c_get_integer_array().contents: array.append(inte)
    return array

  #---------------------------------------------

  def register_cell_magic(self):
    """
    Registers a cell magic in Jupyter notebooks
    Invoke by
    %%tao
    sho lat
    
    """
    from IPython.core.magic import register_cell_magic    
    @register_cell_magic
    def tao(line, cell):
        cell = cell.format(**globals())
        cmds=cell.split('\n')
        output = []
        for c in cmds:
            print('-------------------------')
            print('Tao> '+c)
            res = self.cmd(c)
            for l in res:
                 print(l)
    del tao







