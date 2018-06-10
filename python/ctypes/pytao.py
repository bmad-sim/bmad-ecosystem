import os
from ctypes import CDLL, c_char_p

#--------------------------------------

class pytao:
  # Initialization
  def __init__(self, so_lib = ''):
    if so_lib == '':
      BASE_DIR=os.environ['ACC_ROOT_DIR'] + '/production/lib/'
      if os.path.isfile(BASE_DIR + 'libtao.so'):
        self.so_lib = CDLL(BASE_DIR + 'libtao.so')
      elif os.path.isfile(BASE_DIR + 'libtao.dylib'):
        self.so_lib = CDLL(BASE_DIR + 'libtao.dylib')
      elif os.path.isfile(BASE_DIR + 'libtao.dll'):
        self.so_lib = CDLL(BASE_DIR + 'libtao.dll')
      else:
        raise ValueError ('Shared object library not found in: ' + BASE_DIR + '/production/lib')
    else:
      self.so_lib = CDLL(so_lib)

    self.so_lib.tao_c_out_io_buffer_get_line.restype = c_char_p
    self.so_lib.tao_c_out_io_buffer_reset.restype = None
 
  # Used by init and cmd routines
  def get_output(self):
    n_lines = self.so_lib.tao_c_out_io_buffer_num_lines()
    lines = [self.so_lib.tao_c_out_io_buffer_get_line(i).decode('utf-8') for i in range(1, n_lines+1)]
    self.so_lib.tao_c_out_io_buffer_reset()
    return lines
 
  # Init Tao
  def init(self, cmd, list_out = False):
    self.so_lib.tao_c_init_tao(cmd.encode('utf-8'))
    output = self.get_output()
    if list_out:
      return output
    else:
      for line in output:
        print (line)
 
  # Send a command to Tao and return or print the output
  def cmd(self, cmd, list_out = False):
    self.so_lib.tao_c_command(cmd.encode('utf-8'))
    output = self.get_output()
    if list_out:
      return output
    else:
      for line in output:
        print (line)
