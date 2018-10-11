import os
from ctypes import CDLL, c_char_p

#--------------------------------------

class Tao:
  """
  Class to run and interact with Tao. Requires taolib shared object. 
  
  Usage:

  tao = Tao()
  tao.init('-init tao.init') 
  
  
  """
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
 
    try:
        self.register_cell_magic()
    except:
        print('unable to register cell magic')

  # Used by init and cmd routines
  def get_output(self):
    n_lines = self.so_lib.tao_c_out_io_buffer_num_lines()
    lines = [self.so_lib.tao_c_out_io_buffer_get_line(i).decode('utf-8') for i in range(1, n_lines+1)]
    self.so_lib.tao_c_out_io_buffer_reset()
    return lines
 
  # Init Tao
  def init(self, cmd):
    self.so_lib.tao_c_init_tao(cmd.encode('utf-8'))
    return self.get_output()
 
  # Send a command to Tao and return the output
  def cmd(self, cmd):
    self.so_lib.tao_c_command(cmd.encode('utf-8'))
    return self.get_output()

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

