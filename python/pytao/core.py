import os
import ctypes
import numpy as np
import pytao
from .util import full_path
import tempfile
import shutil


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
        
        
        # Library needs to be set. 
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
        elif not pytao.initialized:
            self.so_lib = ctypes.CDLL(so_lib)
        else: 
            pass
            #print('Tao already initialized.')
            
        
        self.so_lib.tao_c_out_io_buffer_get_line.restype = ctypes.c_char_p
        self.so_lib.tao_c_out_io_buffer_reset.restype = None
        
        # Attributes
        ##self.initialized = False
        
        try:
            self.register_cell_magic()
        except:
            pass
            #print('unable to register cell magic')
            
        if init:
            # Call init
            self.init(init)
                    
        
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
        if not pytao.initialized:
            self.so_lib.tao_c_init_tao(cmd.encode('utf-8'))
            pytao.initialized = True
            return self.get_output()
        else:
            # Reinit
            self.cmd('reinit tao '+cmd)
            pytao.initialized = True
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
        #array = []
        #for inte in self.so_lib.tao_c_get_integer_array().contents: array.append(inte)
        #return array
        #NumPy way:
        # This is a pointer to the scratch space. 
        array = np.ctypeslib.as_array(
            (ctypes.c_int * n).from_address(ctypes.addressof(self.so_lib.tao_c_get_integer_array().contents)))
        # Return a copy
        return np.copy(array)    
    
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
    
        
    
class TaoModel(Tao):
    """
    Base class for setting up a Tao model in a directory. Builds upon the Tao class.
    
    If use_tempdir==True, then the input_file and its directory will be copied to a temporary directory.
    If workdir is given, then this temporary directory will be placed in workdir. 
    
    
    """
    def __init__(self,
                 input_file='tao.init',
                 ploton = True,
                 use_tempdir=True,
                 workdir=None,
                 verbose=True,
                 so_lib='',  # Passed onto Tao superclass
                 auto_configure = True # Should be disables if inheriting. 
                ):
        
        # Save init
        self.original_input_file = input_file
        self.ploton = ploton
        self.use_tempdir = use_tempdir
        self.workdir = workdir
        assert os.path.exists(workdir), 'workdir does not exist: '+workdir
        
        self.verbose=verbose
        self.so_lib=so_lib
    
        # Run control
        self.finished = False
        self.configured = False
        
        if os.path.exists(input_file):
            f = full_path(input_file)
            self.original_path, self.original_input_file = os.path.split(f) # Get original path, filename            
            if auto_configure:
                self.configure()
        else:
            self.vprint('Warning: Input file does not exist. Cannot configure.')
            
    def configure(self):
        
        # Set paths
        if self.use_tempdir:
            
            # Need to attach this to the object. Otherwise it will go out of scope.
            
            self.tempdir = tempfile.TemporaryDirectory(dir=self.workdir)
            # Make yet another directory to overcome the limitations of shutil.copytree
            self.path = full_path(os.path.join(self.tempdir.name, 'tao/'))
            # Copy everything in original_path
            shutil.copytree(self.original_path, self.path, symlinks=True)
        else:
            # Work in place
            self.path = self.original_path           

        self.input_file = os.path.join(self.path, self.original_input_file)                     
        
        self.vprint('Initialized Tao with '+self.input_file)                
            
            
        # Set up Tao library
        super().__init__(init=self.init_line(), so_lib=self.so_lib)                 
            
        self.configured = True

    def init_line(self):
        line = '-init '+self.input_file
        if self.ploton:
            line += ' --noplot'
        else:
            line += ' -noplot'
        return line

    def reinit(self):
        line = 'reinit tao '+self.init_line()
        self.cmd(line)
        self.vprint('Re-initialized with '+line)

    def vprint(self, *args, **kwargs):
        # Verbose print
        if self.verbose:
            print(*args, **kwargs)       
    
    def __str__(self):
        s = 'Tao Model initialized from: '+self.original_path
        s +='\n Working in path: '+self.path
        return s    