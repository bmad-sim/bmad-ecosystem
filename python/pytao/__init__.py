'''
pytao is the python interface to tao.  Contains backend implementations in both
ctypes and pexpect.  The gui package supports a GUI interface to tao, in
place of the tao command line interface, with matplotlib plotting capabilities.
pytao also has some pre-defined constructs for dealing with data from tao
in the util package.
'''
from .tao_pexpect import tao_io
from .tao_ctypes import Tao, TaoModel
from .tao_interface import tao_interface
