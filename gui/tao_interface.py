import sys
import os
sys.path.append(os.environ['ACC_ROOT_DIR']+'/tao/python/tao_pexpect')
from tao_pipe import tao_io

class tao_interface():
  '''
  Provides an interface between the GUI and
  the Tao command line
  '''
  def __init__(self, mode, init_args = "", tao_exe =  "", expect_str = "Tao>"):
    self.mode = mode
    if mode == "pexpect":
      tao_io.__init__(self, init_args, tao_exe, expect_str)
    elif mode == "ctypes":
      pass
      #Start ctypes interface

  def cmd_in(self, cmd_str):
    '''
    Runs cmd_str at the Tao command line and returns the output
    '''
    output = tao_io.cmd_in(self, cmd_str)
    # Scrub output for extra new lines at the end,
    # as well as escape characters
    output = output.replace('\r\n\r\n', '')
    output = output.replace('\x1b[6 q', '')
    return output

  def cmd(self, cmd_str):
    '''
    Runs cmd_str at the Tao command line and prints the output
    '''
    tao_io.cmd(self, cmd_str)
