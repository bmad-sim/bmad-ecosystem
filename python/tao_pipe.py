"""
Class for piping in commands to Tao from python and for grabbing Tao output

The module needs the pexpect module which can be downloaded from:
  <http://sourceforge.net/projects/pexpect/>

Usage:
  import tao_pipe
  pipe = tao_pipe.tao_io("<init_string>")
  pipe.command("<tao_cmd_line>")

"""

import pexpect

class tao_io:

  # tao_exe = Exicutable file with any command line args.

  def __init__(self, tao_exe, expect_str = 'Tao>'):
    self.pipe = pexpect.spawn (tao_exe)
    self.expect_str = expect_str
    self.is_open = True

    try:
      self.pipe.expect (expect_str)
    except pexpect.EOF:
      self.pipe.close
      self.is_open = False

    print self.pipe.before

  # command method will print the output from Tao to the terminal.

  def command (self, cmd_str):
    if not self.is_open: 
      print 'Not connected to Tao...'
      return

    self.pipe.sendline (cmd_str)
    try:
      self.pipe.expect (self.expect_str)
    except pexpect.EOF:
      print self.pipe.before
      self.pipe.close
      self.is_open = False
  
    print self.pipe.after + self.pipe.before
