"""
Class for piping in commands to Tao from python and for grabbing Tao output.

This module needs the pexpect module which can be downloaded from:
  <http://sourceforge.net/projects/pexpect/>

Example:
  import tao_pipe
  pipe = tao_pipe.tao_io("-lat my_lat.bmad")   # Init
  pipe.cmd("show uni")                         # Issue a command & print output at terminal.
  tao_output = pipe.cmd_in("show uni")         # Get the output of a command. No terminal output.
"""

import pexpect
import os
import string
import sys

class tao_io:

  #-----------------------------------------------------------------
  # init_args  = Startup command line args.
  # tao_exe    = Tao executable file name including path.

  def __init__(self, init_args = '', tao_exe = '', timeout = 600, expect_str = 'Tao>'):
    if tao_exe == '': tao_exe = '$ACC_EXE/tao'
    init_string = tao_exe + ' ' + init_args
    init_string = string.Template(init_string).substitute(os.environ) # Expand environmental variables.

    if sys.version[0] == '2':                        # Python 2?
      if expect_str == 'Tao>': expect_str = u'Tao>'  # Python 2 unicode compatability.
      self.pipe = pexpect.spawn (init_string, timeout = timeout)
    else:
      self.pipe = pexpect.spawn (init_string, timeout = timeout, encoding = 'utf-8', codec_errors='replace')

    self.expect_str = expect_str
    self.is_open = True

    try:
      self.pipe.expect (expect_str)
    except pexpect.EOF:
      self.pipe.close
      self.is_open = False

    print (self.pipe.before)

  #-----------------------------------------------------------------
  # tao_io.cmd_in method sends a command to Tao and returns the output string

  def cmd_in (self, cmd_str):

    self.cmd_str = cmd_str

    if not self.is_open:
      print ('Not connected to Tao...')
      return ''

    self.pipe.sendline (cmd_str)
    try:
      self.pipe.expect (self.expect_str)
    except pexpect.EOF:
      self.pipe.close
      self.is_open = False

    # pipe.before is a multiline string. The first line is the input command string.
    # This first line is striped off before it is returned

    return self.pipe.before.partition('\n')[2].strip()

  #-----------------------------------------------------------------
  # tao_io.cmd method calls tao_io.cmd and prints the output, including
  # the command and Tao prompt, to the terminal.
  # Note: self.pipe.after will always be the expect string which is generally 'Tao>'

  def cmd (self, cmd_str):
    self.cmd_in(cmd_str)
    if self.is_open:
      print (self.pipe.after + self.pipe.before)
    else:
      print ('Not connected to Tao...')
