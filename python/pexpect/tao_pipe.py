"""
Class for piping in commands to Tao from python and for grabbing Tao output

The module needs the pexpect module which can be downloaded from:
  <http://sourceforge.net/projects/pexpect/>

Example:
  import tao_pipe
  pipe = tao_pipe.tao_io("../bin/tao -lat my_lat.bmad")  # Init
  pipe.cmd("show top10")                                 # Issue a command
  tao_output = pipe.output                               # Get the output of command
"""

import pexpect

class tao_io:

  #-----------------------------------------------------------------
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

    print (self.pipe.before)

  #-----------------------------------------------------------------
  # tao_io.cmd_in method sends a command to Tao and returns the output string

  def cmd_in (self, cmd_str):

    self.cmd_str = cmd_str

    if not self.is_open: 
      print ('Not connected to Tao...')
      return ""

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
  # Note: self.pipe.after will always be the expect string which is generally "Tao>"

  def cmd (self, cmd_str):
    self.cmd_in(cmd_str)
    if self.is_open:
      print (self.pipe.after + self.pipe.before)
    else:
      print ('Not connected to Tao...')
