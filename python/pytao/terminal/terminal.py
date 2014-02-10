#!/usr/bin/env python

import PyQt4.uic
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QWidget
from PyQt4.QtCore import pyqtSlot

#import pytao.core.pipe
import sys
import os

# Class definition

class TaoTerminal(QWidget):

  def __init__(self, pipe, parent = None):
    super(TaoTerminal, self).__init__(parent)

    # load the ui
    self.ui = PyQt4.uic.loadUi(os.path.dirname(os.path.abspath(__file__))+'/ui/terminal.ui', self) 
    
    # Pipe
    self.__tao_io = pipe
    
    # Single mode
    self.setFocusPolicy(QtCore.Qt.StrongFocus)
    self.singleMode = False

    # Command line history
    self.commandCounter = 0
    self.textHistory = []
    self.commandHistory = []

  # Capture key presses
  def keyPressEvent(self, e):
    print('key pressed!!!')
    if self.singleMode:
      cmd = str(e.text())+'\n'
      print('single mode command: '+cmd)
      x = str(self.pipe.cmd_in(cmd))
      self.textHistory.append('\nTao> '+cmd)
      self.textHistory.append(x)
      self.plainTextEdit.setPlainText(''.join(self.textHistory))
      # Scroll to the end
      self.plainTextEdit.moveCursor(QtGui.QTextCursor.End)
    # Up
    if e.key() == QtCore.Qt.Key_Up:
  	  if (self.commandCounter > 0):
  	    self.commandCounter -= 1
  	  self.commandLine.setText(self.commandHistory[self.commandCounter] )
    # Up
    elif e.key() == QtCore.Qt.Key_Down:
  	  if (self.commandCounter < len(self.commandHistory) -1):
  	    self.commandCounter += 1
  	  self.commandLine.setText(self.commandHistory[self.commandCounter] )
    
    
  # Function for command parsing
  @pyqtSlot()
  def on_commandLine_returnPressed(self):
    cmd = str(self.commandLine.text())
    # Save command history
    self.commandHistory.append(cmd)
    self.commandCounter = len(self.commandHistory)
    print ("Command: " + cmd)
    self.commandLine.clear()
    x = str(self.pipe.cmd_in(cmd))
    # Display
    self.textHistory.append('\nTao> '+cmd+'\n')
    self.textHistory.append(x)
    self.plainTextEdit.setPlainText(''.join(self.textHistory))
    
    # Scroll to the end
    self.plainTextEdit.moveCursor(QtGui.QTextCursor.End)
    # Print in shell as well
    print (x)

#  @pyqtSlot()
#  def on_singleModeButton_clicked(self):
#    print('on_singleModeButton_clicked!')
#    self.singleMode = not(self.singleMode)
#    print('Single Mode: ', self.singleMode)

if __name__ == "__main__":
  print('should not be here')
