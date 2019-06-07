import sys
import os
sys.path.append(os.environ['ACC_ROOT_DIR'] + '/tao/gui')
from parameters import *

class tao_variable():

  def __init__(self, var_name, pipe):
    self.name = var_name
    self.pipe = pipe
    param_list = pipe.cmd_in("python var1 " + self.name)
    param_list = param_list.splitlines()
    self.param_dict = tao_parameter_dict(param_list)

  def update(self, parameter, value):
    '''
    updates this variable's parameter to value
    '''
    self.pipe.cmd_in("set var " + self.name + "|" + parameter + " = " + value)
    self.param_dict[parameter].value = value
