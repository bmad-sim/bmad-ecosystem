from collections import OrderedDict
import string

startup_list = [
  'beam_file;FILE;T;',
  'beam_all_file;FILE;T;',
  'beam_init_position_file;FILE;T;',
  'building_wall_file;FILE;T;',
  'data_file;FILE;T;',
  'hook_init_file;FILE;T;',
  'init_file;FILE;T;',
  'lattice_file;FILE;T;',
  'plot_file;FILE;T;',
  'startup_file;FILE;T;',
  'var_file;FILE;T;',
  'slice_lattice;FILE;T;',
  'disable_smooth_line_calc;LOGIC;T;F',
  'log_startup;LOGIC;T;F',
  'no_stopping;LOGIC;T;F',
  'rf_on;LOGIC;T;F',
]

#-------------------------------------------------

class tao_parameter():

  def __init__(self, param_name, param_type, can_vary, param_value):
    self.name = param_name
    self.type = param_type
    self.can_vary = (can_vary == 'T')

    if param_type == 'STR':
      self.value = param_value
    elif param_type == 'FILE':
      self.value = param_value
    elif param_type == 'DAT_TYPE':
      self.value = param_value
    elif param_type == 'INT':
      try:
        self.value = int(param_value)
      except:
        self.value = None
    elif param_type == 'REAL':
      try:
        self.value = float(param_value)
      except:
        self.value = None
    elif param_type == 'LOGIC':
      self.value = (param_value == 'T')
    elif param_type == 'ENUM':
      self.value = param_value
    elif param_type == 'INUM':
      self.value = int(param_value)
    else:
      print ('UNKNOWN PARAMETER TYPE: ' + param_type)

  def __str__(self):
    return str(self.value)

  def __repr__(self):
    return self.type + ';' + str(self.can_vary) + ';' + str(self.value)

# An item in the parameter list is a string that looks like:
#        'lattice_file;STR;T;bmad.lat'

def tao_parameter_dict(param_list):
    this_dict = OrderedDict()
    for param in param_list:
      v = param.split(';')
      this_dict[v[0]] = tao_parameter(v[0], v[1], v[2], v[3])
    return this_dict

def str_to_tao_param(param_str):
    '''
    Takes a parameter string
    ('lattice_file;STR;T;bmad.lat')
    and returns a tao_parameter
    '''
    v = param_str.split(';')
    return tao_parameter(v[0],v[1],v[2],v[3])

#-------------------------------------------------
param_dict = tao_parameter_dict(startup_list)
