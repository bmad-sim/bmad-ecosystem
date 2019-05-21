from collections import OrderedDict 

startup_list = [
  'beam_file;STR;T;xx',
  'beam_all_file;STR;T;',
  'beam_init_position_file;STR;T;',
  'building_wall_file;STR;T;',
  'data_file;STR;T;',
  'hook_init_file;STR;T;',
  'init_file;STR;T;',
  'lattice_file;STR;T;',
  'plot_file;STR;T;',
  'startup_file;STR;T;',
  'var_file;STR;T;',
  'slice_lattice;STR;T;',
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
    elif param_type == 'INT':
      self.value = int(param_value)
    elif param_type == 'REAL':
      self.value = float(param_value)
    elif param_type == 'LOGIC':
      self.value = (param_value == 'T')
    elif param_type == 'ENUM':
      self.value = param_value
    else:
      print ('UNKNOWN PARAMETER TYPE: ' + param_type)

  def __str__(self):
    return str(self.value)

  def __repr__(self):
    return self.type + ';' + str(self.can_vary) + ';' + str(self.value)

# An item in the parameter list is a string that looks like:
#        'lattice_file;STR;T;bamd.lat'   

def tao_parameter_dict(param_list):
    this_dict = OrderedDict()
    for param in param_list:
      v = param.split(';')
      this_dict[v[0]] = tao_parameter(v[0], v[1], v[2], v[3])
    return this_dict

#-------------------------------------------------

class vertical_param_frame(tk.Frame):

  def __init__(self, param):
    self.label = tk.Label(self, text = param.name)
    self.label.grid(sticky="W", row=0, column=0)

    if param.type == 'STR':
      
    elif param.type == 'INT':
      self.value = int(param_value)
    elif param.type == 'REAL':
      self.value = float(param_value)
    elif param.type == 'LOGIC':
      self.value = (param_value == 'T')
    elif param.type == 'ENUM':
      self.value = param_value


#-------------------------------------------------

