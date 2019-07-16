from parameters import tao_parameter_dict
from tkinter import messagebox

def check_for_changes(tao_list):
  '''
  Takes a list of tk_tao_parameters and returns True
  if any of the items have self.tk_var.get() !=
  self.param.value (i.e. if running tao_set would
  result in at least one set command for the list
  '''
  for item in tao_list:
    #Type casting and validation
    if item.param.type == 'INT':
      try:
        new_val = int(item.tk_var.get())
      except ValueError:
        return True
    elif item.param.type == 'REAL':
      try:
        new_val = float(item.tk_var.get())
      except ValueError:
        return True
    else:
      new_val = item.tk_var.get()
    #Check for any change
    if new_val != item.param.value:
      return True
  return False

def tao_set(tao_list,set_str,pipe, overide=False):
  '''
  Takes a list of tk_tao_parameters and makes a call to tao
  to set the parameters to the values input by the user
  set_str should be "set global ", "set data orbit.x[10]|", or whatever is appropriate
  Use the overide option to run set commands even if no change has been made.
  '''
  # Exit imediately if tao_list is empty
  if tao_list == []:
    return
  # Record the current status of global lattice_calc_on and plot_on
  tao_globals = pipe.cmd_in("python global")
  tao_globals = tao_globals.splitlines()
  tao_globals = tao_parameter_dict(tao_globals)
  lattice_calc_on = str(tao_globals["lattice_calc_on"].value)
  plot_on = str(tao_globals["plot_on"].value)
  # STOP lattice calculation here
  pipe.cmd_in("set global lattice_calc_on = F")
  pipe.cmd_in("set global plot_on = F")
  #Freeze input fields:
  #for item in tao_list:
  #  item.tk_wid.config(state="disabled")
  update_dict = {} #Record of which variables were changed
  for item in tao_list:
    #Type casting and validation
    if item.param.type == 'INT':
      try:
        new_val = int(item.tk_var.get())
      except ValueError:
        messagebox.showwarning("Error",item.param.name + " must be an integer ")
        new_val = item.param.value
    elif item.param.type == 'REAL':
      try:
        new_val = float(item.tk_var.get())
      except ValueError:
        messagebox.showwarning("Error",item.param.name + " must be a real number")
        new_val = item.param.value
    else:
      new_val = item.tk_var.get()
    #Check for any change
    if new_val != item.param.value:
      item.param.value = new_val
      update_dict[item.param.name] = True
    else:
      update_dict[item.param.name] = overide & item.param.can_vary

    #Wait till the end to set lattice_calc_on and plot_on
    if item.param.name == 'lattice_calc_on':
      lattice_calc_on = str(item.param.value)
    elif item.param.name == 'plot_on':
      plot_on = str(item.param.value)
    elif update_dict[item.param.name]:
      #print(set_str + item.param.name + " = " + str(item.param.value))
      if item.param.type == 'STR':
        if len(item.param.value) > 1:
          if ((item.param.value[0] not in ['"', '"'])
              & (item.param.value[-1] != item.param.value[0])):
            set_val = '"' + item.param.value + '"'
          else:
            set_val = item.param.value
        msg = pipe.cmd_in(set_str + item.param.name + " = " + set_val)
      else:
        msg = pipe.cmd_in(set_str + item.param.name + " = " + str(item.param.value))
      #if msg.find("ERROR") != -1:
      if msg != "":
        messagebox.showwarning(item.param.name,msg)
  #Now set lattice_calc_on and plot_on
  pipe.cmd_in("set global plot_on = " + plot_on)
  pipe.cmd_in("set global lattice_calc_on = " + lattice_calc_on)
  #Re-enable input for parameters that can vary
  #for item in tao_list:
  #  if item.param.can_vary:
  #    item.tk_wid.configure(state="normal")
