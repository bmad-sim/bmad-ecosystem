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

def tao_dict_set(tao_dict, set_format, pipe):
    '''
    Runs set commands through the given pipe using the information in tao_dict and set_format
    tao_dict should be a dictionary whose values are themselves dictionaries
    set_format should be a format string that will get formatted with the keys of tao_dict.
    Ex: tao_dict = {1: {'name' : 'x', 'source' : 'lat'}, 2: {'weight': 100}}
    set_format = "set data orbit.x[{}]|"
    In this example, the following set commands will be run:
    set data orbit.x[1]|name = x
    set data orbit.x[1]|source = lat
    set data orbit.x[2]|weight = 100
    Note: Absolutely no input validation is performed
    '''
    #Turn off lattice_calc_on and plot_on
    tao_globals = pipe.cmd_in("python global")
    tao_globals = tao_globals.splitlines()
    tao_globals = tao_parameter_dict(tao_globals)
    lattice_calc_on = str(tao_globals["lattice_calc_on"].value)
    plot_on = str(tao_globals["plot_on"].value)
    # STOP lattice calculation here
    pipe.cmd_in("set global lattice_calc_on = F")
    pipe.cmd_in("set global plot_on = F")
    for key in tao_dict.keys():
        for param in tao_dict[key].keys():
            pipe.cmd_in(set_format.format(str(key)) + str(param) + ' = ' + str(tao_dict[key][param]))
    # Maybe turn lattice calc back on
    pipe.cmd_in("set global plot_on = " + plot_on)
    pipe.cmd_in("set global lattice_calc_on = " + lattice_calc_on)

def tao_set(tao_list,set_str,pipe, overide=False):
    '''
    Takes a list of tk_tao_parameters and makes a call to tao
    to set the parameters to the values input by the user
    set_str should be "set global ", "set data orbit.x[10]|",
    or whatever is appropriate
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
    # Start by unrolling any STRUCTs in tao_list
    unrolled_list = []
    for item in tao_list:
        if item.param.type == 'STRUCT':
            for ttp in item._s:
                unrolled_list.append(ttp)
                # the name to use for setting is struct_name.component_name
                unrolled_list[-1].param.name = item.param.name + '.' + ttp.param._name
        else:
            unrolled_list.append(item)
    tao_list = unrolled_list
    for item in tao_list:
        # Skip items that can't vary
        if not item.param.can_vary:
            continue
        #Type casting and validation
        if item.param.type == 'INT':
            try:
                if item.tk_var.get() == "":
                    continue
                new_val = int(item.tk_var.get())
            except ValueError:
                messagebox.showwarning("Error",item.param.name + " must be an integer ")
                new_val = item.param.value
        elif item.param.type == 'REAL':
            try:
                if item.tk_var.get() == "":
                    continue
                new_val = float(item.tk_var.get())
            except ValueError:
                messagebox.showwarning(
                        "Error",item.param.name + " must be a real number")
                new_val = item.param.value
        else:
            new_val = str(item.tk_var.get()).split(';')[0]
        #Check for any change
        if item.param.type == 'INUM':
            cond = (str(new_val) != str(item.param.value))
        else:
            cond = (new_val != item.param.value)
        if cond:
            item.param.value = new_val
            update_dict[item.param.name] = True
        else:
            update_dict[item.param.name] = overide

        #Wait till the end to set lattice_calc_on and plot_on
        if item.param.name == 'lattice_calc_on':
            lattice_calc_on = str(item.param.value)
        elif item.param.name == 'plot_on':
            plot_on = str(item.param.value)
        elif update_dict[item.param.name]:
            #print(set_str + item.param.name + " = " + str(item.param.value))
            if item.param.type == 'STR':
                #if item.param.value.strip().find(' ') != -1:
                #  if ((item.param.value[0] not in ['"',"'"])
                #            & (item.param.value[-1] != item.param.value[0])):
                #        set_val = '"' + item.param.value + '"'
                #  else:
                #        set_val = item.param.value
                #else:
                set_val = item.param.value
                msg = pipe.cmd_in(set_str + item.param.name + " = " + set_val)
            else:
                msg = pipe.cmd_in(
                        set_str + item.param.name + " = " + str(item.param.value))
            #if msg.find("ERROR") != -1:
            if msg != "":
                messagebox.showwarning(item.param.name,msg)
    #Now set lattice_calc_on and plot_on
    pipe.cmd_in("set global plot_on = " + plot_on)
    pipe.cmd_in("set global lattice_calc_on = " + lattice_calc_on)
    #Re-enable input for parameters that can vary
    #for item in tao_list:
    #  if item.param.can_vary:
    #        item.tk_wid.configure(state="normal")
