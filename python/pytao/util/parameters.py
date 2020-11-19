'''
This module defines the tao_parameter class, which stores the value of a single parameter in tao.
Also defined are a few functions for parsing data into a single tao_parameter or a dictionary of tao_parameters
'''
from collections import OrderedDict
import string

tao_startup_param_list = [
    'beam_file;FILE;T;',
    'beam_init_position_file;FILE;T;',
    'beam_track_data_file;FILE;T;',
    'building_wall_file;FILE;T;',
    'data_file;FILE;T;',
    'hook_init_file;FILE;T;',
    'init_file;FILE;T;',
    'noinit;LOGIC;T;F',
    'lattice_file;FILE;T;',
    'plot_file;FILE;T;',
    'startup_file;FILE;T;',
    'noplot;LOGIC;T;F',
    'var_file;FILE;T;',
    'slice_lattice;FILE;T;',
    'disable_smooth_line_calc;LOGIC;T;F',
    'log_startup;LOGIC;T;F',
    'no_stopping;LOGIC;T;F',
    'rf_on;LOGIC;T;F',
]

#-------------------------------------------------

class tao_parameter():
    '''
    Basic class for holding the properties of a parameter in Tao.

    name:           The name of the parameter
    type:           "STR", "INT", "REAL", "LOGIC", "ENUM", etc...
    can_vary:       Either 'T', 'F', or 'I', indicating whether or not the
                    user may change the value of the paramter. 'I' indicates
                    that the parameter is to be ignored by the gui, except
                    to possibly be displayed as a sub-parameter to another
                    tao_parameter.
    value:          The value held in the parameter, should be of the
                    appropriate type for the specified param_type
                    (or 'T'/'F' for LOGIC)
    NOTE: It is unclear if sub_param is actually ever set by a Tao python command. -- DCS 11/2020
    sub_param:      The name of the sub_parameter associated with this parameter,
                    For example: ele_name has the sub parameter ix_ele.
    '''

    def __init__(self, param_name, param_type, can_vary, param_value, sub_param=None):
        # Enums and inums may have a prefix attached to their name, as in
        # axis^type.  In this case, the prefix is removed from the parameter name
        # and stored in the self.prefix variable
        if (param_type in ['ENUM', 'INUM']) and (param_name.count('^') == 1):
            self.prefix, self.name = param_name.split('^')
        else:
            self.prefix = None
            self.name = param_name
        self.type = param_type
        self.can_vary = (can_vary == 'T')
        self.is_ignored = (can_vary == 'I')
        self.sub_param = sub_param #associated sub_parameter (name)

        if param_type in ['STR', 'FILE', 'DAT_TYPE', 'DAT_TYPE_Z', 'DAT_TYPE_E',
                'REAL_ARR', 'ENUM', 'ENUM_Z', 'STRUCT', 'COMPONENT', 'SPECIES']:
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
        elif param_type == 'INUM':
            try:
                self.value = int(param_value)
            except:
                self.value = None
        else:
            print ('UNKNOWN PARAMETER TYPE: ' + param_type)

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return self.name + ';' + self.type + ';' + str(self.can_vary) + ';' + str(self.value)

    def get_component(self, comp_name):
        '''
        Looks for a component called comp_name in self, and returns the component
        with the matching name
        Returns None if comp_name was not found or if self.type != STRUCT
        '''
        if self.type != 'STRUCT':
            return None
        for param in self.value:
            if param.name == comp_name:
                return param.value
        return None

#

class InvalidParamError(Exception):
    '''
    Provides an exception for when a param_string is improperly formatted
    Examples of improper formatting include: not enough semicolons
    (3 minimum), bad sub-parameter formatting (for STRUCTs), etc
    '''
    pass

# An item in the parameter list is a string that looks like:

def tao_parameter_dict(param_list):
    '''
    Takes a list of strings, each string looks something like: 'param_name;STR;T;abcd'
    and returns a dictionary with keys being the param_name.
    Blank strings will be ignored.
    '''
    this_dict = OrderedDict()
    for param_str in param_list:
        if param_str.strip() == '': continue
        v = param_str.split(';')
        this_dict[v[0]] = str_to_tao_param(param_str)
    return this_dict

#

def str_to_tao_param(param_str):
    '''
    Takes a parameter string (EG: 'param_name;STR;T;abcd')
    and returns a tao_parameter instance
    param_str MUST have at least 3 semicolons
    '''
    v = param_str.split(';')
    if len(v) < 3:
        msg = str(param_str) + " is not a valid param_string (not enough semicolons)"
        raise InvalidParamError(msg)
    sub_param = None #default
    #TEMPORARY FIX
    if (len(v[2]) == 2) & (len(v) == 3):
        v.append(v[2][1])
        v[2] = v[2][0]
    ###
    # Special case: REAL_ARR (unknown length)
    if v[1] == "REAL_ARR":
        arr = []
        for i in range(len(v[3:])):
            x = v[3:][i]
            try:
                arr.append(float(x))
            except:
                if i==len(v[3:])-1: #last item, could be a related parameter name
                    if len(x) > 0:
                        sub_param = x
                else:
                    arr.append(float(0))
        v[3] = arr
    elif v[1] == 'STRUCT':
        n_comp = int(len(v[3:])/3)
        components = v[3:][:3*n_comp]
        c_list = [0]*n_comp
        for n in range(n_comp):
            c_name = components[3*n]
            c_type = components[3*n+1]
            c_val = components[3*n+2]
            c_list[n] = tao_parameter(c_name, c_type, v[2], c_val)
        v[3] = c_list
        if len(v[3:]) % 3 == 1: # one more item than expected --> it is a sub_param
            sub_param = v[-1]
    # Generic case sub_param: name;type;can_vary;value;sub_param
    if len(v) == 5:
        sub_param=v[4]
    return tao_parameter(v[0],v[1],v[2],v[3], sub_param)

#-------------------------------------------------
tao_startup_param_dict = tao_parameter_dict(tao_startup_param_list)
