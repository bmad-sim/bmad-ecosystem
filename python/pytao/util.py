


def parse_tao_lat_ele_list(lines):
    """
    returns mapping of names to index
    
    TODO: if elements are duplicated, this returns only the last one.
    
    Example: 
    ixlist = parse_tao_lat_ele_list(tao.cmd('python lat_ele_list 1@0'))
    """
    ix = {}
    for l in lines:
        index, name = l.split(';')
        ix[name] = int(index)
    return ix




def parse_tao_python_data1(line):
    """
    Parses most common data output from a Tao>python command
    <component_name>;<type>;<is_variable>;<component_value>
    
    and returns a dict
    Example: 
        eta_x;REAL;F;  9.0969865321048662E+00
    parses to:
        {'eta_x':9.0969865321048662E+00}
    
    See: tao_python_cmd.f90
    """
    dat = {}
    name, type, setable, val  = line.split(';')
    if type == 'REAL':
        f = float
    elif type == 'INT':
        f = int
    elif type == 'STR':
        f = str
    elif type == 'ENUM':
        f = str     
    elif type == 'LOGIC':
        f = bool              
    else:
        raise ValueError ('Unknown type: '+type)
    return {name:f(val)}

def parse_tao_python_data(lines):
    """
    returns dict with data
    """
    dat = {}
    for l in lines:
        dat.update(parse_tao_python_data1(l))
    return dat