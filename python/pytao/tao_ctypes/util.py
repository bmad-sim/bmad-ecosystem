import os


def full_path(path):
    """
    Helper function to expand enviromental variables and return the absolute path
    """
    return os.path.abspath(os.path.expandvars(path))

def parse_bool(s):
    x = s.upper()[0]
    if x == 'T':
        return True
    elif x == 'F':
        return False
    else:
        raise ValueError ('Unknown bool: '+s) 


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


def pytype(type):
    """
    Returns python type for type =  REAL, INT, STR, ENUM, LOGIC
    """
    if type == 'REAL':
        f = float
    elif type == 'INT':
        f = int
    elif type == 'STR':
        f = str
    elif type == 'ENUM':
        f = str     
    elif type == 'SPECIES':
        f = str             
    elif type == 'INUM':
        f = int           
    elif type == 'LOGIC':
        f = bool              
    else:
        raise ValueError ('Unknown type: '+type)
    return f

def parse_tao_python_data1(line, clean_key=True):
    """
    Parses most common data output from a Tao>python command
    <component_name>;<type>;<is_variable>;<component_value>
    
    and returns a dict
    Example: 
        eta_x;REAL;F;  9.0969865321048662E+00
    parses to:
        {'eta_x':9.0969865321048662E+00}
    
    If clean key, the key will be cleaned up by replacing '.' with '_' for use as class attributes.
    
    See: tao_python_cmd.f90
    """
    dat = {}
    name, type, setable, val  = line.split(';')
    f = pytype(type)
    if f == bool:
        val = parse_bool(val)
    else:
        val = f(val)  
    if clean_key:
        name = name.replace('.', '_')
        
    return {name:val}

def parse_tao_python_data(lines, clean_key=True):
    """
    returns dict with data
    """
    dat = {}
    for l in lines:
        dat.update(parse_tao_python_data1(l, clean_key))
    return dat
    
    
    
def simple_lat_table(tao, ix_universe=1, ix_branch=0, which='model', who='twiss'):
    """
    Takes the tao object, and returns columns of parameters associated with lattice elements
     "which" is one of:
       model
       base
       design
     and "who" is one of:
       general         ! ele%xxx compnents where xxx is "simple" component (not a structure nor an array, nor allocatable, nor pointer).
       parameters      ! parameters in ele%value array
       multipole       ! nonzero multipole components.
       floor           ! floor coordinates.
       twiss           ! twiss parameters at exit end.
       orbit           ! orbit at exit end.
     Example:
    
    
    """
    # Form list of ele names
    cmd = 'python lat_ele_list '+str(ix_universe)+'@'+str(ix_branch)
    lines = tao.cmd(cmd)
    # initialize 
    ele_table = {}
    for x in lines:
        ix, name = x.split(';')
        # Single element information
        cmd = 'python lat_ele1 '+str(ix_universe)+'@'+str(ix_branch)+'>>'+str(ix)+'|'+which+' '+who
        lines2=tao.cmd(cmd)
        # Parse, setting types correctly
        ele = parse_tao_python_data(lines2)
        # Add name and index
        ele['name'] = name
        ele['ix_ele'] = int(ix)
        
        # Add data to columns 
        for key in ele:
            if key not in ele_table:
                ele_table[key] = [ele[key]]
            else:
                ele_table[key].append(ele[key])
        
        # Stop at the end ele
        if name == 'END': 
            break
    return ele_table
