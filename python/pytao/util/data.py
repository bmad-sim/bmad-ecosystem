'''

'''





def create_d2_data_cmd(d2_name, d1_names, upper_bounds=[], lower_bounds=None):
    """
    Helper function to create a d2_data with several d1_data.
    
    Returns the command string for use in Tao. 
    
    Example usage:
        cmd = create_d2_data_cmd('X', ['A', 'B'], [2,4] )
        tao.cmd(cmd)
    """

    n = len(d1_names)
    if not lower_bounds:
        lower_bounds = n*[1]
    
    cmd = f'python data_d2_create {d2_name}^^{n}'
    for name, lower_bound, upper_bound in zip(d1_names, lower_bounds, upper_bounds):
        cmd+= f'^^{name}^^{lower_bound}^^{upper_bound}'
        
    return cmd



def create_datum_cmd(
 
    datum_name='d2.d1[1]',
    data_type='beta.a',
    ele_ref_name='',
    ele_start_name='',
    ele_name='end',
    merit_type='target',
    meas=0,
    good_meas=False,
    ref=0,
    good_ref=False,
    weight=0,
    good_user=True,
    data_source='',
    eval_point='end',
    s_offset=0,
    ix_bunch=0,
    invalid_value=0,
    spin_n0_x=0,
    spin_n0_y=0,
    spin_n0_z=0 
                ):
    """
    Helper function to create a datum, with defaults.
    
    Returns the command string for use in Tao. 
    
    Example usage:
        cmd=create_datum_cmd('X.A[1]', meas=1, good_meas=True) 
        tao.cmd(cmd)
        tao.cmd('python data_set_design_value') # Needed to initialize
        
    """
    # Split for readablity 
    cmd = 'python datum_create '\
            f'{datum_name}^^{data_type}^^{ele_ref_name}^^{ele_start_name}^^{ele_name}^^{merit_type}'\
            f'^^{meas}^^{good_meas}^^{ref}^^{good_ref}^^{weight}^^{good_user}^^{data_source}^^{eval_point}'\
            f'^^{s_offset}^^{ix_bunch}^^{invalid_value}^^{spin_n0_x}^^{spin_n0_y}^^{spin_n0_z}'
    
    return cmd




def create_variable_cmd(
    var_name='',
    ele_name='',
    attribute='',
    universes='1',
    weight=0,
    step=0,
    low_lim=0,
    high_lim=0,
    merit_type='target',
    good_user=True,
    key_bound=False,
    key_delta=1
):
    """
    TODO: this doesn't work. 
    
    Create a single variable
    Command syntax:
    python var_create {var_name}^^{ele_name}^^{attribute}^^{universes}^^{weight}^^{step}^^{low_lim}^^{high_lim}^^
                                                                         {merit_type}^^{good_user}^^{key_bound}^^{key_delta}
    
    """
    
    
    # Split for readablity 
    cmd = 'python var_create '\
        f'{var_name}^^{ele_name}^^{attribute}^^{universes}'\
        f'^^{weight}^^{step}^^{low_lim}^^{high_lim}'\
        f'^^{merit_type}^^{good_user}^^{key_bound}^^{key_delta}'
    return cmd