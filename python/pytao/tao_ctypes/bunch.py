
from pytao.tao_ctypes import util
import numpy as np


def get_bunch_stats(tao, ele=0, ix_universe=1, ix_branch=0, ix_bunch=1, verbose=False):
    """
    Gets bunch statistics by calling:
    
    python bunch1 {ix_universe}@{ix_branch}>>{ix_ele}|{which} {ix_bunch}
    
    Returns a dict
    
    """
    cmd = f'python bunch1 {ix_universe}@{ix_branch}>>{ele}|model {ix_bunch}'
    if verbose:
        print(cmd)
    stats = tao.cmd(cmd)
    
    # Form into dict
    stats=util.parse_tao_python_data(stats)
    
    return stats




def get_bunch_data(tao, ele=0, ix_universe=1, ix_branch=0, ix_bunch=1, verbose=False):
    """
    Returns bunch data in openPMD-beamphysics format/notation.
    
    Note that Tao's 'write beam' will also write a proper h5 file in this format.
    
    Expected usage:
        data = get_bunch_data(ta0, ele='end')
        from pmd_beamphysics import ParticleGroup
        P = ParicleGroup(data=data)
    
    """
    
    # Get species
    stats = get_bunch_stats(tao, ele=ele, ix_universe=ix_universe, ix_branch=ix_branch, ix_bunch=ix_bunch)
    species = stats['species']
    
    
    dat = {}
    for key in ['x', 'px', 'y', 'py',  't', 'pz', 'p0c', 'charge', 'state']:
        
        cmd = f'python bunch1 {ix_universe}@{ix_branch}>>{ele}|model {ix_bunch} {key}'
        if verbose:
            print(cmd)
            
        # rename this one
        if key == 'charge':
            key = 'weight'
            
        if key == 'state':    
            key = 'status'
            dat[key] = tao.cmd_integer(cmd)
        else:
            dat[key] = tao.cmd_real(cmd)
    
    
    # Remove normalizations
    p0c = dat.pop('p0c')
    
    # px from Bmad is px/p0c
    
    # pz from Bmad is delta = p/p0c -1. 
    # pz = sqrt( (delta+1)**2 -px**2 -py**2)*p0c
    dat['pz'] = np.sqrt((dat['pz'] + 1)**2 - dat['px']**2 - dat['py']**2) * p0c
    dat['px'] = dat['px']*p0c
    dat['py'] = dat['py']*p0c

    # z = 0 by definition
    dat['z'] = np.full(len(dat['x']), 0)
        
    dat['species'] = species.lower()
    
    return dat