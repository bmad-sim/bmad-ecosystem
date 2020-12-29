# Make superimosed markers for extra beam saving
import numpy as np

def make_markers(slist, filename=None, ref=None):
    """
    Makes markers relative to ref ele.
    
    If filename is given, the lines will be written to ta file. 
    
    """
    lines = []
    
    if ref:
        ref = f'ref={ref},'
    else:
        ref = ''
        
    for i, s in enumerate(slist):
        line = f'm_{i:03}: marker, superimpose, {ref} offset = {s}'
        lines.append(line)
        
    if filename:
        with open(filename, 'w') as f:
            for line in lines:
                f.write(line+'\n')
    return lines