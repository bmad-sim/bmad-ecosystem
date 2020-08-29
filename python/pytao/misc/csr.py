# Functions to handle CSR
import numpy as np

def parse_csr_wake(file, prepend_index=True):
    """
    Parses CSR wake file and returns a dict of elements with
    
    
    ix_ele: index of the element
    
    s_positions: np.array of s possitions correponding to the steps in data
    
    labels: list of string labels for columns
    
    data: np.array of (step, bin, column)
        column 0: Z
        column 1: Charge/Meter
        column 2: CSR_Kick/m
        column 3: I_CSR/m
        column 4: S_Source
    """
    with open(file) as f:
        eles = {}
        line = f.readline()

        data = None # Start 
        while line:
        
            
            # Check for new ele
            if line.startswith('!------------------------------------------------------------'):
                line = f.readline()    
                _, ix, ele_name = line.split()
                ix = int(ix)
                
                # Add index to key for uniqueness
                if prepend_index:
                    key = f'{ix}:{ele_name}'
                else:
                    key = ele_name
                
                eles[key] = {'data':[], 's_positions':[], 'ix_ele':ix,
                                  'ele_name':ele_name,
                                  'labels':['z',
                                            'charge_per_meter',
                                            'csr_kick_per_meter',
                                            'I_csr_per_meter',
                                            's_source']}
                # Current ele
                ele = eles[key]
                
                # point to current data
                ele_steps= ele['data'] = []
                
              #  print('new ele', ele_name)
                line = f.readline()
                continue
            #print(line)
            
            # Check for new step
            if line.startswith('!#-----------------------------'):
                
                line = f.readline()
                step = int(line.split('index:')[-1])
                line = f.readline()
                s_position = float(line.split()[-1])
                f.readline() # header
                
                # new data 
                data = []
                ele_steps.append(data)
                ele['s_positions'].append(s_position)
               # print('data: ', data)
                
               # print('new step:', step)
                line = f.readline()
                continue
                
           # print(line)    
            data.append(line)
                
            
            
            line = f.readline()
                        
    # Cast to numpy
    for _, ele in eles.items():
        ele['data'] = np.array([np.loadtxt(data) for data in ele['data']])
        ele['s_positions'] = np.array(ele['s_positions'] )
    return eles



def process_csr_wake_data(data):
    """
    Process csr_wake data by concatenating the data. 
    
    Returns a dict of np.array with labels and shapes:
    
        s_position          (steps,)      
        ix_ele              (steps,)
        z                   (steps, bins)
        charge_per_meter    (steps, bins)
        csr_kick_per_meter  (steps, bins)
        I_csr_per_meter     (steps, bins)
        s_source            (steps, bins)
    
    and a dict:
        ele_name_from_ix
    that returns the ele_name from an ix_ele for convenience. 
    
    
    """
    
    # Sort by ix_ele
    skeys = sorted(data, key=lambda n: data[n]['ix_ele'] )
    labels = data[skeys[0]]['labels'].tolist()
    #return labels

    output = {}
    # Concatenate all
    dat = np.concatenate([data[key]['data'] for key in skeys ])
    output['s_position'] = np.concatenate([data[key]['s_positions'] for key in skeys ])
    
    # Fill out ix_ele
    output['ix_ele'] = np.concatenate([np.full(len(data[key]['s_positions']),  data[key]['ix_ele']) for key in skeys])
    
    # ix-> name mapping
    output['ele_name_from_ix'] = {data[key]['ix_ele']:key for key in data }
    
    # Split columns
    for i, key in enumerate(labels):
        output[key] = dat[:,:,i]
    
    return output
    

    
def csr_wake_stats_at_step(data, step=0):
    """
    Calculates stats from a single step of processed csr_wake data.
    
    Useful for plotting.
    
    Returns tuple:
    
        avkick, stdkick, avz, stdz
    
    """
 
    z = data['z'][step,:]

    
    zmax = z.max()
    zmin = z.min()
    nz = len(z)
    dz = z.ptp()/(nz-1)
    density = data['charge_per_meter'][step,:]
    dsum = np.sum(density)
    qtot = dsum*dz
    density = density/dsum
 
    
    avz = np.sum(z*density)
    stdz = np.sqrt(np.sum(z**2*density) - avz**2)
    
    kick = data['csr_kick_per_meter'][step,:]
    avkick = np.sum(kick*density)
    stdkick = np.sqrt(np.sum( kick**2*density) - avkick**2)
    
    return avkick, stdkick, avz, stdz    


def write_csr_wake_data_h5(h5, data, name=None):
    """
    Write parsed csr_wake data to an open h5 handle
    """
    if name:
        h5 = h5.create_group(name)
    
    for ele in data:
        g = h5.create_group(ele)
        for k, v in data[ele].items():
            # Write np.array as datasets, otherwise to attributes
            if isinstance(v, np.ndarray):
                g[k] = v
            else:
                g.attrs[k] = v
                
def read_csr_wake_data_h5(h5, name=None):
    """
    Read csr_wake data from h5 file.
    
    See: write_csr_wake_data_h5
    """
    
    
    if name:
        h5 = h5[name]
    dat = {}
    for g in h5:
        dat[g] = {}
        ele = dat[g]
        # Read all datasets and attributes
        for key in h5[g]:
            ele[key] = h5[g][key][:]
        ele.update(dict(h5[g].attrs))
    return dat



