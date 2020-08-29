from .csr import csr_wake_stats_at_step

import matplotlib.pyplot as plt
import numpy as np


def plot_csr_wake(data, step=0, **kwargs):
    """
    Plots a step from process CSR wake data
    
    
    A simple interactive plot can be made with:
    
        from ipywidgets import interact
        def plot1(step=0):
            plot_csr_wake(data, step=step)
        nstep = len(data['z'])
        interact(plot1, step=(0, nstep-1, 1) )
    
    """
    
    fig, ax = plt.subplots(**kwargs)

    
    z = data['z'][step,:]
    density =  data['charge_per_meter'][step,:]
    kick = data['csr_kick_per_meter'][step,:]
    
    # Full bounds for plot
    zmin = data['z'].min()
    zmax = data['z'].max()
    
    kmin = np.min(data['csr_kick_per_meter'])
    kmax = np.max(data['csr_kick_per_meter'])
    
    #nz = len(z)
    #dz = z.ptp()/(nz-1)
    #dsum = np.sum(density)
    #qtot = dsum*dz 
    ## Normalize for stats
    #density = density/dsum
            
    ax.set_xlabel('z (µm)')
    ax.set_ylabel('CSR Kick/m')
    
    
    
    zscale = 1e6 # m -> µm
    ax.set_xlim(zmin*zscale, zmax*zscale)
    ax.set_ylim(kmin, kmax)
    
   
    ax2 = ax.twinx()
    ax2.set_ylabel('density')
    ax2.fill(z*zscale, density, color='grey', alpha=0.5)
    
    ax.plot(z*zscale, kick, color='blue')
    
    
    
def plot_csr_stats(data, **kwargs):
    """
    Plots csr wake stats along s from processed csr_wake data
    """
    
    s_pos = data['s_position']
    nsteps = len(data['z'])
    
    stats = np.array([csr_wake_stats_at_step(data, step=i) for i in range(nsteps)])
    
    fig, ax = plt.subplots(**kwargs)
    
    ax.set_xlabel('s (m)')
    ax.set_ylabel('CSR Kick/m')
    ax.plot(s_pos, stats[:,0], color='black', label='Average Wake')
    ax.plot(s_pos, stats[:,1], color='red', label='std Wake')
    
    ax2 = ax.twinx()
    ax2.plot(s_pos, stats[:,3]*1e15/299792458, color='blue', label='$\sigma_z/c (fs)$')
    ax2.set_ylabel('$\sigma_z/c (fs)$')
    ax.legend()

    
    
 
