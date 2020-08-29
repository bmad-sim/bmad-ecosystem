from pytao import run_tao
from pytao.misc.csr import parse_csr_wake, write_csr_wake_data_h5

from .tools import fingerprint

from h5py import File

import os


def evaluate_tao(settings,
                     run_commands=['set global track_type=beam'],
                     expressions=['lat::orbit.x[end]'],
                     input_file=None,
                     ploton=False,
                    
                     beam_archive_path=None,
                     archive_csr_wake=False,
                     workdir=None,
                     so_lib='',
                     verbose=False):
    """
    
    settings: dict of set_command:value where set_command is a string.
        Example:
            'global:track_type':'beam'
                will issue command:
            set global track_type = beam
    
    run_commands: list of command strings that will be executed.
    
    expressions: list of expression strings that will be used to form the output.
    
    beam_archive_path: if given, the all of the saved beams will be written to
        a file named by a fingerprint (hash) of the inputs into path beam_archive_path.
        This uses the command:
            write beam -at *
        which writes ALL of the bunches that are saved using the beam_saved_at list in beam_init.
        
    archive_csr_wake: if given, will look for csr_wake.dat, parse, and archive to the h5 file above.
    
    Returns a dict of expression:value, according to the expressions above, as well as 
        beam_archive if a  beam_archive_path was given.
    
    
    
    Example:
    
    evaluate_tao(settings={}, 
                   input_file=tao.init', 
                   run_commands=['set global track_type=beam'],
                   expressions = ['lat::orbit.x[FF.PIP02A]', 'beam::norm_emit.x[end]'],
                   ploton=False, 
                   beam_archive_path = '.')
                   
    Returns:
    
        {'lat::orbit.x[FF.PIP02A]': 0.0,
         'beam::norm_emit.x[end]': 9.9982321550206e-07,
         'beam_archive': /path/to/bmad_beam_7fd6d30ac45a3d8c0d45112f4b569dee.h5'}
    
    
    
    See: run_tao
    
    """

    M = run_tao(settings=settings,
                run_commands=run_commands,
                input_file=input_file,
                ploton=ploton,
                workdir=workdir,
                so_lib=so_lib,
                verbose=verbose)
    
    
    output = {}
    
    for expression in expressions:
        try:
            val = M.evaluate(expression)
        except:
            print(f'error with {expression}')
            val = None
        output[expression] = val    
    
    if beam_archive_path:
        ff = fingerprint({'input_file':input_file, 'settings':settings})
        beam_archive_path = os.path.expandvars(beam_archive_path)
        beam_archive = os.path.abspath(os.path.join(beam_archive_path, f'bmad_beam_{ff}'+'.h5'))
        if verbose:
            print('Archiving beam to', beam_archive)
        M.cmd(f'write beam -at * {beam_archive}')    
        output['beam_archive'] = beam_archive
        

        # Reopen and attach settings
        assert os.path.exists(beam_archive), 'No archive was written. Perhaps there was no beam?'
        
        with File(beam_archive, 'r+') as h5:
            # Input
            g = h5.create_group('input')
            g.attrs['input_file'] = input_file
            
            #g.attrs['input_file'] = input_file
            
            
            # Settings
            g = h5.create_group('settings')
            for k, v in settings.items():
                g.attrs[k] = v
        
            g = h5.create_group('expressions')
            for k, v in output.items():
                if v:
                    g.attrs[k] = v
                    
            # CSR wake
            csr_wake_file = os.path.join(M.path, 'csr_wake.dat')

            if archive_csr_wake and os.path.exists(csr_wake_file):
                csr_wake_data = parse_csr_wake(csr_wake_file)
                write_csr_wake_data_h5(h5, csr_wake_data, name='csr_wake')
                

        
    
    return output


