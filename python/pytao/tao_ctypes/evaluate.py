from pytao import run_tao

from .tools import fingerprint

from h5py import File

import os


def evaluate_tao(settings,
                     run_commands=['set global track_type=beam'],
                     expressions=['lat::orbit.x[end]'],
                     input_file=None,
                     ploton=False,
                    
                     beam_archive_path=None,
                     workdir=None,
                     so_lib='',
                     verbose=False):
    """
    
    
    
    Expressions is a list of expressions that will be used to form the output
    
    
    beam::n_particle_loss[end]
    
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

        
    
    return output


