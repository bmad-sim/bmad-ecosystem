from pathlib import Path
import numpy as np
from pytao import Tao

def assert_dicts_allclose(d1, d2, rtol=1e-5, atol=1e-8):
    """
    Assert that all matching array entries in two particle dicts are close.

    Parameters
    ----------
    d1 : dict
        First particle dictionary.
    d2 : dict
        Second particle dictionary.
    rtol : float
        Relative tolerance for np.allclose.
    atol : float
        Absolute tolerance for np.allclose.

    Raises
    ------
    AssertionError
        If any entry mismatches or types differ.
    """
    assert d1.keys() == d2.keys(), f"Keys differ: {d1.keys()} vs {d2.keys()}"
    
    for key in d1:
        v1 = d1[key]
        v2 = d2[key]
        
        if isinstance(v1, np.ndarray):
            assert isinstance(v2, np.ndarray), f"Type mismatch for key '{key}': ndarray vs {type(v2)}"
            assert v1.shape == v2.shape, f"Shape mismatch for key '{key}': {v1.shape} vs {v2.shape}"
            assert np.allclose(v1, v2, rtol=rtol, atol=atol), f"Array values differ for key '{key}'"
        else:
            assert v1 == v2, f"Scalar values differ for key '{key}': {v1} vs {v2}"



def test_sr_wakes_onoff():
    lattice_file = Path(__file__).parent / "lat.bmad"
    startup_file = Path(__file__).parent / "tao.startup"

    
    tao = Tao(lattice_file=lattice_file, startup_file=startup_file, noplot=True)
   
    tao.cmd('set bmad sr_wakes_on = F')
    d1 = tao.bunch_data('end')
    
    tao.cmd('set bmad sr_wakes_on = T')
    d2 = tao.bunch_data('end')

    assert_dicts_allclose(d1, d2)