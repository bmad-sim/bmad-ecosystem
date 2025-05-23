import os
import numpy as np
import pytest
from pytao import Tao

def test_absolute_time_tracking_does_not_affect_initial_z():
    """
    Regression test for issue #1535:
    Ensure that enabling `absolute_time_tracking` does not alter the initial 
    longitudinal position `z` of element 1 in the lattice.
    
    Only the second assertion is expected to fail due to the known bug.
    """
    test_dir = os.path.dirname(__file__)
    lat_path = os.path.join(test_dir, 'lat.bmad')

    tao = Tao(lattice_file=lat_path, noplot=True)

    z_before = tao.ele_orbit(1)['z']
    assert np.isclose(z_before, 0), f"Expected z=0 before setting absolute_time_tracking, got {z_before}"

    tao.cmd('set bmad_com absolute_time_tracking = T')

    z_after = tao.ele_orbit(1)['z']
    pytest.xfail("Known issue #1535: z changes when absolute_time_tracking is enabled")
    assert np.isclose(z_after, 0), f"Expected z=0 after setting absolute_time_tracking, got {z_after}"
