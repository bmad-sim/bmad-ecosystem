import pytest
import numpy as np
from pytao import Tao

@pytest.mark.xfail(reason="Known issue #1535: z-coordinate changes unexpectedly")
def test_absolute_time_tracking_does_not_affect_z():
    """
    Regression test for issue #1535:
    Ensure that enabling `absolute_time_tracking` does not alter the initial 
    longitudinal position `z` of element 1 in the lattice.
    
    This checks that switching `bmad_com.absolute_time_tracking` from False (default) 
    to True does not affect the reported `z` position of the first element.
    """
    # Initialize Tao without plotting
    tao = Tao(lattice_file='lat.bmad', noplot=True)

    # Get initial z before changing time tracking setting
    z_before = tao.ele_orbit(1)['z']
    assert np.isclose(z_before, 0), f"Expected z=0 before setting absolute_time_tracking, got {z_before}"

    # Enable absolute time tracking
    tao.cmd('set bmad_com[0].absolute_time_tracking = T')

    # Get z after setting absolute time tracking
    z_after = tao.ele_orbit(1)['z']
    assert np.isclose(z_after, 0), f"Expected z=0 after setting absolute_time_tracking, got {z_after}"
