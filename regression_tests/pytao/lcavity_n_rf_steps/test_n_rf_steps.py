from pathlib import Path
from pytao import SubprocessTao
import numpy as np
import pytest


def test_beam_energy_n_rf_steps():
    """
    Regression test for issue #1619:
    Ensure that reference particle energy is changed correctly for multipass cavities using bmad-standard tracking with
    `n_rf_steps` set.
    """
    lat_path = Path(__file__).parent / "lat1.bmad"
    assert lat_path.is_file(), f"Lattice file not found: {lat_path}"

    # Check beam energy
    with SubprocessTao(lattice_file=str(lat_path), noplot=True) as tao:
        np.testing.assert_allclose(
            tao.ele_gen_attribs(r"cav\1")["E_TOT"],
            tao.ele_gen_attribs(r"cav\1")["E_TOT_START"] + tao.ele_gen_attribs(r"cav\1")["VOLTAGE"],
            err_msg="Expected beam energy at end of multipass cavity with n_rf_steps != 0 to be start energy plus voltage."
        )


def test_ref_time_n_rf_steps():
    """
    Regression test for issue #1642:
    Ensure that reference time is changed correctly for multipass cavities using bmad-standard tracking with
    `n_rf_steps` set.
    """
    lat_path = Path(__file__).parent / "lat1.bmad"
    assert lat_path.is_file(), f"Lattice file not found: {lat_path}"

    # Check beam energy
    with SubprocessTao(lattice_file=str(lat_path), noplot=True) as tao:
        assert tao.ele_gen_attribs(r"cav\1")["DELTA_REF_TIME"] > 0
        assert tao.ele_gen_attribs(r"cav\2")["DELTA_REF_TIME"] > 0    


def test_segfault_multiple_tracking_mode_n_rf_steps():
    """
    Regression test for issue #1629:
    Ensure that tao doesn't crash when loading lattice with multiple tracking modes in multipass cavity with n_rf_steps != 0.
    """
    lat_path = Path(__file__).parent / "lat2.bmad"
    assert lat_path.is_file(), f"Lattice file not found: {lat_path}"

    # Attempt to read basic parameter from lattice without tao crashing
    with SubprocessTao(lattice_file=str(lat_path), noplot=True) as tao:
        np.testing.assert_allclose(
            tao.ele_gen_attribs(r"beginning")["E_TOT"],
            10e6,
            err_msg="Could not read starting energy of lattice."
        )


# Run on eight t_offsets ignoring t_offset = 0.0 and t_offset = RF_PERIOD
@pytest.mark.parametrize("t_offset", np.linspace(0, 1/1e9, 10)[1:-1])
def test_n_rf_steps_patch(t_offset):
    """
    Regression test for issue #1693
    Confirm second pass cavity timing does not depend on T_OFFSET of patch element placed before.
    """
    lat_path = Path(__file__).parent / "lat3.bmad"
    assert lat_path.is_file(), f"Lattice file not found: {lat_path}"

    def ptot(orbit):
        return (1+orbit['pz']) * orbit['p0c']
    
    # Check beam energy
    with SubprocessTao(lattice_file=str(lat_path), noplot=True) as tao:
        # Grab energy of unmodified lattice
        energy1 = ptot(tao.ele_orbit(r"cav\2"))

        # Apply t_offset to the patch
        tao.cmd(f"set ele p t_offset = {t_offset}")
        
        # Grab second pass energy and compare
        np.testing.assert_allclose(
            energy1,
            ptot(tao.ele_orbit(r"cav\2")),
            err_msg="Expected second pass energy to not vary with PATCH[T_OFFSET]."
        )
