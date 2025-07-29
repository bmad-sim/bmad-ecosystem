from pathlib import Path
import numpy as np
from pytao import SubprocessTao


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
    Regression test for issue #1619:
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
