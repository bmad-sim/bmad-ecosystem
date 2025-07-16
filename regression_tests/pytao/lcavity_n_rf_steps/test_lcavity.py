from pathlib import Path
import numpy as np
from pytao import SubprocessTao


def test_beam_energy_n_rf_steps():
    """
    Regression test for issue #1619:
    Ensure that reference particle energy is changed correctly for multipass cavities using bmad-standard tracking with
    `n_rf_steps` set.
    """
    lat_path = Path(__file__).parent / "lat.bmad"
    assert lat_path.is_file(), f"Lattice file not found: {lat_path}"

    # Check beam energy
    with SubprocessTao(lattice_file=str(lat_path), noplot=True) as tao:
        np.testing.assert_allclose(
            tao.ele_gen_attribs(r"cav\1")["E_TOT"],
            tao.ele_gen_attribs(r"cav\1")["E_TOT_START"] + tao.ele_gen_attribs(r"cav\1")["VOLTAGE"],
            err_msg="Expected beam energy at end of multipass cavity with n_rf_steps != 0 to be start energy plus voltage."
        )
