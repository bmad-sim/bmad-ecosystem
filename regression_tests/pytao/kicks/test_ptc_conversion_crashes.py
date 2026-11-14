"""
Lattice configurations that used to terminate Tao while converting an element to PTC (issue #2245).

1. A zero length element with an electric multipole: ele_to_fibre gave a fatal error and the caller
   then crashed.  An electric field over zero length gives no kick, so the field is now ignored.
2. A zero length traveling wave lcavity: PTC stopped with "TWCAVITY MUST HAVE A LENGTH".
3. A planar wiggler with an electric multipole: the wiggler was converted to a PTC sbend, and the
   wiggler setup then wrote to unallocated wiggler data.
"""

from __future__ import annotations

import math

import pytest
from pytao import SubprocessTao

HEADER = """\
parameter[particle] = electron
parameter[geometry] = open
beginning[e_tot] = 1e9
beginning[beta_a] = 10
beginning[beta_b] = 10
"""

CASES = {
    "sextupole_zero": ("sextupole, l = 0, scale_multipoles = F", "b0_elec", 20),
    "kicker_zero": ("kicker, l = 0, scale_multipoles = F", "b0_elec", 20),
    "rfcavity_zero": ("rfcavity, l = 0, rf_frequency = 500e6, voltage = 0", "b0_elec", 20),
    "lcavity_zero_tw": (
        "lcavity, l = 0, rf_frequency = 500e6, voltage = 0, cavity_type = traveling_wave", "hkick", 1e-8),
    "wiggler_planar": (
        "wiggler, l = 0.5, l_period = 0.1, field_calc = planar_model, b_max = 1e-4, scale_multipoles = F",
        "b0_elec", 20),
}


@pytest.mark.parametrize("method", ["symp_lie_ptc", "taylor"])
@pytest.mark.parametrize("case", CASES, ids=list(CASES))
def test_ptc_conversion_does_not_crash(tmp_path, case: str, method: str):
    definition, attr, value = CASES[case]
    lattice_file = tmp_path / "lat.bmad"
    lattice_file.write_text(HEADER + f"ELE: {definition}, tracking_method = {method}\nlat: line = (ELE)\nuse, lat\n")
    tao = SubprocessTao(lattice_file=str(lattice_file), noplot=True)
    try:
        tao.cmd("set particle_start x = 5e-3")
        tao.cmd("set particle_start y = 3e-3")
        tao.cmd(f"set ele ELE {attr} = {value}")
        orb = tao.ele_orbit("ELE")
        assert orb["state"] == "Alive"
        assert all(math.isfinite(orb[c]) for c in ("x", "px", "y", "py", "z", "pz"))
    finally:
        tao.close_subprocess()
