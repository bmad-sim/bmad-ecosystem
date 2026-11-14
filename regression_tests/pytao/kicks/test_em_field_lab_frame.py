"""
Check the lab-frame field that em_field_calc returns for a tilted element (issue #2239).

hkick and vkick are lab-frame kicks, so their field must not change with the tilt.  Magnetic and
electric multipoles are body-frame quantities, so their lab-frame field rotates with the tilt.
The field is evaluated with "pipe em_field", which calls em_field_calc with local_ref_frame = F.
"""

from __future__ import annotations

import math

import pytest
from pytao import SubprocessTao

LATTICE = """\
beginning[e_tot] = 5.8e9
beginning[beta_a] = 10
beginning[beta_b] = 10
parameter[particle] = electron
parameter[geometry] = open
QQ: quadrupole, l = 0.6, scale_multipoles = F
lat: line = (QQ)
use, lat
"""

TILT = 0.3


@pytest.fixture(scope="module")
def tao(tmp_path_factory: pytest.TempPathFactory):
    lattice_file = tmp_path_factory.mktemp("em_field") / "lat.bmad"
    lattice_file.write_text(LATTICE)
    t = SubprocessTao(lattice_file=str(lattice_file), noplot=True)
    yield t
    t.close_subprocess()


def field(tao: SubprocessTao, attr: str, value: float, tilt: float) -> tuple[float, ...]:
    tao.cmd(f"set ele QQ tilt = {tilt}")
    tao.cmd(f"set ele QQ {attr} = {value}")
    try:
        line = tao.cmd("pipe em_field QQ|model 0 0 0.3 0")[0]
    finally:
        tao.cmd(f"set ele QQ {attr} = 0")
        tao.cmd("set ele QQ tilt = 0")
    return tuple(float(v) for v in line.split(";"))  # Bx, By, Bz, Ex, Ey, Ez


def rotate(vx: float, vy: float, angle: float) -> tuple[float, float]:
    c, s = math.cos(angle), math.sin(angle)
    return vx * c - vy * s, vx * s + vy * c


@pytest.mark.parametrize("attr, value", [("hkick", 2e-6), ("vkick", -3e-6)])
def test_hvkick_field_ignores_tilt(tao: SubprocessTao, attr: str, value: float):
    f0 = field(tao, attr, value, 0.0)
    f1 = field(tao, attr, value, TILT)
    scale = max(abs(v) for v in f0)
    assert scale > 0
    assert f1 == pytest.approx(f0, abs=1e-12 * scale)


@pytest.mark.parametrize("attr, value, ix", [("b0", 1e-6, 0), ("a0", 1e-6, 0), ("b0_elec", 1e5, 3), ("a0_elec", 1e5, 3)])
def test_multipole_field_rotates_with_tilt(tao: SubprocessTao, attr: str, value: float, ix: int):
    f0 = field(tao, attr, value, 0.0)
    f1 = field(tao, attr, value, TILT)
    want = rotate(f0[ix], f0[ix + 1], TILT)
    scale = math.hypot(f0[ix], f0[ix + 1])
    assert scale > 0
    assert f1[ix : ix + 2] == pytest.approx(want, abs=1e-12 * scale)
