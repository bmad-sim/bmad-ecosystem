"""
Check the lab-frame field that em_field_calc returns for a tilted element (issue #2239).

For historical reasons, hkick and vkick are not rotated by the element tilt (ref_tilt for a bend), except
for kicker and elseparator elements. So for the elements here their field must not change with the tilt.
Magnetic and electric multipoles are body-frame quantities, so their lab-frame field rotates with the tilt.
The field is evaluated with "pipe em_field", which calls em_field_calc with local_ref_frame = F.

Each field is the difference between the field with the attribute set and with it zero, which removes
the bend and wiggler fields.
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
BB: sbend, l = 0.6, g = 0.01, scale_multipoles = F
BX: sbend, l = 0.6, g = 0.01, scale_multipoles = F, exact_multipoles = vertically_pure
WW: wiggler, l = 1.0, field_calc = planar_model, b_max = 0.1, l_period = 0.1, tracking_method = runge_kutta,
    scale_multipoles = F
lat: line = (QQ, BB, BX, WW)
use, lat
"""

TILT = 0.3

# Element name, tilt attribute, s position (body coords) at which to evaluate the field.
ELEMENTS = [
    ("QQ", "tilt", 0.3),
    ("BB", "ref_tilt", 0.3),
    ("BX", "ref_tilt", 0.3),
    ("WW", "tilt", 0.5),
]
ELEMENT_IDS = ["quadrupole", "sbend", "sbend_exact_multipoles", "planar_wiggler"]


@pytest.fixture(scope="module")
def tao(tmp_path_factory: pytest.TempPathFactory):
    lattice_file = tmp_path_factory.mktemp("em_field") / "lat.bmad"
    lattice_file.write_text(LATTICE)
    t = SubprocessTao(lattice_file=str(lattice_file), noplot=True)
    yield t
    t.close_subprocess()


def em_field(tao: SubprocessTao, ele: str, x: float, s: float) -> list[float]:
    line = tao.cmd(f"pipe em_field {ele}|model {x} 0 {s} 0")[0]
    return [float(v) for v in line.split(";")]  # Bx, By, Bz, Ex, Ey, Ez


def field(tao: SubprocessTao, ele: str, tilt_attr: str, s: float, attr: str, value: float,
          tilt: float, settings: tuple[str, ...] = (), x: float = 0) -> list[float]:
    """Lab-frame field of attr at (x, 0, s) with the given tilt and extra element settings."""
    tao.cmd(f"set ele {ele} {tilt_attr} = {tilt}")
    for setting in settings:
        tao.cmd(f"set ele {ele} {setting}")
    try:
        f_base = em_field(tao, ele, x, s)
        tao.cmd(f"set ele {ele} {attr} = {value}")
        f_attr = em_field(tao, ele, x, s)
    finally:
        tao.cmd(f"set ele {ele} {attr} = 0")
        for setting in settings:
            tao.cmd(f"set ele {ele} {setting.split('=')[0]} = 0")
        tao.cmd(f"set ele {ele} {tilt_attr} = 0")
    return [b - a for a, b in zip(f_base, f_attr)]


def rotate(vx: float, vy: float, angle: float) -> tuple[float, float]:
    c, s = math.cos(angle), math.sin(angle)
    return vx * c - vy * s, vx * s + vy * c


@pytest.mark.parametrize("ele, tilt_attr, s", ELEMENTS, ids=ELEMENT_IDS)
@pytest.mark.parametrize("attr, value", [("hkick", 2e-6), ("vkick", -3e-6)])
def test_hvkick_field_ignores_tilt(tao: SubprocessTao, ele: str, tilt_attr: str, s: float, attr: str, value: float):
    f0 = field(tao, ele, tilt_attr, s, attr, value, 0.0)
    f1 = field(tao, ele, tilt_attr, s, attr, value, TILT)
    scale = max(abs(v) for v in f0)
    assert scale > 0
    assert f1 == pytest.approx(f0, abs=1e-9 * scale)


@pytest.mark.parametrize("ele, tilt_attr, s", ELEMENTS, ids=ELEMENT_IDS)
@pytest.mark.parametrize("attr, value, ix", [("b0", 1e-6, 0), ("a0", 1e-6, 0), ("b0_elec", 1e5, 3), ("a0_elec", 1e5, 3)])
def test_multipole_field_rotates_with_tilt(tao: SubprocessTao, ele: str, tilt_attr: str, s: float, attr: str,
                                           value: float, ix: int):
    f0 = field(tao, ele, tilt_attr, s, attr, value, 0.0)
    f1 = field(tao, ele, tilt_attr, s, attr, value, TILT)
    want = rotate(f0[ix], f0[ix + 1], TILT)
    scale = math.hypot(f0[ix], f0[ix + 1])
    assert scale > 0
    assert f1[ix : ix + 2] == pytest.approx(want, abs=1e-9 * scale)


@pytest.mark.parametrize("ele", ["BB", "BX"], ids=["sbend", "sbend_exact_multipoles"])
def test_bend_multipole_field_follows_offset(tao: SubprocessTao, ele: str):
    """At the bend center, the b1 field of a bend offset by x_offset at lab x + x_offset equals the
    field of the unoffset bend at x."""
    f0 = field(tao, ele, "ref_tilt", 0.3, "b1", 0.1, 0.0, x=0.02)
    f1 = field(tao, ele, "ref_tilt", 0.3, "b1", 0.1, 0.0, settings=("x_offset = 0.01",), x=0.03)
    scale = math.hypot(f0[0], f0[1])
    assert scale > 0
    assert f1 == pytest.approx(f0, abs=1e-9 * scale)
