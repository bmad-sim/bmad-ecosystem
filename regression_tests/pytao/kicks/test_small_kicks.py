"""
Check that small kicks and multipole strengths give the expected first-order change in momentum.

For each element type that accepts kick or multipole attributes, a particle with a transverse
offset is tracked with and without a small perturbation of one attribute.  The change in exit
momentum is compared with the thin-lens expectation from the Bmad manual:

    magnetic:  (q L / P0) (By + i Bx) = (b_n + i a_n) (x + i y)^n
    electric:  Ex - i Ey = (b_en - i a_en) (x + i y)^n,  dp/ds = q E / (beta P0 c)

hkick and vkick are lab-frame kicks except on kicker and elseparator elements, where they rotate
with the tilt.  Multipoles are in the body frame; a tilt T multiplies (b_n + i a_n) by exp(-i (n+1) T).

Each element sits alone in its own lattice so nothing downstream alters the measured momentum, and
each element type is tested with every tracking method it accepts.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import pytest
from pytao import SubprocessTao
from pytao.errors import TaoCommandError

E_TOT = 1e9
M_ELECTRON = 0.51099895e6
P0C = math.sqrt(E_TOT**2 - M_ELECTRON**2)
BETA = P0C / E_TOT
CHARGE = -1

X0, Y0, PZ0 = 5e-3, 3e-3, 0.01
KICK = 1e-8  # Size of every expected kick. Small enough that second-order terms are below 1e-6.

STD = ["bmad_standard", "symp_lie_ptc", "runge_kutta", "time_runge_kutta", "linear", "taylor"]
STD_SLB = STD + ["symp_lie_bmad"]
ALL_METHODS = STD_SLB


@dataclass(frozen=True)
class Element:
    name: str
    definition: str
    length: float
    methods: list[str]
    tilts: tuple[float, ...] = (0.0, 0.3)
    body_kicks: bool = False  # hkick/vkick rotate with the tilt instead of staying in the lab frame.
    kick_sign: int = 1  # elseparator kicks are defined for a positive particle.
    rtol: float = 0.0  # Overrides the per-method tolerance when larger.


WIG = "b_max = 1e-4, num_steps = 100, ds_step = 0.005"
# By default an element's a_n/b_n are scaled by its main strength (k1 L, k2 L, ...), which is zero here.
NS = "scale_multipoles = F"

ELEMENTS = [
    Element("ecollimator", "ecollimator, l = 0.5", 0.5, STD),
    Element("rcollimator", "rcollimator, l = 0.5", 0.5, STD),
    Element("monitor", "monitor, l = 0.5", 0.5, STD),
    Element("instrument", "instrument, l = 0.5", 0.5, STD),
    Element("pipe", "pipe, l = 0.5", 0.5, STD),
    Element("quadrupole", f"quadrupole, l = 0.5, {NS}", 0.5, STD_SLB),
    Element("sextupole", f"sextupole, l = 0.5, {NS}", 0.5, STD),
    Element("octupole", f"octupole, l = 0.5, {NS}", 0.5, STD),
    Element("solenoid", f"solenoid, l = 0.5, {NS}", 0.5, STD_SLB),
    Element("sol_quad", f"sol_quad, l = 0.5, {NS}", 0.5, STD_SLB),
    Element("sbend", f"sbend, l = 0.5, {NS}", 0.5, STD, tilts=(0.0,)),
    Element("thick_multipole", "thick_multipole, l = 0.5", 0.5, STD),
    Element("kicker", f"kicker, l = 0.5, {NS}", 0.5, STD, body_kicks=True),
    Element("hkicker", f"hkicker, l = 0.5, {NS}", 0.5, STD, body_kicks=True),
    Element("vkicker", f"vkicker, l = 0.5, {NS}", 0.5, STD, body_kicks=True),
    Element("elseparator", f"elseparator, l = 0.5, {NS}", 0.5, STD, body_kicks=True, kick_sign=CHARGE),
    Element("rfcavity", "rfcavity, l = 0.5, rf_frequency = 500e6, voltage = 0", 0.5, STD),
    Element("lcavity", "lcavity, l = 0.5, rf_frequency = 500e6, voltage = 0", 0.5, STD),
    Element("lcavity_old", "lcavity, l = 0.5, rf_frequency = 500e6, voltage = 0, n_rf_steps = 0", 0.5, STD),
    Element("crab_cavity", f"crab_cavity, l = 0.5, rf_frequency = 500e6, voltage = 0, {NS}", 0.5,
            ["bmad_standard", "symp_lie_ptc", "runge_kutta", "time_runge_kutta", "linear"]),
    Element("wiggler_planar", f"wiggler, l = 0.5, l_period = 0.1, field_calc = planar_model, {WIG}, {NS}", 0.5, STD_SLB,
            rtol=1e-4),
    Element("wiggler_helical", f"wiggler, l = 0.5, n_period = 5, field_calc = helical_model, {WIG}, {NS}", 0.5,
            ["bmad_standard", "runge_kutta", "time_runge_kutta", "linear", "symp_lie_bmad"], rtol=1e-4),
    Element("quadrupole_zero", f"quadrupole, l = 0, {NS}", 0.0, STD),
    Element("sextupole_zero", f"sextupole, l = 0, {NS}", 0.0, STD),
    Element("kicker_zero", f"kicker, l = 0, {NS}", 0.0, STD, body_kicks=True),
    Element("rfcavity_zero", "rfcavity, l = 0, rf_frequency = 500e6, voltage = 0", 0.0, STD),
    Element("lcavity_zero", "lcavity, l = 0, rf_frequency = 500e6, voltage = 0, cavity_type = traveling_wave", 0.0, STD),
]


@dataclass(frozen=True)
class Attribute:
    name: str
    kind: str  # kick_x, kick_y, mag_b, mag_a, elec_b, elec_a
    n: int = 0
    delta: float = KICK
    to_pole: Callable[[float, float], float] = lambda d, L: d  # (delta, length) -> b_n or a_n


R = abs(complex(X0, Y0))
L0 = 0.5
EF = L0 / (BETA * P0C)  # Electric kick per (V/m) over the standard length.

ATTRIBUTES = [
    Attribute("hkick", "kick_x"),
    Attribute("vkick", "kick_y"),
    Attribute("kick", "kick_x"),
    Attribute("b0", "mag_b", 0),
    Attribute("a0", "mag_a", 0),
    Attribute("b1", "mag_b", 1, KICK / R),
    Attribute("a1", "mag_a", 1, KICK / R),
    Attribute("b2", "mag_b", 2, KICK / R**2),
    Attribute("a2", "mag_a", 2, KICK / R**2),
    Attribute("b3", "mag_b", 3, KICK / R**3),
    Attribute("a3", "mag_a", 3, KICK / R**3),
    Attribute("k1", "mag_b", 1, KICK / R / L0, lambda d, L: d * L),
    Attribute("k2", "mag_b", 2, 2 * KICK / R**2 / L0, lambda d, L: d * L / 2),
    Attribute("k3", "mag_b", 3, 6 * KICK / R**3 / L0, lambda d, L: d * L / 6),
    Attribute("dg", "mag_b", 0, KICK / L0, lambda d, L: d * L),
    Attribute("b0_elec", "elec_b", 0, KICK / EF),
    Attribute("a0_elec", "elec_a", 0, KICK / EF),
    Attribute("b1_elec", "elec_b", 1, KICK / EF / R),
    Attribute("a1_elec", "elec_a", 1, KICK / EF / R),
]

TOLERANCE = {
    "bmad_standard": 1e-5,
    "symp_lie_ptc": 1e-5,
    "linear": 1e-5,
    "taylor": 1e-5,
    "symp_lie_bmad": 1e-5,
    "runge_kutta": 1e-4,
    "time_runge_kutta": 1e-4,
}

PTC = ("symp_lie_ptc", "taylor")
RK = ("runge_kutta", "time_runge_kutta")
ALL = tuple(ALL_METHODS)
MAG = ("a0", "b0", "a1", "b1", "a2", "b2", "a3", "b3")
ELEC = ("a0_elec", "b0_elec", "a1_elec", "b1_elec")
HV = ("hkick", "vkick")
CAVITIES = ("rfcavity", "lcavity", "lcavity_old")

# (elements, methods, attributes, reason). Any of these kills the Tao process, so they are skipped
# rather than marked xfail. None matches everything.
CRASHES: list = []

# (elements, methods, attributes, reason[, tilts]). Known wrong answers, marked xfail(strict) so a fix
# is noticed. None matches everything.
KNOWN = [
    (("kicker", "hkicker", "vkicker"), PTC, ("b0", "a1", "b1", "a2", "b2", "a3", "b3"),
     "PTC kicker: b0 has the wrong sign and orders 1 to 3 are applied at half strength"),
    (("kicker_zero",), PTC, ("b0",), "PTC kicker: b0 has the wrong sign"),
    (("elseparator",), PTC, ("b0",), "PTC: multipoles in an elseparator are not supported"),
    (CAVITIES + ("rfcavity_zero",), PTC, HV, "PTC cavity applies hkick/vkick over l_active instead of l"),
    (("crab_cavity",), PTC, MAG + HV, "PTC crab_cavity ignores kicks and multipoles"),
    (CAVITIES, ALL, ("a1_elec", "b1_elec"), "electric quadrupole in a cavity gives no kick"),
    (CAVITIES, RK + PTC, ("a0_elec", "b0_elec"), "electric dipole in a cavity is applied over l_active instead of l"),
    (("solenoid",), ("bmad_standard", "linear", "symp_lie_bmad"), ("b1",), "solenoid drops its b1 multipole"),
    (("quadrupole_zero",), ("bmad_standard", "linear"), ("b1",), "zero length quadrupole drops its b1 multipole"),
    (("quadrupole_zero", "sextupole_zero"), PTC, MAG + HV, "PTC zero length element ignores kicks and multipoles"),
    (("wiggler_planar", "wiggler_helical"), RK, MAG + ELEC, "em_field_calc periodic wiggler field has no multipoles"),
    (("wiggler_planar",), PTC, MAG + HV, "PTC periodic wiggler ignores kicks and multipoles"),
    (("wiggler_planar",), PTC, ELEC, "PTC wiggler cannot have electric multipoles so they are ignored"),
    (("lcavity_zero",), PTC, MAG + HV, "PTC zero length lcavity ignores kicks and multipoles"),
    (("wiggler_planar", "wiggler_helical"), ("symp_lie_bmad",), ("b1",), "symp_lie_bmad wiggler drops its b1 multipole"),
]


def lookup(table, ele: Element, method: str, attr: Attribute, tilt: float) -> str | None:
    for elements, methods, attrs, reason, *tilts in table:
        if ele.name in elements and method in methods and (attrs is None or attr.name in attrs):
            if not tilts or tilt in tilts[0]:
                return reason
    return None


def expected_kick(attr: Attribute, ele: Element, tilt: float) -> tuple[float, float]:
    if attr.kind in ("kick_x", "kick_y"):
        kx, ky = (attr.delta, 0.0) if attr.kind == "kick_x" else (0.0, attr.delta)
        kx, ky = ele.kick_sign * kx, ele.kick_sign * ky
        if not ele.body_kicks:
            return kx, ky
        c, s = math.cos(tilt), math.sin(tilt)
        return kx * c - ky * s, kx * s + ky * c

    z = complex(X0, Y0)
    pole = attr.to_pole(attr.delta, ele.length)
    phase = cmath.exp(-1j * (attr.n + 1) * tilt)
    if attr.kind in ("mag_b", "mag_a"):
        w = (pole if attr.kind == "mag_b" else 1j * pole) * phase * z**attr.n
        return -w.real, w.imag

    w = (pole if attr.kind == "elec_b" else -1j * pole) * phase * z**attr.n
    f = CHARGE * ele.length / (BETA * P0C)
    return f * w.real, -f * w.imag


LATTICE = """\
parameter[particle] = electron
parameter[geometry] = open
beginning[e_tot] = {e_tot}
beginning[beta_a] = 10
beginning[beta_b] = 10
bmad_com[rel_tol_adaptive_tracking] = 1e-10
bmad_com[abs_tol_adaptive_tracking] = 1e-13
ELE: {definition}
lat: line = (ELE)
use, lat
"""


class TaoElement:
    def __init__(self, ele: Element, lattice_file: Path):
        self.ele = ele
        self.tao = SubprocessTao(lattice_file=str(lattice_file), noplot=True)
        self.tao.cmd(f"set particle_start x = {X0}")
        self.tao.cmd(f"set particle_start y = {Y0}")
        self.tao.cmd(f"set particle_start pz = {PZ0}")
        self.state: tuple[str, float] | None = None
        self.base: dict[tuple[str, float], tuple[float, float]] = {}
        self.has: dict[str, bool] = {}

    def has_attribute(self, name: str) -> bool:
        if name not in self.has:
            try:
                self.tao.cmd(f"set ele ELE {name} = 0")
            except TaoCommandError:
                self.has[name] = False
            else:
                self.has[name] = True
        return self.has[name]

    def momentum(self) -> tuple[float, float]:
        orb = self.tao.ele_orbit("ELE")
        return orb["px"], orb["py"]

    def configure(self, method: str, tilt: float) -> tuple[float, float]:
        if self.state != (method, tilt):
            self.tao.cmd(f"set ele ELE tracking_method = {method}")
            if tilt != 0 or self.ele.tilts != (0.0,):
                self.tao.cmd(f"set ele ELE tilt = {tilt}")
            self.state = (method, tilt)
        if (method, tilt) not in self.base:
            self.base[(method, tilt)] = self.momentum()
        return self.base[(method, tilt)]

    def perturbed(self, name: str, value: float) -> tuple[float, float]:
        self.tao.cmd(f"set ele ELE {name} = {value}")
        try:
            return self.momentum()
        finally:
            self.tao.cmd(f"set ele ELE {name} = 0")


@pytest.fixture(scope="module", params=ELEMENTS, ids=lambda e: e.name)
def tao_element(request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory):
    ele: Element = request.param
    lattice_file = tmp_path_factory.mktemp("kicks") / f"{ele.name}.bmad"
    lattice_file.write_text(LATTICE.format(e_tot=E_TOT, definition=ele.definition))
    te = TaoElement(ele, lattice_file)
    yield te
    te.tao.close_subprocess()


@pytest.mark.parametrize("tilt", [0.0, 0.3])
@pytest.mark.parametrize("method", ALL_METHODS)
@pytest.mark.parametrize("attr", ATTRIBUTES, ids=lambda a: a.name)
def test_small_kick(request: pytest.FixtureRequest, tao_element: TaoElement, method: str, tilt: float, attr: Attribute):
    ele = tao_element.ele
    if method not in ele.methods:
        pytest.skip(f"{method} is not valid for {ele.name}")
    if tilt not in ele.tilts:
        pytest.skip(f"{ele.name} has no tilt")
    if ele.length == 0 and attr.kind.startswith("elec"):
        pytest.skip("an electric field over zero length gives no kick")
    crash = lookup(CRASHES, ele, method, attr, tilt)
    if crash:
        pytest.skip(crash)
    if not tao_element.has_attribute(attr.name):
        pytest.skip(f"{ele.name} has no {attr.name}")
    if method == "linear" and attr.n > 1:
        pytest.skip("linear tracking cannot represent a nonlinear kick at an offset")
    if attr.name == "kick" and ele.name == "vkicker":
        attr = Attribute("kick", "kick_y")

    want_x, want_y = expected_kick(attr, ele, tilt)
    scale = max(abs(want_x), abs(want_y))
    if scale == 0:
        pytest.skip("no kick expected")

    known = lookup(KNOWN, ele, method, attr, tilt)
    if known:
        request.applymarker(pytest.mark.xfail(reason=known, strict=True))

    px0, py0 = tao_element.configure(method, tilt)
    px1, py1 = tao_element.perturbed(attr.name, attr.delta)
    got_x, got_y = px1 - px0, py1 - py0

    tol = max(TOLERANCE[method], ele.rtol) * scale
    assert got_x == pytest.approx(want_x, abs=tol) and got_y == pytest.approx(want_y, abs=tol), (
        f"{ele.name} {method} tilt={tilt} {attr.name}={attr.delta:.4g}: "
        f"dpx={got_x:.6e} (expected {want_x:.6e}), dpy={got_y:.6e} (expected {want_y:.6e})"
    )
