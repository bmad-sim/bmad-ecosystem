"""
Tests of the "change" command: switch placement, the "@", "d", and "%" prefixes
written with and without a space, and the precision of a percent change.
Also tests that the single mode saved value "key_val0" of a variable is settable.
"""

import shutil
from pathlib import Path

import pytest
from pytao import SubprocessTao, TaoCommandError

ROOT = Path(__file__).resolve().parent
INPUT_FILES = ("change_cmd.bmad", "change_cmd.init")

VAR = "quad[1]"  # Varies the K1 of all the QF quads.
DESIGN_K1 = 0.3


@pytest.fixture
def tao(tmp_path, monkeypatch):
    for name in INPUT_FILES:
        shutil.copy(ROOT / name, tmp_path / name)
    monkeypatch.chdir(tmp_path)
    with SubprocessTao(init_file="change_cmd.init", noplot=True) as tao:
        yield tao


def qf_k1(tao) -> float:
    return tao.ele_gen_attribs("1")["K1"]  # Element 1 is the first QF.


def var_param(tao, name: str) -> tuple[str, str]:
    """Return (settable_flag, value) of a "pipe var" parameter."""
    for line in tao.cmd(f"pipe var {VAR}"):
        pname, _type, settable, value = line.split(";", 3)
        if pname == name:
            return settable, value.strip()
    raise KeyError(name)


def set_k1(tao, k1: float):
    tao.cmd(f"set ele qf k1 = {k1}")


# Switches before and after the subcommand word. The manual puts them before.


@pytest.mark.parametrize(
    "cmd", ["change -update element qf k1 0.01", "change element -update qf k1 0.01"]
)
def test_update_switch(tao, cmd):
    set_k1(tao, 0.3)
    tao.cmd(cmd)
    assert qf_k1(tao) == pytest.approx(0.31, abs=1e-14)


@pytest.mark.parametrize(
    "cmd", [f"change -silent variable {VAR} 0.01", f"change variable -silent {VAR} 0.01"]
)
def test_silent_switch(tao, cmd):
    set_k1(tao, 0.3)
    out = "\n".join(tao.cmd(cmd))
    assert "Old_Model" not in out
    assert qf_k1(tao) == pytest.approx(0.31, abs=1e-14)


@pytest.mark.parametrize(
    "cmd", ["change -listing tune 0 0.01", "change tune -listing 0 0.01"]
)
def test_listing_switch(tao, cmd):
    out = tao.cmd(cmd)
    assert any(line.split()[:1] == ["QF"] for line in out)
    assert qf_k1(tao) != pytest.approx(DESIGN_K1, abs=1e-6)


@pytest.mark.parametrize(
    "cmd", ["change -mask Q* tune 0 0.01", "change tune -mask Q* 0 0.01"]
)
def test_mask_switch(tao, cmd):
    # The mask is a veto list so vetoing all the quads must give an error.
    with pytest.raises(TaoCommandError, match="CANNOT FIND A QUAD"):
        tao.cmd(cmd)


@pytest.mark.parametrize(
    "cmd", ["change -branch 0 z_tune 0.0", "change z_tune -branch 0 0.0"]
)
def test_branch_switch(tao, cmd):
    out = "\n".join(tao.cmd(cmd))
    assert "Change who" not in out


@pytest.mark.parametrize(
    "cmd",
    [
        "change -update -silent element qf k1 0.01",
        "change -silent element -update qf k1 0.01",
    ],
)
def test_two_switches(tao, cmd):
    set_k1(tao, 0.3)
    tao.cmd(cmd)
    assert qf_k1(tao) == pytest.approx(0.31, abs=1e-14)


# Prefixes "@", "d", "%" with and without a space. The start K1 of 0.25 differs from
# the design value so that a change that is ignored is caught. The tolerance is tight
# enough to catch a single precision percent change (which is off by ~2e-9 here).

START_K1 = 0.25
CHANGES = [
    ("0.01", 0.26),
    ("-0.01", 0.24),
    ("@0.4", 0.4),
    ("@ 0.4", 0.4),
    ("d0.02", DESIGN_K1 + 0.02),
    ("d 0.02", DESIGN_K1 + 0.02),
    ("d -0.01", DESIGN_K1 - 0.01),
    ("%50", 0.375),
    ("% 50", 0.375),
    ("%-10", 0.225),
    ("% -10", 0.225),
]


@pytest.mark.parametrize("change, k1", CHANGES)
def test_change_element_prefix(tao, change, k1):
    set_k1(tao, START_K1)
    tao.cmd(f"change element qf k1 {change}")
    assert qf_k1(tao) == pytest.approx(k1, abs=1e-14)


@pytest.mark.parametrize("change, k1", CHANGES)
def test_change_variable_prefix(tao, change, k1):
    set_k1(tao, START_K1)
    tao.cmd(f"change variable {VAR} {change}")
    assert float(var_param(tao, "model_value")[1]) == pytest.approx(k1, abs=1e-14)
    assert qf_k1(tao) == pytest.approx(k1, abs=1e-14)


def test_percent_change_is_double_precision(tao):
    """The reported case: 0.32 changed by 50% must be 0.48 and not 0.4799999964."""
    set_k1(tao, 0.32)
    tao.cmd("change element qf K1 %50")
    assert qf_k1(tao) == pytest.approx(0.48, abs=1e-14)


def test_key_val0_settable(tao):
    tao.cmd(f"set var {VAR}|key_val0 = 0.123")
    settable, value = var_param(tao, "key_val0")
    assert settable == "T"
    assert float(value) == pytest.approx(0.123, abs=1e-15)


# A switch in the last of the 8 command words (or "-mask <list>" as the last two) must not
# be found again after the switch removal shifts the words down. A hang here fails via the timeout.


@pytest.mark.parametrize(
    "cmd",
    [
        "change element qf k1 0.01 -update -update -update -silent",
        "change -silent -update -update element qf k1 0.01 -update",
        "change tune 0 0.01 -listing -listing -listing -mask XX*",
        "change tune 0 0.01 -listing -listing -listing -branch 0",
    ],
)
def test_switch_in_last_word_terminates(tao, cmd):
    set_k1(tao, 0.3)
    with tao.timeout(30):
        tao.cmd(cmd)
    if "element" in cmd:
        assert qf_k1(tao) == pytest.approx(0.31, abs=1e-14)
    else:
        assert qf_k1(tao) != pytest.approx(0.3, abs=1e-6)
