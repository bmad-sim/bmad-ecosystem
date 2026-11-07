"""
Tests of the "set data" and "set var" commands and their agreement with the
settable flags reported by "pipe data" and "pipe var".
"""

import shutil
from pathlib import Path

import pytest
from pytao import SubprocessTao, TaoCommandError

ROOT = Path(__file__).resolve().parent
INPUT_FILES = ("set_data_var.bmad", "set_data_var.init")

DATUM = "1@twiss.beta[1]"
VAR = "quad[1]"


@pytest.fixture
def tao(tmp_path, monkeypatch):
    for name in INPUT_FILES:
        shutil.copy(ROOT / name, tmp_path / name)
    monkeypatch.chdir(tmp_path)
    with SubprocessTao(init_file="set_data_var.init", noplot=True) as tao:
        yield tao


def pipe_params(tao, what: str, who: str) -> dict[str, tuple[str, str]]:
    """Return {name: (settable_flag, value)} from "pipe data" or "pipe var"."""
    params = {}
    for line in tao.cmd(f"pipe {what} {who}"):
        name, _type, settable, value = line.split(";", 3)
        # Universe prefixed names like "1^ix_bunch" and type prefixed names like
        # "data^merit_type" are set using the bare name.
        params[name.split("^")[-1]] = (settable, value.strip())
    return params


def set_ok(tao, cmd: str) -> bool:
    try:
        tao.cmd(cmd, raises=True)
    except TaoCommandError:
        return False
    return True


@pytest.mark.parametrize("what, who", [("data", DATUM), ("var", VAR)])
def test_settable_flags_match_set_command(tao, what, who):
    """
    Every parameter flagged settable (T) must be accepted by "set", and every
    parameter flagged not settable (F) must be refused.
    Each parameter is set to its current value so the state is not changed.
    """
    wrong = []
    for name, (settable, value) in pipe_params(tao, what, who).items():
        if settable not in ("T", "F") or value == "":
            continue
        accepted = set_ok(tao, f"set {what} {who}|{name} = {value}")
        if accepted != (settable == "T"):
            wrong.append(f"{name}: flag={settable} but set {'accepted' if accepted else 'refused'}")
    assert not wrong, "\n".join(wrong)


@pytest.mark.parametrize("component", ["ele_name", "ele_start_name", "ele_ref_name"])
def test_failed_ele_name_set_leaves_datum_unchanged(tao, component):
    """
    "qd" matches several elements so the set must fail, and the datum must be unchanged.
    """
    before = pipe_params(tao, "data", DATUM)
    assert not set_ok(tao, f"set data {DATUM}|{component} = qd")
    assert pipe_params(tao, "data", DATUM) == before


def test_failed_ix_ele_set_leaves_data_unchanged(tao):
    """An out of range element index set on several datums must not modify any of them."""
    who = "1@twiss.beta[1:2]"
    before = [pipe_params(tao, "data", f"1@twiss.beta[{i}]") for i in (1, 2)]
    assert not set_ok(tao, f"set data {who}|ix_ele = 9999")
    after = [pipe_params(tao, "data", f"1@twiss.beta[{i}]") for i in (1, 2)]
    assert after == before


def test_ele_name_set(tao):
    """A valid element name set updates both the name and the element index."""
    assert set_ok(tao, f"set data {DATUM}|ele_name = 3")
    params = pipe_params(tao, "data", DATUM)
    assert params["ele_name"][1] == "B"
    assert params["ix_ele"][1] == "3"
    assert params["exists"][1] == "T"
