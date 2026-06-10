"""Regression tests for `pipe ele:floor` on girder_lord elements."""

from pathlib import Path

import numpy as np
import pytest
from pytao import SubprocessTao


LAT_PATH = Path(__file__).parent / "lat.bmad"
LAT_MULTIPASS_PATH = Path(__file__).parent / "lat_multipass.bmad"


@pytest.fixture(scope="module")
def tao():
    assert LAT_PATH.is_file(), f"Lattice file not found: {LAT_PATH}"
    with SubprocessTao(lattice_file=str(LAT_PATH), noplot=True) as t:
        yield t


@pytest.fixture(scope="module")
def tao_multipass():
    assert LAT_MULTIPASS_PATH.is_file(), f"Lattice file not found: {LAT_MULTIPASS_PATH}"
    with SubprocessTao(lattice_file=str(LAT_MULTIPASS_PATH), noplot=True) as t:
        yield t


def test_girder_floor_beginning_does_not_crash(tao):
    """Regression test: `pipe ele:floor <girder> beginning` used to segfault.

    For a girder_lord the previous implementation called
    `pointer_to_next_ele(ele, -1)`, which returns null for any lord that is
    not a super_lord, and then dereferenced the result.
    """
    result = tao.ele_floor("g1", where="beginning")
    # Floor position of the upstream end of the first slave (p1) is the
    # lattice origin.
    np.testing.assert_allclose(result["Reference"], [0, 0, 0, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(result["Actual"], [0, 0, 0, 0, 0, 0], atol=1e-12)


def test_girder_floor_end_matches_last_slave(tao):
    """`pipe ele:floor <girder> end` returns the downstream end of the last slave."""
    girder_end = tao.ele_floor("g1", where="end")
    last_slave_end = tao.ele_floor("p2", where="end")
    np.testing.assert_allclose(girder_end["Reference"], last_slave_end["Reference"], atol=1e-12)
    np.testing.assert_allclose(girder_end["Actual"], last_slave_end["Actual"], atol=1e-12)
    # Total length p1+p2 = 3 m, lattice is along z.
    np.testing.assert_allclose(girder_end["Reference"], [0, 0, 3, 0, 0, 0], atol=1e-12)


def test_girder_floor_beginning_matches_first_slave(tao):
    """`pipe ele:floor <girder> beginning` returns the upstream end of the first slave."""
    girder_begin = tao.ele_floor("g1", where="beginning")
    first_slave_begin = tao.ele_floor("p1", where="beginning")
    np.testing.assert_allclose(
        girder_begin["Reference"], first_slave_begin["Reference"], atol=1e-12
    )
    np.testing.assert_allclose(girder_begin["Actual"], first_slave_begin["Actual"], atol=1e-12)


def test_girder_floor_center(tao):
    """`pipe ele:floor <girder> center` returns the girder's own floor position."""
    result = tao.ele_floor("g1", where="center")
    # Center of a girder is the midpoint of the slave span: z = (0 + 3) / 2 = 1.5.
    np.testing.assert_allclose(result["Reference"], [0, 0, 1.5, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(result["Actual"], [0, 0, 1.5, 0, 0, 0], atol=1e-12)


def test_girder_floor_multipass_does_not_crash(tao_multipass):
    """Regression test: `pipe ele:floor <girder> beginning` used to segfault
    when the girder's slaves are multipass_lord elements.

    pointer_to_slave(girder, 1) returns the multipass_lord (which is not in
    the tracking part of the lattice), so pointer_to_next_ele returned null.
    The fix descends through nested lord chains via find_element_ends.
    """
    begin = tao_multipass.ele_floor("g1", where="beginning")
    end = tao_multipass.ele_floor("g1", where="end")
    # First pass of the multipass line is at z in [0, 3].
    np.testing.assert_allclose(begin["Reference"], [0, 0, 0, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(end["Reference"], [0, 0, 3, 0, 0, 0], atol=1e-12)
