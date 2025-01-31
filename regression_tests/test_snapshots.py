from __future__ import annotations
import numpy as np
import pathlib
import pytest
from pytao import Tao, SubprocessTao

from typing import NamedTuple

import conftest

SNAPSHOTS = conftest.TESTS_ROOT / "snapshots"


class Tolerance(NamedTuple):
    atol: float
    rtol: float


class SnapshotPaths(NamedTuple):
    lat_list: pathlib.Path


examples = pytest.mark.parametrize(
    ("example_name",),
    [
        pytest.param("cesr"),
        # pytest.param("cbeta_cell"),
        # pytest.param("cbeta_ffag"),
        # pytest.param("csr_beam_tracking"),
        # pytest.param("custom_tao_with_measured_data"),
        # pytest.param("driving_terms"),
        # pytest.param("dynamic_aperture"),
        # pytest.param("erl"),
        # pytest.param("fodo"),
        # pytest.param("multi_turn_orbit"),
        # pytest.param("optics_matching"),
        # pytest.param("space_charge"),
        # pytest.param("x_axis_param_plot"),
    ],
)


lat_list_to_tolerance = {
    "ele.a.beta": Tolerance(rtol=1e-7, atol=0),
    "ele.b.beta": Tolerance(rtol=1e-7, atol=0),
    # TODO what else?
}


def load_example(example_name: str) -> Tao:
    return SubprocessTao(
        init_file=conftest.TAO_EXAMPLES_ROOT / example_name / "tao.init",
        noplot=True,
    )


def get_snapshot_paths(example_name: str) -> SnapshotPaths:
    example_snapshots = SNAPSHOTS / example_name
    example_snapshots.mkdir(parents=True, exist_ok=True)

    return SnapshotPaths(
        lat_list=example_snapshots / "lat_list.npz",
    )


def snapshot_lat_list(tao: Tao) -> dict[str, np.ndarray]:
    return {
        who: tao.lat_list(elements="*", who=who, flags="-array_out -track_only")
        for who in lat_list_to_tolerance
    }


@examples
def test_lat_list_against_snapshot(example_name: str) -> None:
    tao = load_example(example_name)
    paths = get_snapshot_paths(example_name)

    expected_lat_list = np.load(paths.lat_list)
    lat_list = snapshot_lat_list(tao)

    for key, current_value in lat_list.items():
        tol = lat_list_to_tolerance[key]
        print(f"Checking {key} with tolerance {tol}")
        np.testing.assert_allclose(
            current_value,
            expected_lat_list[key],
            atol=tol.atol,
            rtol=tol.rtol,
        )


def update_snapshots() -> None:
    example_name: str

    for param_set in examples.args[1]:
        example_name, *_ = param_set.values
        print(f"Updating snapshot for: {example_name}")

        paths = get_snapshot_paths(example_name)

        tao = load_example(example_name)
        lat_list = snapshot_lat_list(tao)

        np.savez(paths.lat_list, **lat_list)


if __name__ == "__main__":
    update_snapshots()
