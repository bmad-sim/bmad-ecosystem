from __future__ import annotations

import contextlib
import os
import pprint
from dataclasses import dataclass
from pathlib import Path

import pytest
from pytao import TaoStartup

INIT_ROOT = Path(__file__).resolve().parent


@pytest.fixture(params=[1, 10, 100, 250])
def path_prefix_length(request: pytest.FixtureRequest):
    return request.param


@pytest.fixture
def path_prefix(path_prefix_length: int):
    return "x" * path_prefix_length


@dataclass
class TaoInitLog:
    init_file: list[str] | None = None
    lattice_info_file: list[str] | None = None
    bmad_file: list[str] | None = None
    beam_file: list[str] | None = None
    variable_file: list[str] | None = None
    data_file: list[str] | None = None
    plotting_file: list[str] | None = None
    other_plotting_file: list[str] | None = None

    @classmethod
    def from_output(cls, output: str):
        prefix_to_attr = {
            "*Init: Opening Init File: ": "init_file",
            "*Init: Opening Lattice Info File: ": "lattice_info_file",
            "Reading Bmad file: ": "bmad_file",
            "*Init: Opening Beam File: ": "beam_file",
            "*Init: Opening Variable File: ": "variable_file",
            "*Init: Opening Data File: ": "data_file",
            "*Init: Opening Plotting File: ": "plotting_file",
            "Init: Opening another plotting file: ": "other_plotting_file",
        }

        values = {}

        for line in output.splitlines():
            for prefix, attr in prefix_to_attr.items():
                if line.startswith(prefix):
                    value = line.removeprefix(prefix).strip()
                    values[attr] = [*values.get(attr, ()), value]

                    break
        return cls(**values)


@contextlib.contextmanager
def init_with_log(workdir: Path, startup: TaoStartup):
    prev_wd = Path.cwd()
    init_fn = workdir / "tao_init.log"
    startup.log_startup = True
    startup.noplot = True

    for fn in (
        startup.beam_file,
        startup.beam_init_position_file,
        startup.building_wall_file,
        startup.data_file,
        startup.hook_init_file,
        startup.init_file,
        startup.lattice_file,
        startup.plot_file,
        startup.startup_file,
        startup.var_file,
    ):
        if fn:
            Path(workdir / fn).touch()

    try:
        os.chdir(workdir)
        with startup.run(use_subprocess=True) as tao:
            yield tao, TaoInitLog.from_output(init_fn.read_text())
    finally:
        os.chdir(prev_wd)


def test_can_initialize(tmp_path: Path, path_prefix: str):
    startup = TaoStartup(lattice_file=INIT_ROOT / "lat1.bmad")
    startup.beam_file = path_prefix  # relative to the tao workdir
    startup.beam_init_position_file = startup.beam_file
    startup.building_wall_file = startup.beam_file
    startup.data_file = f"{path_prefix}_data"
    startup.hook_init_file = startup.beam_file
    startup.init_file = startup.beam_file
    startup.plot_file = startup.beam_file
    startup.startup_file = startup.beam_file
    startup.var_file = startup.beam_file
    with init_with_log(tmp_path, startup) as (_tao, log):
        assert log.beam_file == [startup.beam_file]
        assert log.data_file == [startup.data_file]
        assert log.bmad_file == [str(startup.lattice_file)]
        pprint.pprint(log)
