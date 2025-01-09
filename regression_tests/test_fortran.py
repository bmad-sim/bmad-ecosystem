from __future__ import annotations
import os
import shlex
import math
import pathlib
import pytest
import sys
from typing import NamedTuple

from conftest import TESTS_ROOT

# List all Fortran tests here.
# If a test should be skipped, add the `marks` kwarg to the `pytest.param` as follows:
#    marks=pytest.mark.skip(reason="Insert reason for the skip here")
# If a test is expected to fail, use:
#    marks=pytest.mark.xfail
fortran_tests = pytest.mark.parametrize(
    ("test_name",),
    [
        pytest.param("abs_time_test"),
        pytest.param("analysis_test"),
        pytest.param("aperture_test"),
        pytest.param("autoscale_test"),
        pytest.param("backwards_time_track_test"),
        pytest.param("bbu_test"),
        pytest.param("beam_test"),
        pytest.param("bookkeeper_test"),
        pytest.param(
            "cathode_sc_test", marks=pytest.mark.skip(reason="TODO reason here")
        ),
        pytest.param("cesr_test"),
        pytest.param("closed_orbit_test"),
        pytest.param("coord_test"),
        pytest.param("csr_and_space_charge_test"),
        pytest.param("envelope_ibs_test"),
        pytest.param("fringe_test"),
        pytest.param("geometry_test"),
        pytest.param("girder_test"),
        pytest.param("hdf5_test"),
        pytest.param("long_term_tracking_test"),
        pytest.param("mat6_calc_method_test"),
        pytest.param("match_test"),
        pytest.param("mode3_test"),
        pytest.param("mode_flip_test"),
        pytest.param("multipass_test"),
        pytest.param("multipole_test"),
        pytest.param("nonlin_test"),
        pytest.param("parse_test"),
        pytest.param("particle_species_test"),
        pytest.param("patch_test"),
        pytest.param("photon_test"),
        pytest.param("ptc_test"),
        pytest.param("radiation_test"),
        pytest.param("reverse_test", marks=pytest.mark.skip(reason="TODO reason here")),
        pytest.param("sad_test"),
        pytest.param("sim_utils_test"),
        pytest.param("slice_test"),
        pytest.param("space_charge_test"),
        pytest.param("spin_general_test"),
        pytest.param("superimpose_test"),
        pytest.param("synrad3d_test"),
        pytest.param("tao_test"),
        pytest.param("taylor_test"),
        pytest.param("time_runge_kutta_test"),
        pytest.param("tracking_method_test"),
        pytest.param(
            "twiss_test",
            marks=pytest.mark.skip(reason="TODO reason here"),
        ),
        pytest.param("wake_test"),
        pytest.param("wall3d_test"),
        pytest.param("write_bmad_test"),
        pytest.param("write_foreign_test"),
        pytest.param("xraylib_test"),
    ],
)


class Line(NamedTuple):
    lineno: int

    name: str
    type: str  # type: Literal["STR", "ABS", "REL", "VEC_REL"]
    tolerance: float | None
    values: list[str] | list[float]

    def __str__(self) -> str:
        return f"[line {self.lineno}] {self.name} type={self.type} tolerance={self.tolerance} values={self.values}"

    @classmethod
    def from_string(cls, line: str, lineno: int = 0) -> Line:
        """
        Create a Line object from a string.

        Parameters
        ----------
        line : str
            The line of text to parse.
        lineno : int, optional
            The line number in the source file, by default 0.

        Returns
        -------
        Line

        Raises
        ------
        ValueError
            If the line is blank or a comment line.
        """
        if not line:
            raise ValueError("Blank line")
        elif line[0] == "!":
            raise ValueError("Comment line")
        name, type_, *values = shlex.split(line)

        if type_ == "STR":
            tolerance = None
        else:
            values = [float(value) for value in values]
            tolerance = values.pop(0)

        return Line(
            name=name,
            type=type_,
            values=values,
            tolerance=tolerance,
            lineno=lineno,
        )


def load_output_file(fn: pathlib.Path) -> list[Line]:
    """
    Load and parse a file to extract relevant lines.

    This function reads a file specified by `fn`, processes each line, and
    returns a list of `Line` objects. It skips blank lines and those prefixed
    with '!'. If a line cannot be parsed into a `Line` object, an exception is
    raised.

    Parameters
    ----------
    fn : pathlib.Path
        The path to the file that needs to be loaded and parsed.

    Returns
    -------
    list of Line

    Raises
    ------
    ValueError
        If a line cannot be parsed into a `Line` object.
    """
    lines = []
    with open(fn) as fp:
        for lineno, line in enumerate(fp.read().splitlines(), start=1):
            line = line.strip()
            if not line:
                continue  # Skip blank line
            elif line[0] == "!":
                continue  # Skip comment line
            try:
                lines.append(Line.from_string(line, lineno=lineno))
            except Exception:
                raise ValueError(f"Cannot parse line from '{fn}':{lineno} {line!r}")
    return lines


def run_fortran_code(
    bmad_bin_dir: pathlib.Path,
    test_path: pathlib.Path,
):
    output_fn = test_path / "output.now"
    output_fn.unlink(missing_ok=True)

    program = test_path.name

    run_py = test_path / "run.py"
    if run_py.exists():
        print("     Found run.py. Running this script with python3.")
        if os.system(f"{sys.executable} run.py {bmad_bin_dir}/") != 0:
            raise RuntimeError("`run.py` returned a non-zero exit code")

    else:
        program = bmad_bin_dir / program
        print(f"     Running program: {program}")

        if not program.exists():
            raise FileNotFoundError(f"Test binary does not exist: {program}")

        if os.system(str(program)) != 0:
            raise RuntimeError(f"`{program}` returned a non-zero exit code")


def check_real_line(now: Line, correct: Line) -> bool:
    """
    Checks if two lines of real number values are approximately equal within a tolerance.

    Parameters
    ----------
    now : Line
        The current line of values to check.
    correct : Line
        The line of correct - or expected - values to compare against.

    Returns
    -------
    bool
        True if the lines are approximately equal within the specified tolerance, False otherwise.

    Notes
    -----
    The method considers relative and vector relative tolerances depending on the type of the line.

    Examples
    --------
    >>> now = Line(values=[1.0, 2.0, 3.0], type="REL", tolerance=0.01, lineno=0)
    >>> correct = Line(values=[1.01, 2.01, 3.0], type="REL", tolerance=0.01, lineno=0)
    >>> check_real_line(now, correct)
    True
    """
    assert now.tolerance is not None

    if len(now.values) != len(correct.values):
        raise RuntimeError(
            f'Number of components in "output.now" line: {now}\n'
            f'Does not match number in "output.correct:  {correct}'
        )

    vec_amp = 0
    bad_at = -1
    bad_diff_val = 0
    bad_abs_val = 0

    if now.type == "VEC_REL":
        vec_amp = 0
        for ix, (now1, correct1) in enumerate(list(zip(now.values, correct.values))):
            assert isinstance(now1, float)
            assert isinstance(correct1, float)
            vec_amp += ((abs(now1) + abs(correct1)) / 2) ** 2
        vec_amp = math.sqrt(vec_amp)

    for ix, (now1, correct1) in enumerate(zip(now.values, correct.values), start=1):
        assert isinstance(now1, float)
        assert isinstance(correct1, float)

        diff_val = abs(now1 - correct1)
        abs_val = (abs(now1) + abs(correct1)) / 2
        if now.type == "REL":
            factor = abs_val
        elif now.type == "VEC_REL":
            factor = vec_amp
        else:
            factor = 1.0

        # TODO: replace with math.isclose?
        if diff_val > factor * now.tolerance and diff_val > bad_diff_val:
            bad_at = ix
            bad_diff_val = diff_val
            bad_abs_val = abs_val

            if len(now.values) != 1:
                print(f"Regression test failed for datum number: {ix}")

    if bad_at > -1:
        print("")
        print(f"Regression test local failure: {now.name}")

        print(
            f'   Data from "output.now":     {now.values}\n'
            f'   Data from "output.correct": {correct.values}\n'
            f"   Diff: {bad_diff_val} Diff/Val: {abs(bad_diff_val) / bad_abs_val}"
        )
        return False
    return True


def compare_lines(now: Line, correct: Line) -> bool:
    """
    Compares two Line objects to determine if they match based on their type and values.

    Parameters
    ----------
    now : Line
        The current Line object to be compared.
    correct : Line
        The correct - or expected - Line object to be compared against.

    Returns
    -------
    bool
        True if the lines match based on their type and values.
        False if this should be marked as an error to be reported after all
        other tests are done.  These are referred to as 'local' errors.

    Raises
    ------
    RuntimeError
        If the expected types are different.
    ValueError
        If the type of the `now` line is not recognized as one of the valid types.
    """
    if now.type != correct.type:
        raise RuntimeError(
            f'Identification string for a line in "output.now":    {now.type}'
            f'Does not match corresponding ID in "output.correct": {correct.type}'
        )

    if now.type == "STR":
        assert now.values == correct.values
        return True
    if now.type in {"ABS", "REL", "VEC_REL"}:
        return check_real_line(now, correct)
    raise ValueError(
        f'Bad data ID string in "output.now" file: {now.type!r}\n'
        f"Should be one of: STR, REL, VEC_REL, or ABS."
    )


@fortran_tests
def test_fortran(bmad_bin: pathlib.Path, test_name: str) -> None:
    """
    Runs a Fortran test and compares its output to the expected output.

    Parameters
    ----------
    test_name : str
        The name of the test to run, which corresponds to a directory
        in the TESTS_ROOT containing the necessary test files.
    """
    test_path = TESTS_ROOT / test_name

    if not test_path.exists():
        raise FileNotFoundError(f"Test path does not exist: {test_path}")

    os.chdir(test_path)

    run_fortran_code(bmad_bin, test_path)

    # Compare the output of the program "output.now" to the expected output "output.correct"
    now_lines = load_output_file(test_path / "output.now")
    correct_lines = load_output_file(test_path / "output.correct")
    num_local_failures = 0

    for now_line, correct_line in zip(now_lines, correct_lines):
        if not compare_lines(now_line, correct_line):
            num_local_failures += 1

    if len(now_lines) > len(correct_lines):
        raise ValueError("More lines in output.now than in output.correct")
    if len(now_lines) > len(correct_lines):
        raise ValueError("More lines in output.correct than in output.now")
    if num_local_failures:
        raise ValueError("One or more local failures were found")

    print(f"{test_name} 'output.now' and 'output.correct' match.")
