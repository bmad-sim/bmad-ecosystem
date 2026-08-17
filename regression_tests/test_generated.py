"""
Check that checked-in generated source files are up to date with respect to their inputs.

Generated tables that are maintained by hand drift. The point of generating
tao/code/tao_attrib_resolve_mod.f90 from the structure definitions is that adding a component
to, say, bmad_common_struct cannot silently leave the resolver behind. This test is what
enforces that: it regenerates in memory and fails if the result differs from the checked-in
file. The fix when it fails is to run the generator and commit the result.

Unlike the Fortran tests, this needs no build, only Python.
"""

from __future__ import annotations

import subprocess
import sys

from conftest import BMAD_REPO_ROOT

GENERATORS = [
    "tao/scripts/generate_attrib_tables.py",
]


def test_generated_files_are_up_to_date():
    for script in GENERATORS:
        result = subprocess.run(
            [sys.executable, script, "--check"],
            cwd=BMAD_REPO_ROOT,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"{script} --check failed:\n{result.stdout}\n{result.stderr}"
        )
