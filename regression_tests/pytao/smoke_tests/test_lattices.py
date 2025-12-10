from pathlib import Path
from pytao import SubprocessTao
import pytest


# Get all lattice files from the lattices subdirectory
lattices_dir = Path(__file__).parent / "lattices"
lattice_files = sorted(lattices_dir.glob("*.bmad"))


@pytest.mark.parametrize("lattice_file", lattice_files, ids=lambda p: p.name)
def test_load_lattice(lattice_file):
    """Test that each lattice file can be loaded without errors."""
    with SubprocessTao(lattice_file=lattice_file, noplot=True):
        pass