from pathlib import Path
from pytao import SubprocessTao
import pytest


# Get all lattice files from the lattices subdirectory
lattices_dir = Path(__file__).parent / "lattices"
lattice_files = sorted(lattices_dir.glob("*.bmad")) + sorted(lattices_dir.glob("*.lat"))


@pytest.mark.parametrize("lat_path", lattice_files, ids=lambda p: p.name)
def test_load_lattice(lat_path):
    """Test that each lattice file can be loaded without errors."""
    with SubprocessTao(lattice_file=lat_path, noplot=True):
        pass