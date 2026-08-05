import shutil
from pathlib import Path

from pytao import SubprocessTao

BOOKKEEPING_ROOT = Path(__file__).resolve().parent
INPUT_FILES = ("stale_bookkeeping.bmad", "stale_bookkeeping.init")


def test_no_stale_bookkeeping_flags_during_optimization(tmp_path, monkeypatch):
    """
    Regression test for stale bookkeeping flags reported by lattice_bookkeeper
    when the optimizer changes the length of a closed lattice.
    """
    for name in INPUT_FILES:
        shutil.copy(BOOKKEEPING_ROOT / name, tmp_path / name)

    # The init file refers to the lattice by a relative path, and the optimizer
    # writes var1.out into the working directory when it finishes.
    monkeypatch.chdir(tmp_path)

    with SubprocessTao(init_file="stale_bookkeeping.init", noplot=True) as tao:
        tao.cmd("run", raises=False)
        stale = [msg for msg in tao.last_messages if "Stale bookkeeping" in msg.message]

    assert not stale, "\n".join(str(msg) for msg in stale)
