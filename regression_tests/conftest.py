import pathlib

import pytest

TESTS_ROOT = pathlib.Path(__file__).resolve().parent
BMAD_REPO_ROOT = TESTS_ROOT.parent
DEBUG_BIN_PATH = pathlib.Path("..") / "debug" / "bin"
PRODUCTION_BIN_PATH = pathlib.Path("..") / "production" / "bin"

BMAD_DOC_ROOT = BMAD_REPO_ROOT / "bmad-doc"
TAO_EXAMPLES_ROOT = BMAD_DOC_ROOT / "tao_examples"

DEFAULT_BIN_DIR = (
    PRODUCTION_BIN_PATH if PRODUCTION_BIN_PATH.is_dir() else DEBUG_BIN_PATH
)


def pytest_addoption(parser: pytest.Parser):
    parser.addoption(
        "--bmad-bin",
        action="store",
        default=str(DEFAULT_BIN_DIR),
        help="Bmad binary directory, which should include 'tao' and the other regression tests",
    )


@pytest.fixture(scope="module")
def bmad_bin(request: pytest.FixtureRequest) -> pathlib.Path:
    bin = pathlib.Path(str(request.config.getoption("--bmad-bin"))).resolve()
    if not bin.is_dir():
        raise ValueError(f"--bmad-bin path {bin} is not a directory")
    return bin
