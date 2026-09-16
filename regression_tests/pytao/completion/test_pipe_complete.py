from pathlib import Path

import pytest
from pytao import SubprocessTao

CESR_INIT = Path(__file__).resolve().parents[2] / "pipe_test" / "cesr" / "tao.init"


@pytest.fixture(scope="module")
def tao():
    with SubprocessTao(init_file=str(CESR_INIT), noplot=True) as tao:
        yield tao


def complete(tao: SubprocessTao, line: str) -> tuple[str, str, list[str]]:
    """
    Run ``pipe complete`` on a partial command line.

    Returns
    -------
    tuple[str, str, list[str]]
        The partial word being completed, the context (LIST, FILE, or NONE),
        and the candidate list.
    """
    out = tao.cmd(f'pipe complete "{line}"')
    if isinstance(out, str):
        out = [out]
    word, context = out[0].split(";")
    return word, context, out[1:]


def test_empty_line_lists_all_commands(tao):
    word, context, matches = complete(tao, "")
    assert word == ""
    assert context == "LIST"
    assert "show" in matches
    assert "set" in matches
    assert len(matches) >= 49


def test_unique_command_prefix(tao):
    word, context, matches = complete(tao, "sho")
    assert word == "sho"
    assert context == "LIST"
    assert matches == ["show"]


def test_show_subcommands(tao):
    word, context, matches = complete(tao, "show ")
    assert word == ""
    assert context == "LIST"
    assert "element" in matches
    assert "lattice" in matches


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        ("set gl", ["global"]),
        ("show el", ["element"]),
    ],
)
def test_subcommand_prefix(tao, line, expected):
    _, context, matches = complete(tao, line)
    assert context == "LIST"
    assert matches == expected


def test_pipe_subcommands(tao):
    word, context, matches = complete(tao, "pipe lat_")
    assert word == "lat_"
    assert context == "LIST"
    assert matches
    assert all(m.startswith("lat_") for m in matches)
    assert "lat_ele_list" in matches


def test_help_pipe_subcommands(tao):
    _, context, matches = complete(tao, "help pipe comp")
    assert context == "LIST"
    assert matches == ["complete"]


def test_element_names(tao):
    word, context, matches = complete(tao, "show element Q")
    assert word == "Q"
    assert context == "LIST"
    assert matches
    assert all(m.upper().startswith("Q") for m in matches)


def test_data_names(tao):
    _, context, matches = complete(tao, "use data ")
    assert context == "LIST"
    assert matches


def test_var_names(tao):
    _, context, matches = complete(tao, "veto var ")
    assert context == "LIST"
    assert matches


def test_call_completes_file_names(tao):
    word, context, matches = complete(tao, "call ")
    assert word == ""
    assert context == "FILE"
    assert matches == []


def test_unknown_context(tao):
    _, context, matches = complete(tao, "xyzzy plugh ")
    assert context == "NONE"
    assert matches == []


def test_trailing_blank_differs_from_no_blank(tao):
    word_no_blank, _, matches_no_blank = complete(tao, "show")
    word_blank, _, matches_blank = complete(tao, "show ")
    assert word_no_blank == "show"
    assert word_blank == ""
    assert matches_no_blank == ["show"]
    assert "element" in matches_blank


def test_help_pipe_complete_documented(tao):
    out = tao.cmd("help pipe complete")
    if isinstance(out, str):
        out = [out]
    assert any("pipe complete" in line for line in out)
