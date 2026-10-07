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
    fields: dict[str, str] = {}
    matches: list[str] = []
    for item in out:
        name, _type, _settable, value = item.split(";", 3)
        if name.startswith("match["):
            matches.append(value)
        else:
            fields[name] = value
    return fields["word"], fields["context"], matches


def test_output_parses_with_generic_parameter_list_parser(tao):
    """
    The output must parse with PyTao's generic parser so that bindings regenerated
    against a PyTao without a dedicated ``parse_complete`` still work.
    """
    from pytao.util.parsers import parse_tao_python_data

    data = parse_tao_python_data(tao.cmd('pipe complete "sho"'), clean_key=False)
    assert data["word"] == "sho"
    assert data["context"] == "LIST"
    assert data["match[1]"] == "show"


def test_empty_line_lists_all_commands(tao):
    word, context, matches = complete(tao, "")
    assert word == ""
    assert context == "LIST"
    assert "show" in matches
    assert "set" in matches
    assert len(matches) >= 46
    # Undocumented internal commands are accepted by the parser but not offered.
    assert "debug" not in matches


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


def test_set_global_lists_struct_components(tao):
    _, context, matches = complete(tao, "set global ")
    assert context == "LIST"
    assert "n_opti_cycles" in matches
    assert "quiet" in matches
    assert len(matches) > 30


@pytest.mark.parametrize(
    ("line", "expected_match"),
    [
        ("set global track_t", "track_type"),
        ("set global phase_u", "phase_units"),
        ("set bmad_com max_aperture_l", "max_aperture_limit"),
        ("set beam_init n_par", "n_particle"),
        ("set space_charge_com ds_track_s", "ds_track_step"),
        ("set ptc_com exact_mo", "exact_model"),
        ("set beam track_s", "track_start"),
        ("set plot_page tit", "title"),
        ("set element Q01W k", "k1"),
    ],
)
def test_set_component_names(tao, line, expected_match):
    _, context, matches = complete(tao, line)
    assert context == "LIST"
    assert expected_match in matches


@pytest.mark.parametrize(
    ("line", "expected_match"),
    [
        ("show -app", "-append"),
        ("show lattice -orb", "-orbit"),
        ("show lat -6d", "-6d_radiation_integrals"),
        ("show element Q01W -floor", "-floor_coords"),
        ("set -up", "-update"),
        ("change -si", "-silent"),
        ("place -no_b", "-no_buffer"),
        ("pipe -nop", "-noprint"),
    ],
)
def test_switch_completion(tao, line, expected_match):
    _, context, matches = complete(tao, line)
    assert context == "LIST"
    assert expected_match in matches


def test_unknown_switch_context_offers_nothing(tao):
    _, context, matches = complete(tao, "show wave -")
    assert matches == []


def test_show_data_names(tao):
    _, context, matches = complete(tao, "show data orbit.")
    assert context == "LIST"
    assert "orbit.x" in matches


def test_show_var_names(tao):
    _, context, matches = complete(tao, "show var quad")
    assert context == "LIST"
    assert "quad_k1" in matches


def test_place_regions_then_templates(tao):
    _, context, regions = complete(tao, "place ")
    assert context == "LIST"
    assert regions
    _, context, templates = complete(tao, f"place {regions[0]} ")
    assert context == "LIST"
    assert templates
    assert regions[0] not in templates


def test_show_plot_names(tao):
    _, context, matches = complete(tao, "show plot ")
    assert context == "LIST"
    assert matches


@pytest.mark.parametrize(
    ("line", "expected_match"),
    [
        ("set element Q01W tracking_method = ", "runge_kutta"),
        ("set element Q01W tracking_method = run", "runge_kutta"),
        ("set element Q01W tracking_method=r", "tracking_method=runge_kutta"),
        ("set element Q01W tracking_method =b", "=bmad_standard"),
        ("set element Q01W field_master = ", "T"),
        ("set global track_type = ", "beam"),
        ("set global rf_on = ", "F"),
        ("set bmad_com radiation_damping_on = ", "T"),
        ("set beam_init random_engine = q", "quasi"),
    ],
)
def test_set_value_completion(tao, line, expected_match):
    _, context, matches = complete(tao, line)
    assert context == "LIST"
    assert expected_match in matches


def test_set_value_glued_form_keeps_attribute_prefix(tao):
    word, context, matches = complete(tao, "set element Q01W tracking_method=r")
    assert word == "tracking_method=r"
    assert context == "LIST"
    assert matches
    assert all(match.startswith("tracking_method=") for match in matches)


def test_set_numeric_attribute_has_no_value_candidates(tao):
    _, context, matches = complete(tao, "set element Q01W k1 = ")
    assert context == "LIST"
    assert matches == []


@pytest.mark.parametrize(
    ("line", "expected_match"),
    [
        ("set element quad::* ", "k1"),
        ("set element quad::* spin_fringe_on = ", "T"),
        ("set element quad::* tracking_method = ", "runge_kutta"),
        ("set element Q0*W fringe_t", "fringe_type"),
        ("set element 1:10 ", "l"),
        ("set element 1@Q01W k", "k1"),
    ],
)
def test_element_selector_syntax(tao, line, expected_match):
    _, context, matches = complete(tao, line)
    assert context == "LIST"
    assert expected_match in matches


def test_mixed_selector_offers_only_common_attributes(tao):
    _, context, matches = complete(tao, "set element * ")
    assert context == "LIST"
    assert "l" in matches
    assert "k1" not in matches


def test_key_prefix_completes_elements_of_that_type(tao):
    word, context, matches = complete(tao, "set element quad::Q0")
    assert word == "quad::Q0"
    assert context == "LIST"
    assert "quad::Q01W" in matches
    assert all(match.startswith("quad::Q0") for match in matches)


def test_element_types_offered_as_selectors(tao):
    _, context, matches = complete(tao, "set element qu")
    assert context == "LIST"
    assert "quadrupole::" in matches


def test_only_element_types_present_in_lattice_are_offered(tao):
    _, _, matches = complete(tao, "set element ")
    types = [match for match in matches if match.endswith("::")]
    assert "quadrupole::" in types
    assert "crystal::" not in types
    assert "beginning_ele::" not in types


def test_universe_prefix_kept_on_element_candidates(tao):
    _, context, matches = complete(tao, "set element 1@Q0")
    assert context == "LIST"
    assert "1@Q01W" in matches


def test_set_element_attributes_require_known_element(tao):
    _, context, matches = complete(tao, "set element NO_SUCH_ELE ")
    assert context == "LIST"
    assert matches == []


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
