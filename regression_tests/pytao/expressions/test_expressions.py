from __future__ import annotations

import pathlib
from typing import cast

import numpy as np
import pytest
from pytao import SubprocessTao, Tao

EXPR_ROOT = pathlib.Path(__file__).resolve().parent


def rms(vec):
    vec = np.asarray(vec)
    return np.sqrt(np.sum(vec**2) / len(vec) - np.mean(vec) ** 2)


basic_arithmetic_tests = [
    pytest.param("3 + 4", 7, id="3 + 4"),
    pytest.param("10 - 2", 8, id="10 - 2"),
    pytest.param("6 * 5", 30, id="6 * 5"),
    pytest.param("20 / 4", 5, id="20 / 4"),
    pytest.param("2^3", 8, id="2^3"),
    pytest.param("-5", -5, id="-5"),
    pytest.param("--5", 5, id="--5"),
]

rounding_tests = [
    pytest.param("floor(2.9)", 2, id="floor_pos"),
    pytest.param("floor(-2.1)", -3, id="floor_neg"),
    pytest.param("ceiling(2.1)", 3, id="ceiling_pos"),
    pytest.param("ceiling(-2.9)", -2, id="ceiling_neg"),
    pytest.param("nint(2.6)", 3, id="nint_round_up"),
    pytest.param("nint(2.4)", 2, id="nint_round_down"),
    pytest.param("int(2.9)", 2, id="int_trunc_pos"),
    pytest.param("int(-2.9)", -2, id="int_trunc_neg"),
    pytest.param("sign(-5.5)", -1, id="sign_neg"),
    pytest.param("sign(5.5)", 1, id="sign_pos"),
    pytest.param("sign(0)", 0, id="sign_zero"),
]

vector_tests = [
    pytest.param("[1, 2, 3] + [4, 5, 6]", [5, 7, 9], id="[1,2,3]+[4,5,6]"),
    pytest.param("[1, 2, 3] * [4, 5, 6]", [4, 10, 18], id="[1,2,3]*[4,5,6]"),
    pytest.param("[4, 5, 6] / [1, 2, 3]", [4, 2.5, 2], id="[4,5,6]/[1,2,3]"),
    pytest.param(
        "[2, 2, 2] ^ [3, 2, 1]", np.power([2, 2, 2], [3, 2, 1]), id="pow_vectors"
    ),
    pytest.param("[1, 2, 3] * 2", [2, 4, 6], id="[1,2,3]*2"),
    pytest.param("2 + [1, 2, 3]", [3, 4, 5], id="2+[1,2,3]"),
    pytest.param("-[1, 2, 3]", [-1, -2, -3], id="unary_minus_vec"),
    pytest.param("[1, 2, 3] + [10]", [11, 12, 13], id="vec_plus_vec_len_1"),
    pytest.param("[10] + [1, 2, 3]", [11, 12, 13], id="vec_len_1_plus_vec"),
    pytest.param("-[1, -2, 3]", [-1, 2, -3], id="unary_minus_vec"),
    pytest.param("abs(-[1, -2, 3])", [1, 2, 3], id="abs_unary_minus_vec"),
    pytest.param("sqrt([4, 16] * 4)", [4, 8], id="sqrt_vec_broadcast"),
]

standard_math_tests = [
    pytest.param("sin(0.5)", np.sin(0.5), id="sin(0.5)"),
    pytest.param("cos(0.5)", np.cos(0.5), id="cos(0.5)"),
    pytest.param("tan(0.5)", np.tan(0.5), id="tan(0.5)"),
    pytest.param("asin(0.5)", np.asin(0.5), id="asin(0.5)"),
    pytest.param("acos(0.5)", np.acos(0.5), id="acos(0.5)"),
    pytest.param("atan(0.5)", np.atan(0.5), id="atan(0.5)"),
    pytest.param("sinh(0.5)", np.sinh(0.5), id="sinh(0.5)"),
    pytest.param("cosh(0.5)", np.cosh(0.5), id="cosh(0.5)"),
    pytest.param("tanh(0.5)", np.tanh(0.5), id="tanh(0.5)"),
    pytest.param("sqrt(4)", 2, id="sqrt(4)"),
    pytest.param("sqrt(-1)", np.nan, id="sqrt(-1)"),
    pytest.param("exp(1)", np.exp(1), id="exp(1)"),
    pytest.param(f"log({np.e})", 1, id="log(e)"),
    pytest.param("abs(-10.5)", 10.5, id="abs(-10.5)"),
    pytest.param("sqrt([4, 9, 16])", np.sqrt([4, 9, 16]), id="sqrt_vec"),
    pytest.param("abs([-1, 2, -3])", np.abs([-1, 2, -3]), id="abs_vec"),
    pytest.param("cot(0.5)", 1 / np.tan(0.5), id="cot"),
    pytest.param("sec(0.5)", 1 / np.cos(0.5), id="sec"),
    pytest.param("csc(0.5)", 1 / np.sin(0.5), id="csc"),
    pytest.param("coth(0.5)", 1 / np.tanh(0.5), id="coth"),
    pytest.param("acosh(2.0)", np.arccosh(2.0), id="acosh"),
    pytest.param("asinh(2.0)", np.arcsinh(2.0), id="asinh"),
    pytest.param("atanh(0.5)", np.arctanh(0.5), id="atanh"),
]

two_arg_function_tests = [
    pytest.param("modulo(5.5, 2)", 1.5, id="modulo(5.5, 2)"),
    # Not a thing:
    # pytest.param("modulo([5.5, 6.5], 2)", [1.5, 0.5], id="modulo_vec"),
    pytest.param("atan2(1, 1)", np.atan2(1, 1), id="atan2(1, 1)"),
    # Not a thing:
    # pytest.param(
    #     "atan2([1, 0, -1], [1, 1, 1])",
    #     np.arctan2([1, 0, -1], [1, 1, 1]),
    #     id="atan2_vec",
    # ),
]

reduction_function_tests = [
    pytest.param("min([1, 2, 3, 4, 5])", 1, id="min"),
    pytest.param("max([1, 2, 3, 4, 5])", 5, id="max"),
    pytest.param("sum([1, 2, 3, 4, 5])", 15, id="sum"),
    pytest.param("mean([1, 2, 3, 4, 5])", 3, id="mean"),
    pytest.param("average([1, 2, 3, 4, 5])", 3, id="average"),
    pytest.param("rms([1, 2, 3, 4, 5])", rms([1, 2, 3, 4, 5]), id="rms"),
]

special_function_tests = [
    pytest.param("factorial(5)", 120, id="factorial(5)"),
    pytest.param("factorial([3, 4])", [6, 24], id="factorial_vec"),
    pytest.param("sinc(0.5)", np.sin(0.5) / 0.5, id=f"sinc({0.5})"),
]

physics_tests = [
    pytest.param("mass_of(C)", 1.1187803082139433e10, id="mass_electron_ev"),
    pytest.param("anomalous_moment_of(C)", 0.0, id="anomalous_moment_of"),
    pytest.param("charge_of(C-)", -1.0, id="charge_of_neg"),
    pytest.param("charge_of(C--)", -2.0, id="charge_of_neg2"),
    pytest.param("charge_of(C+)", 1.0, id="charge_of_pos"),
    pytest.param("charge_of(C++)", 2.0, id="charge_of_pos2"),
]

misc_tests = [
    pytest.param(
        "360 * modulo(2.2 - 0.1, 1)",
        360 * ((2.2 - 0.1) % 1.0),
        id="nesting_scalar",
    ),
    # Modulo doesn't accept vectors
    # pytest.param(
    #     "10 * modulo([2.2, 3.5] - [0.1, 0.4], 1)",
    #     10 * (np.array([2.1, 3.1]) % 1),
    #     id="nesting_vec",
    # ),
]

EXPRESSION_TEST_CASES = [
    *basic_arithmetic_tests,
    *rounding_tests,
    *vector_tests,
    *standard_math_tests,
    *two_arg_function_tests,
    *reduction_function_tests,
    *special_function_tests,
    *physics_tests,
    *misc_tests,
]


@pytest.fixture(scope="module")
def tao():
    with SubprocessTao(lattice_file=str(EXPR_ROOT / "lat1.bmad"), noplot=True) as tao:
        yield tao


@pytest.mark.parametrize("expression, expected", EXPRESSION_TEST_CASES)
def test_expressions(
    tao: Tao, expression: str, expected: float | list | np.ndarray
) -> None:
    result = cast(np.ndarray, tao.evaluate(expression))
    np.testing.assert_allclose(result, expected, rtol=1e-7, atol=1e-8)


def test_ran(tao: Tao) -> None:
    val = np.asarray(tao.evaluate("ran()"))[0]
    assert 0 <= val <= 1


def test_species_of(tao: Tao) -> None:
    tao.evaluate("species_of(C)")  # ?


@pytest.mark.parametrize(
    "param",
    [
        "",
        "0",
        "1",
    ],
)
def test_ran_gauss_parameter_smoke(tao: Tao, param: str) -> None:
    val = np.asarray(tao.evaluate(f"ran_gauss({param})"))[0]
    print("value=", val)
    # Just checking it doesn't crash


if __name__ == "__main__":
    exit(pytest.main(["-v", __file__]))
