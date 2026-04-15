#!/usr/bin/env python3
"""
Tao init-file preprocessor.

Expands `&symbolic_number` namelist definitions by substituting symbol
references with their numeric values throughout the file. This enables
symbolic numbers to be used in positions where the Fortran namelist
reader cannot handle symbolic names (e.g., the real-valued `meas` field
in datum definitions).

Invoked automatically by tao; end users do not run this script directly.

Usage:
    python3 tao_preprocess.py INPUT_FILE OUTPUT_FILE [--symbols-file SYMBOLS]

- INPUT_FILE:   path to the original tao init file
- OUTPUT_FILE:  path to write the preprocessed file (same content if no
                substitution is needed)
- --symbols-file SYMBOLS: optional path to a JSON file holding symbols
                from previous invocations. Updated on each run so
                subsequent invocations inherit earlier definitions.

Substitution rules (applied only when a symbolic value was successfully
evaluated by this preprocessor):
  - Inside single- or double-quoted strings: NOT substituted
  - Immediately preceded by '%' (e.g. `datum%meas`): NOT substituted
  - Followed by '=' (e.g. `foo = ...`): NOT substituted
  - Followed by '(' (e.g. `foo(1)`): NOT substituted
  - Whole-word matches only (not substring)
  - Matching is case-insensitive
  - `&symbolic_number` namelist blocks are left intact in the output so
    Fortran's expression evaluator continues to see the definitions.
"""

import argparse
import json
import math
import os
import re
import sys


# Safe namespace for evaluating &symbolic_number values.
_MATH_NAMES = {
    'pi': math.pi, 'e': math.e,
    'sqrt': math.sqrt, 'exp': math.exp, 'log': math.log, 'log10': math.log10,
    'sin': math.sin, 'cos': math.cos, 'tan': math.tan,
    'asin': math.asin, 'acos': math.acos, 'atan': math.atan, 'atan2': math.atan2,
    'sinh': math.sinh, 'cosh': math.cosh, 'tanh': math.tanh,
    'abs': abs, 'min': min, 'max': max,
}


def _eval_value(expr_text, known_symbols):
    """Evaluate the RHS of a &symbolic_number. Returns float or None on failure.

    Accepts literal numbers, basic arithmetic, math functions (pi, sqrt, etc.),
    and references to previously-defined symbols. Complex tao-specific
    expressions (e.g. ones referring to lattice data) are not handled here and
    will be left for Fortran to evaluate via its expression evaluator.
    """
    namespace = dict(_MATH_NAMES)
    for name, value in known_symbols.items():
        namespace[name.upper()] = value
        namespace[name.lower()] = value
    # Convert Fortran's '**' and 'd'/'D' exponent to Python form.
    expr_py = expr_text.strip()
    expr_py = re.sub(r'(\d)[dD]([+\-]?\d)', r'\1e\2', expr_py)
    try:
        result = eval(expr_py, {'__builtins__': {}}, namespace)
        return float(result)
    except Exception:
        return None


# Matches a full &symbolic_number NAME = VALUE / block (may span lines).
_SYMBOLIC_NUMBER_RE = re.compile(
    r'&symbolic_number\s+(\w+)\s*=\s*([^/]*?)\s*/',
    re.IGNORECASE | re.DOTALL,
)


def collect_symbols(text, known_symbols):
    """Parse &symbolic_number blocks and return newly resolved symbols.

    Updates known_symbols in place and returns it.
    """
    for match in _SYMBOLIC_NUMBER_RE.finditer(text):
        name = match.group(1).upper()
        value_expr = match.group(2).strip()
        value = _eval_value(value_expr, known_symbols)
        if value is not None:
            known_symbols[name] = value
    return known_symbols


def _is_word_char(ch):
    return ch.isalnum() or ch == '_'


def substitute_symbols(text, symbols):
    """Substitute symbolic names with their numeric values in `text`.

    Respects quote boundaries and avoids substituting in contexts where
    a numeric value is not expected. See module docstring for rules.
    """
    if not symbols:
        return text

    # Sort names longest-first so a longer name matches before a shorter
    # prefix of it (e.g. 'TARGET_BETA_A' before 'TARGET_BETA').
    names_upper = sorted(symbols.keys(), key=len, reverse=True)
    values_str = {name: _format_value(symbols[name]) for name in names_upper}

    out = []
    i = 0
    n = len(text)
    in_quote = False
    quote_ch = ''

    while i < n:
        ch = text[i]

        # Track quote state.
        if not in_quote and ch in ("'", '"'):
            in_quote = True
            quote_ch = ch
            out.append(ch)
            i += 1
            continue
        if in_quote and ch == quote_ch:
            in_quote = False
            out.append(ch)
            i += 1
            continue
        if in_quote:
            out.append(ch)
            i += 1
            continue

        # Try to match a known symbol at this position.
        matched = False
        for name in names_upper:
            name_len = len(name)
            if i + name_len > n:
                continue
            if text[i:i + name_len].upper() != name:
                continue
            # Word boundaries.
            if i > 0 and _is_word_char(text[i - 1]):
                continue
            end = i + name_len
            if end < n and _is_word_char(text[end]):
                continue
            # Context: skip component access (preceded by '%').
            if i > 0 and text[i - 1] == '%':
                continue
            # Context: skip namelist-var assignments and calls/arrays.
            j = end
            while j < n and text[j] == ' ':
                j += 1
            if j < n and text[j] in ('=', '('):
                continue
            out.append(values_str[name])
            i = end
            matched = True
            break

        if not matched:
            out.append(ch)
            i += 1

    return ''.join(out)


def _format_value(v):
    """Format a float value as a Fortran-parseable numeric literal."""
    return f'{v:.15e}'


def main():
    ap = argparse.ArgumentParser(description='Preprocess a tao init file.')
    ap.add_argument('input', help='input file path')
    ap.add_argument('output', help='output file path')
    ap.add_argument('--symbols-file', default=None,
                    help='JSON file with symbols from prior invocations '
                         '(read and updated)')
    args = ap.parse_args()

    known_symbols = {}
    if args.symbols_file and os.path.exists(args.symbols_file):
        try:
            with open(args.symbols_file) as f:
                known_symbols = json.load(f)
        except Exception:
            known_symbols = {}

    try:
        with open(args.input, 'r') as f:
            text = f.read()
    except OSError as exc:
        print(f'tao_preprocess: cannot read {args.input}: {exc}',
              file=sys.stderr)
        sys.exit(1)

    collect_symbols(text, known_symbols)
    new_text = substitute_symbols(text, known_symbols)

    with open(args.output, 'w') as f:
        f.write(new_text)

    if args.symbols_file:
        with open(args.symbols_file, 'w') as f:
            json.dump(known_symbols, f)


if __name__ == '__main__':
    main()
