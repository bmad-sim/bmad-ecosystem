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
- --verbose:    print a one-line substitution summary to stderr

Parsing details:
  - Fortran `!` comments (outside quoted strings) are stripped before
    scanning for `&symbolic_number` namelists, so commented-out
    definitions are correctly ignored.
  - The namelist terminator `/` is distinguished from the division
    operator by context: a `/` that is immediately adjacent to a digit,
    identifier character, or `)` on both sides is treated as division.
    The terminator is the first `/` that does NOT look like division.

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


def _strip_comments(text):
    """Strip Fortran !-comments (outside quoted strings) from text.

    Returns a new string with the same line structure but with comment
    portions replaced by spaces (to preserve column positions).
    """
    out = []
    in_quote = False
    quote_ch = ''
    for ch in text:
        if not in_quote and ch in ("'", '"'):
            in_quote = True
            quote_ch = ch
            out.append(ch)
        elif in_quote and ch == quote_ch:
            in_quote = False
            out.append(ch)
        elif not in_quote and ch == '!':
            # Replace rest of line with spaces up to the newline.
            out.append(' ')  # replace the '!' itself
            # We'll handle the rest via a different approach below.
            # Actually, let's do this line-by-line instead.
            pass
        else:
            out.append(ch)
    # The char-by-char approach above doesn't handle "rest of line" well.
    # Redo with a line-by-line approach.
    lines = text.split('\n')
    result = []
    for line in lines:
        new_line = []
        in_q = False
        q_ch = ''
        for c in line:
            if not in_q and c in ("'", '"'):
                in_q = True
                q_ch = c
                new_line.append(c)
            elif in_q and c == q_ch:
                in_q = False
                new_line.append(c)
            elif not in_q and c == '!':
                # Rest of line is comment; stop
                break
            else:
                new_line.append(c)
        result.append(''.join(new_line))
    return '\n'.join(result)


def _is_expr_char(ch):
    """Return True if ch could be part of an expression operand (digit,
    letter, underscore, closing paren, or dot)."""
    return ch.isalnum() or ch in ('_', ')', '.')


def _find_namelist_terminator(text, start):
    """Find the terminator `/` for a namelist starting at `start`.

    The terminator is the first `/` (outside quotes) that is NOT a division
    operator. A `/` is considered division if it has expression characters
    on both sides (ignoring whitespace).

    Returns the index of the terminator `/`, or -1 if not found.
    """
    i = start
    n = len(text)
    in_quote = False
    quote_ch = ''
    while i < n:
        ch = text[i]
        if not in_quote and ch in ("'", '"'):
            in_quote = True
            quote_ch = ch
            i += 1
            continue
        if in_quote:
            if ch == quote_ch:
                in_quote = False
            i += 1
            continue
        if ch == '/':
            # Check if this looks like division: expr_char on left AND right
            # (skipping whitespace).
            left_ok = False
            j = i - 1
            while j >= start and text[j] in (' ', '\t'):
                j -= 1
            if j >= start and _is_expr_char(text[j]):
                left_ok = True

            right_ok = False
            k = i + 1
            while k < n and text[k] in (' ', '\t'):
                k += 1
            if k < n and (text[k].isalnum() or text[k] in ('_', '(', '+', '-', '.')):
                right_ok = True

            if left_ok and right_ok:
                # This is division, skip it.
                i += 1
                continue
            else:
                return i
        i += 1
    return -1


def collect_symbols(text, known_symbols):
    """Parse &symbolic_number blocks and return newly resolved symbols.

    Strips Fortran comments before scanning so that commented-out
    definitions are ignored. Updates known_symbols in place and returns it.
    """
    stripped = _strip_comments(text)
    # Find each &symbolic_number header
    header_re = re.compile(r'&symbolic_number\s+(\w+)\s*=\s*', re.IGNORECASE)
    pos = 0
    while True:
        m = header_re.search(stripped, pos)
        if not m:
            break
        value_start = m.end()
        # Find the terminator `/` starting from after `&symbolic_number`
        term_idx = _find_namelist_terminator(stripped, m.start())
        if term_idx < 0:
            # No terminator found; skip
            pos = value_start
            continue
        value_expr = stripped[value_start:term_idx].strip()
        name = m.group(1).upper()
        value = _eval_value(value_expr, known_symbols)
        if value is not None:
            known_symbols[name] = value
        pos = term_idx + 1
    return known_symbols


def _is_word_char(ch):
    return ch.isalnum() or ch == '_'


def substitute_symbols(text, symbols):
    """Substitute symbolic names with their numeric values in `text`.

    Uses a single compiled regex alternation for O(n) scanning instead of
    O(n * num_symbols). Respects quote boundaries and avoids substituting
    in contexts where a numeric value is not expected. See module docstring
    for rules.
    """
    if not symbols:
        return text, 0

    # Build a single regex alternation of all symbol names (case-insensitive,
    # word-boundary delimited).
    names_upper = sorted(symbols.keys(), key=len, reverse=True)
    values_str = {name: _format_value(symbols[name]) for name in names_upper}
    # Escape names for regex safety (they should be \w+ but be safe).
    pattern = r'\b(' + '|'.join(re.escape(n) for n in names_upper) + r')\b'
    symbol_re = re.compile(pattern, re.IGNORECASE)

    # Pre-compute quote regions to skip.
    quote_ranges = []
    in_quote = False
    quote_ch = ''
    qstart = 0
    for i, ch in enumerate(text):
        if not in_quote and ch in ("'", '"'):
            in_quote = True
            quote_ch = ch
            qstart = i
        elif in_quote and ch == quote_ch:
            in_quote = False
            quote_ranges.append((qstart, i))

    def _in_quotes(pos, end):
        for qs, qe in quote_ranges:
            if pos >= qs and end <= qe + 1:
                return True
            if qs > end:
                break
        return False

    sub_count = 0

    def _replacer(m):
        nonlocal sub_count
        start = m.start()
        end = m.end()
        name = m.group(1).upper()

        # Skip if inside quotes.
        if _in_quotes(start, end):
            return m.group(0)
        # Skip if preceded by '%'.
        if start > 0 and text[start - 1] == '%':
            return m.group(0)
        # Skip if followed by '=' or '(' (after optional spaces).
        j = end
        n = len(text)
        while j < n and text[j] == ' ':
            j += 1
        if j < n and text[j] in ('=', '('):
            return m.group(0)

        sub_count += 1
        return values_str[name]

    result = symbol_re.sub(_replacer, text)
    return result, sub_count


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
    ap.add_argument('--verbose', action='store_true', default=False,
                    help='print substitution summary to stderr')
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

    symbols_before = len(known_symbols)
    collect_symbols(text, known_symbols)
    new_text, sub_count = substitute_symbols(text, known_symbols)

    if args.verbose:
        new_syms = len(known_symbols) - symbols_before
        print(f'tao_preprocess: {len(known_symbols)} symbols '
              f'({new_syms} new), {sub_count} substitutions',
              file=sys.stderr)

    with open(args.output, 'w') as f:
        f.write(new_text)

    if args.symbols_file:
        with open(args.symbols_file, 'w') as f:
            json.dump(known_symbols, f)


if __name__ == '__main__':
    main()
