#!/usr/bin/env python3
"""
Tao init-file preprocessor -- pure text substitution.

Substitutes symbolic names with numeric values that were pre-evaluated by
Fortran's tao_evaluate_expression. This script does NOT evaluate expressions
itself; it only consumes symbols exported by Fortran.

The --symbols-file is the ONLY source of symbol definitions. It must be a
JSON file containing a dict of ``{"NAME": number, ...}`` with uppercase
keys. If the file is absent or empty, no substitution is performed and the
output is identical to the input.

Usage:
    python3 tao_preprocess.py INPUT_FILE OUTPUT_FILE [--symbols-file SYMBOLS] [--verbose]

- INPUT_FILE:   path to the original tao init file
- OUTPUT_FILE:  path to write the preprocessed file (same content if no
                substitution is needed)
- --symbols-file SYMBOLS: path to a JSON file of symbols pre-evaluated by
                Fortran. Read only; never written back.
- --verbose:    print a one-line substitution summary to stderr

Substitution rules:
  - Inside single- or double-quoted strings: NOT substituted
  - Immediately preceded by '%' (e.g. ``datum%meas``): NOT substituted
  - Followed by '=' (e.g. ``foo = ...``): NOT substituted
  - Followed by '(' (e.g. ``foo(1)``): NOT substituted
  - Whole-word matches only (not substring)
  - Matching is case-insensitive
"""

import argparse
import json
import os
import re
import sys


def _is_word_char(ch):
    return ch.isalnum() or ch == '_'


def _format_value(v):
    """Format a float value as a Fortran-parseable numeric literal."""
    return f'{v:.15e}'


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


def main():
    ap = argparse.ArgumentParser(description='Preprocess a tao init file.')
    ap.add_argument('input', help='input file path')
    ap.add_argument('output', help='output file path')
    ap.add_argument('--symbols-file', default=None,
                    help='JSON file with symbols pre-evaluated by Fortran '
                         '(read only; never written back)')
    ap.add_argument('--verbose', action='store_true', default=False,
                    help='print substitution summary to stderr')
    args = ap.parse_args()

    # Load symbols from file (if provided and parseable).
    symbols = {}
    if args.symbols_file and os.path.exists(args.symbols_file):
        try:
            with open(args.symbols_file) as f:
                symbols = json.load(f)
        except Exception:
            symbols = {}

    # Read input.
    try:
        with open(args.input, 'r') as f:
            text = f.read()
    except OSError as exc:
        print(f'tao_preprocess: cannot read {args.input}: {exc}',
              file=sys.stderr)
        sys.exit(1)

    # Substitute and write output.
    new_text, sub_count = substitute_symbols(text, symbols)

    if args.verbose:
        print(f'tao_preprocess: {len(symbols)} symbols, '
              f'{sub_count} substitutions',
              file=sys.stderr)

    with open(args.output, 'w') as f:
        f.write(new_text)


if __name__ == '__main__':
    main()
