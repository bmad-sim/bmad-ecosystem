#!/usr/bin/env python3
"""
Parser for Fortran derived type ("struct") definitions as written in the Bmad ecosystem.

This module is deliberately narrow: it understands the declaration styles actually used in
bmad_struct.f90, tao_struct.f90, tao_input_struct.f90, quick_plot_struct.f90, etc.  It is not
a general Fortran parser.

The output is a dict of Struct objects keyed by lower case type name.  It is consumed by
generate_attrib_tables.py which emits the name -> component pointer resolver code.

Declaration styles handled:
    real(rp) :: a = 0, b(3) = 0, mat(6,6) = 0
    integer :: n = 0
    integer(8) :: big = 0
    character(40) :: name = ''
    logical :: on = .false.
    complex(rp) :: vec(8) = 0
    type (twiss_struct) :: x = twiss_struct()
    type (coord_struct), allocatable :: orb(:)
    real(rp), pointer :: p => null()
    real(rp) c0, c1, n_exp                        ! Old style, no "::"
    type (twiss_struct) a, b, c                   ! Old style, no "::"
"""

import re

# Intrinsic types recognized at the start of a declaration.

INTRINSIC_TYPES = ('real', 'integer', 'logical', 'character', 'complex', 'double')

RE_TYPE_START = re.compile(r'^\s*type\s*(?:,[^:]*::)?\s*([a-z]\w*)\s*$', re.I)
RE_TYPE_END = re.compile(r'^\s*end\s*type', re.I)
RE_CONTAINS = re.compile(r'^\s*contains\s*$', re.I)


class Component:
  """One component of a derived type."""

  def __init__(self, name, base_type, kind, dims, allocatable, pointer, init, comment):
    self.name = name                  # Lower case component name.
    self.base_type = base_type        # 'real', 'integer', 'logical', 'character', 'complex', 'type'
    self.kind = kind                  # 'rp', '8', '40', '*', 'twiss_struct', ... or '' if none.
    self.dims = dims                  # List of dimension strings, eg ['6'] or ['0:6'] or [':', ':'].
    self.allocatable = allocatable
    self.pointer = pointer
    self.init = init                  # Initializer text or '' if none.
    self.comment = comment

  @property
  def rank(self):
    return len(self.dims)

  @property
  def is_struct(self):
    return self.base_type == 'type'

  @property
  def is_deferred(self):
    """True if the shape is not known at compile time (allocatable or pointer array)."""
    return any(d == ':' for d in self.dims)

  def bounds(self, idim):
    """Return (lo, hi) as strings for dimension idim, or None if deferred."""
    d = self.dims[idim]
    if d == ':': return None
    if ':' in d:
      lo, hi = d.split(':', 1)
      return (lo.strip(), hi.strip())
    return ('1', d.strip())

  def __repr__(self):
    d = '(' + ','.join(self.dims) + ')' if self.dims else ''
    k = '(' + self.kind + ')' if self.kind else ''
    return f'{self.base_type}{k} :: {self.name}{d}'


class Struct:
  """One derived type definition."""

  def __init__(self, name, source_file):
    self.name = name                  # Lower case type name.
    self.source_file = source_file
    self.components = []

  def __repr__(self):
    return f'<Struct {self.name}: {len(self.components)} components>'


def strip_comment(line):
  """Remove a trailing Fortran comment, respecting quoted strings."""
  out = []
  quote = None
  for ch in line:
    if quote:
      out.append(ch)
      if ch == quote: quote = None
    elif ch in ('"', "'"):
      quote = ch
      out.append(ch)
    elif ch == '!':
      return ''.join(out), line[len(out) + 1:].strip()
    else:
      out.append(ch)
  return ''.join(out), ''


def join_continuations(lines):
  """Join Fortran '&' continuation lines. Returns list of (text, comment) pairs."""
  out = []
  pending = ''
  pending_comment = ''
  for raw in lines:
    code, comment = strip_comment(raw.rstrip('\n'))
    code = code.rstrip()
    if pending:
      code = pending + ' ' + code.lstrip().lstrip('&').lstrip()
      comment = pending_comment or comment
    if code.endswith('&'):
      pending = code[:-1].rstrip()
      pending_comment = comment
      continue
    pending = ''
    pending_comment = ''
    out.append((code, comment))
  if pending: out.append((pending, pending_comment))
  return out


def split_top_level(text, sep=','):
  """Split on `sep` that is not nested inside (), [], or quotes."""
  parts = []
  depth = 0
  quote = None
  cur = []
  for ch in text:
    if quote:
      cur.append(ch)
      if ch == quote: quote = None
      continue
    if ch in ('"', "'"):
      quote = ch
      cur.append(ch)
      continue
    if ch in '([':
      depth += 1
    elif ch in ')]':
      depth -= 1
    if ch == sep and depth == 0:
      parts.append(''.join(cur))
      cur = []
      continue
    cur.append(ch)
  parts.append(''.join(cur))
  return [p.strip() for p in parts]


def match_paren(text, start):
  """Given text[start] == '(', return index just past the matching ')'."""
  depth = 0
  for i in range(start, len(text)):
    if text[i] == '(': depth += 1
    elif text[i] == ')':
      depth -= 1
      if depth == 0: return i + 1
  raise ValueError('Unbalanced parenthesis in: ' + text)


def parse_type_spec(text):
  """
  Parse the leading type specification of a declaration.
  Returns (base_type, kind, rest_of_line) or None if this is not a declaration.
  """
  m = re.match(r'^\s*([a-z]\w*)', text, re.I)
  if not m: return None
  word = m.group(1).lower()
  pos = m.end()

  if word == 'type':
    # "type (foo_struct)" -- the space is optional.
    m2 = re.match(r'\s*\(', text[pos:])
    if not m2: return None
    lp = pos + m2.end() - 1
    rp = match_paren(text, lp)
    kind = text[lp + 1:rp - 1].strip().lower()
    return ('type', kind, text[rp:])

  if word not in INTRINSIC_TYPES: return None

  if word == 'double':
    # "double precision"
    m2 = re.match(r'\s+precision\b', text[pos:], re.I)
    if not m2: return None
    return ('real', 'rp', text[pos + m2.end():])

  kind = ''
  m2 = re.match(r'\s*\(', text[pos:])
  if m2:
    lp = pos + m2.end() - 1
    rp = match_paren(text, lp)
    kind = text[lp + 1:rp - 1].strip()
    # Normalize "kind=rp" / "len=40" forms.
    kind = re.sub(r'^\s*(kind|len)\s*=\s*', '', kind, flags=re.I).strip()
    pos = rp
  return (word, kind, text[pos:])


def parse_entity(text, base_type, kind):
  """
  Parse one entity from a declaration entity-list, eg "mat(6,6) = 0" or "p => null()".
  Returns (name, dims, init).
  """
  init = ''
  # Split off the initializer. "=>" must be checked before "=".
  ix = None
  depth = 0
  quote = None
  for i, ch in enumerate(text):
    if quote:
      if ch == quote: quote = None
      continue
    if ch in ('"', "'"):
      quote = ch
      continue
    if ch in '([': depth += 1
    elif ch in ')]': depth -= 1
    elif ch == '=' and depth == 0:
      ix = i
      break
  if ix is not None:
    init = text[ix:].lstrip('=>').lstrip('=').strip()
    text = text[:ix].strip()

  dims = []
  m = re.match(r'^\s*([a-z]\w*)\s*(\(.*\))?\s*$', text, re.I)
  if not m:
    raise ValueError(f'Cannot parse entity: {text!r}')
  name = m.group(1).lower()
  if m.group(2):
    dims = split_top_level(m.group(2)[1:-1])
  return (name, dims, init)


def parse_declaration(code, comment):
  """Parse one declaration line into a list of Components, or [] if not a declaration."""
  spec = parse_type_spec(code)
  if not spec: return []
  base_type, kind, rest = spec

  allocatable = False
  pointer = False

  # Consume attribute list up to "::" if present.
  if '::' in rest:
    attr_text, _, entity_text = rest.partition('::')
    for attr in split_top_level(attr_text):
      a = attr.strip().lower()
      if a == 'allocatable': allocatable = True
      elif a == 'pointer': pointer = True
      elif a.startswith('dimension'):
        pass    # Handled below via explicit per-entity dims; not used in this codebase.
      elif a in ('', 'private', 'public', 'target', 'save', 'contiguous'): pass
      else: pass
  else:
    # Old style: no "::" and hence no attributes.
    if rest.lstrip().startswith(','): return []
    entity_text = rest

  entity_text = entity_text.strip()
  if not entity_text: return []

  comps = []
  for ent in split_top_level(entity_text):
    if not ent: continue
    name, dims, init = parse_entity(ent, base_type, kind)
    comps.append(Component(name, base_type, kind, dims, allocatable, pointer, init, comment))
  return comps


def parse_file(file_name, structs=None):
  """Parse all derived type definitions in a file. Returns dict name -> Struct."""
  if structs is None: structs = {}
  with open(file_name, errors='replace') as f:
    lines = join_continuations(f.readlines())

  cur = None
  in_contains = False
  for code, comment in lines:
    if cur is None:
      m = RE_TYPE_START.match(code)
      if m:
        cur = Struct(m.group(1).lower(), file_name)
        in_contains = False
      continue
    if RE_TYPE_END.match(code):
      structs[cur.name] = cur
      cur = None
      continue
    if RE_CONTAINS.match(code):
      in_contains = True
      continue
    if in_contains or not code.strip(): continue
    try:
      cur.components.extend(parse_declaration(code, comment))
    except ValueError as e:
      raise ValueError(f'{file_name}: in type {cur.name}: {e}')
  return structs


def parse_files(file_names):
  structs = {}
  for fn in file_names:
    parse_file(fn, structs)
  return structs


if __name__ == '__main__':
  import sys
  st = parse_files(sys.argv[1:])
  print(f'Parsed {len(st)} structs.')
  for name in sorted(st):
    print(f'  {name}: {len(st[name].components)} components')
