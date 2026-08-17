#!/usr/bin/env python3
"""
Generate tao_attrib_resolve_mod.f90: recursive resolvers mapping a structure component name
string onto a pointer to that component.

Run from the root of the bmad-ecosystem repo:

    python3 tao/scripts/generate_attrib_tables.py

The output file is checked into the repo.  Regenerating must be a no-op unless a structure
definition changed; the tao_attrib_test regression test enforces this.

Design notes:
  - Only components reachable from the root structures without traversing a `pointer` component
    are emitted.  Following pointers pulls in the whole of bmad plus the PTC types (fibre,
    layout) and is meaningless for setting values from a string in any case.
  - Components that cannot be resolved to a pointer (rank >= 2, complex, or pointer components)
    get an explicit runtime error case rather than being silently omitted, so that a user typing
    such a name gets a real explanation instead of "no such component".
"""

import os
import sys
import textwrap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fortran_struct_parser as fsp

# Files holding the structure definitions, and the module that exports each.

SOURCE_MODULES = [
  ('bmad/modules/bmad_struct.f90', 'bmad_struct'),
  ('tao/code/tao_struct.f90', 'tao_struct'),
  ('tao/code/tao_input_struct.f90', 'tao_input_struct'),
  ('sim_utils/plot/quick_plot_struct.f90', 'quick_plot_struct'),
  ('sim_utils/optimizers/opti_de_mod.f90', 'opti_de_mod'),
  ('sim_utils/geodesic_lm/geodesic_lm.f90', 'geodesic_lm'),
  ('sim_utils/math/spline_mod.f90', 'spline_mod'),
  ('bmad/space_charge/csr_and_space_charge_mod.f90', 'csr_and_space_charge_mod'),
]

# Root structures. Everything reachable from these (not through pointers) gets a resolver.
# The first group is what the `set` command needs, the second is the init file input structures.

ROOT_STRUCTS = [
  'tao_global_struct',
  'bmad_common_struct',
  'space_charge_common_struct',
  'beam_init_struct',
  'tao_plot_page_struct',
  'opti_de_param_struct',
  'geodesic_lm_param_struct',
  'tao_ele_shape_struct',
  'qp_point_struct',
  'aperture_param_struct',

  'tao_d2_data_input',
  'tao_d1_data_input',
  'tao_datum_input',
  'tao_v1_var_input',
  'tao_var_input',
  'tao_region_input',
  'tao_place_input',
  'tao_curve_input',
  'tao_graph_input',
  'tao_plot_input',
  'tao_design_lat_input',
  'tao_key_input',
  'tao_plot_page_input',
  'tao_ele_shape_input',
  'tao_building_wall_point_struct',
  'tao_shape_pattern_point_struct',
  'qp_line_struct',
]

OUT_FILE = 'tao/code/tao_attrib_resolve_mod.f90'

# Map a parsed base type to the tao_ptr_struct scalar and rank 1 component names.

PTR_COMPONENT = {
  'real': ('r', 'r1'),
  'integer': ('i', 'i1'),
  'logical': ('l', 'l1'),
  'character': ('str', 'str1'),
}


def resolver_name(struct_name):
  return 'tao_res_' + struct_name


def compute_closure(structs):
  """Transitive closure of ROOT_STRUCTS over non-pointer derived type components."""
  missing = []
  seen = set()
  stack = list(ROOT_STRUCTS)
  while stack:
    name = stack.pop()
    if name in seen: continue
    if name not in structs:
      missing.append(name)
      continue
    seen.add(name)
    for c in structs[name].components:
      if not c.is_struct or c.pointer: continue
      stack.append(c.kind)
  return sorted(seen), sorted(set(missing))


def emit_case(comp, structs, unsupported):
  """Emit the select-case branch for one component. Returns a list of source lines."""
  name = comp.name
  out = [f"case ('{name}')"]

  # Components that cannot be turned into a pointer.

  reason = None
  if comp.pointer:
    reason = 'IS A POINTER COMPONENT AND CANNOT BE SET'
  elif comp.rank >= 2:
    reason = 'IS A RANK 2 OR HIGHER ARRAY WHICH IS NOT SUPPORTED'
  elif comp.base_type == 'complex':
    reason = 'IS OF TYPE COMPLEX WHICH IS NOT SUPPORTED'
  elif comp.is_struct and comp.kind not in structs:
    reason = 'HAS A TYPE WHOSE DEFINITION WAS NOT FOUND'

  if reason:
    unsupported.append(f'{name}: {reason}')
    out.append(f"  err = .true.; why = 'COMPONENT " + name.upper() + " ' // &")
    out.append(f"          '{reason}'")
    return out

  is_alloc = comp.allocatable

  if comp.is_struct:
    sub = resolver_name(comp.kind)
    if comp.rank == 0:
      out.append('  if (tao_res_struct_bad(head, has_sub, rest, err, why)) return')
      out.append(f'  call {sub} (obj%{name}, rest, ptr, err, why)')
    else:
      if is_alloc:
        out.append(f'  if (tao_res_alloc_bad(head, allocated(obj%{name}), err, why)) return')
        lo, hi = f'lbound(obj%{name},1)', f'ubound(obj%{name},1)'
      else:
        lo, hi = comp.bounds(0)
      out.append(f'  if (tao_res_struct_array_bad(head, rest, has_sub, isub, {lo}, {hi}, err, why)) return')
      out.append(f'  call {sub} (obj%{name}(isub), rest, ptr, err, why)')
    return out

  pscal, parr = PTR_COMPONENT[comp.base_type]

  if comp.rank == 0:
    if is_alloc:
      out.append(f'  if (tao_res_alloc_bad(head, allocated(obj%{name}), err, why)) return')
    out.append('  if (tao_res_scalar_bad(head, has_sub, rest, err, why)) return')
    out.append(f'  ptr%{pscal} => obj%{name}')
    return out

  # Rank 1.

  if is_alloc:
    out.append(f'  if (tao_res_alloc_bad(head, allocated(obj%{name}), err, why)) return')
    lo, hi = f'lbound(obj%{name},1)', f'ubound(obj%{name},1)'
  else:
    lo, hi = comp.bounds(0)
  out.append(f'  if (tao_res_array_bad(head, rest, has_sub, isub, {lo}, {hi}, err, why)) return')
  out.append(f'  if (has_sub) then')
  out.append(f'    ptr%{pscal} => obj%{name}(isub)')
  out.append(f'  else')
  out.append(f'    ptr%{parr} => obj%{name}')
  out.append(f'  endif')
  return out


def flatten(name, structs, stack=()):
  """
  Flattened list of ultimate components of a structure, in declaration order.

  This is the order a namelist read assigns to when given a whole structure value list such as
    ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
  Arrays are expanded element by element and nested structures are expanded in place.

  Returns a list of (access_path, base_type) or None if the structure cannot be flattened,
  which happens when it has an allocatable or deferred shape component. Those are runtime
  structures that are never the target of a positional assignment from an init file.
  """
  if name in stack: return None
  if name not in structs: return None

  out = []
  for c in structs[name].components:
    if c.pointer: continue
    if c.allocatable or c.is_deferred: return None
    if c.base_type == 'complex': return None

    # Element subscripts for this component, in Fortran storage order.
    subs = [None]
    for i in range(c.rank):
      b = c.bounds(i)
      if b is None: return None
      try:
        lo, hi = int(b[0]), int(b[1])
      except ValueError:
        return None
      subs = [(s, j) for s in subs for j in range(lo, hi + 1)] if subs != [None] \
             else [(None, j) for j in range(lo, hi + 1)]

    if c.rank == 0:
      paths = ['%' + c.name]
    else:
      paths = ['%' + c.name + '(' + str(j) + ')' for (_s, j) in subs]

    for p in paths:
      if c.is_struct:
        sub = flatten(c.kind, structs, stack + (name,))
        if sub is None: return None
        for (sp, st) in sub:
          out.append((p + sp, st))
      else:
        out.append((p, c.base_type))

  return out


def emit_slot_resolver(struct, structs):
  """Emit the by-slot resolver used for positional (whole structure) assignment."""
  name = struct.name
  flat = flatten(name, structs)

  lines = []
  lines.append('!--------------------------------------------------------------------------')
  lines.append('!+')
  lines.append(f'! Subroutine {resolver_name(name)}_slot (obj, name, i_slot, ptr, err, why)')
  lines.append('!')
  lines.append('! Pointer to the i_slot-th ultimate component, in declaration order, of the thing')
  lines.append(f'! that name selects within a {name}.')
  lines.append('!')
  lines.append('! This is the order a namelist read uses when something is given a positional value')
  lines.append('! list. A blank name means the whole structure, as in')
  lines.append('!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2')
  lines.append('! and a non blank name selects a component, which may itself be a structure, as in')
  lines.append('!     plot_page%border = 0, 0, 0, 0, "%PAGE"')
  lines.append('!')
  lines.append('! err is set when i_slot is past the last component, which is how the caller')
  lines.append('! detects that too many values were supplied.')
  lines.append('!-')
  lines.append('')
  lines.append(f'recursive subroutine {resolver_name(name)}_slot (obj, name, i_slot, ptr, err, why)')
  lines.append('')
  lines.append(f'type ({name}), target :: obj')
  lines.append('type (tao_ptr_struct) ptr')
  lines.append('')
  lines.append('character(*), intent(in) :: name')
  lines.append('character(*), intent(out) :: why')
  lines.append('character(:), allocatable :: head, rest')
  lines.append('')
  lines.append('integer, intent(in) :: i_slot')
  lines.append('integer isub')
  lines.append('')
  lines.append('logical, intent(out) :: err')
  lines.append('logical has_sub')
  lines.append('')
  lines.append('!')
  lines.append('')
  lines.append('ptr = tao_ptr_struct()')
  lines.append('err = .false.')
  lines.append("why = ''")
  lines.append('')
  lines.append('! A blank name means the whole structure.')
  lines.append('')
  lines.append("if (name == '') then")

  if flat is None:
    lines.append('  ! This structure has an allocatable or deferred shape component so it has no')
    lines.append('  ! fixed component ordering and cannot be assigned positionally.')
    lines.append('  err = .true.')
    lines.append(f"  why = 'STRUCTURE {name.upper()} CANNOT BE SET FROM A POSITIONAL VALUE LIST'")
    lines.append('  return')
  else:
    lines.append('  select case (i_slot)')
    for i, (path, btype) in enumerate(flat, start=1):
      pscal, _ = PTR_COMPONENT[btype]
      lines.append(f'  case ({i});  ptr%{pscal} => obj{path}')
    lines.append('  case default')
    lines.append('    err = .true.')
    lines.append(f"    why = 'TOO MANY VALUES. {name.upper()} HAS {len(flat)} COMPONENTS'")
    lines.append('  end select')
    lines.append('  return')

  lines.append('endif')
  lines.append('')
  lines.append('! A non blank name selects a component to descend into.')
  lines.append('')
  lines.append('call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)')
  lines.append('if (err) return')
  lines.append('')
  lines.append('select case (head)')
  lines.append('')

  seen = set()
  for c in struct.components:
    if c.name in seen: continue
    seen.add(c.name)
    if c.pointer or c.base_type == 'complex' or c.rank >= 2: continue
    if c.is_struct and c.kind not in structs: continue

    lines.append(f"case ('{c.name}')")

    if c.is_struct:
      sub = resolver_name(c.kind) + '_slot'
      if c.rank == 0:
        lines.append('  if (tao_res_struct_bad_for_slot(head, has_sub, err, why)) return')
        lines.append(f'  call {sub} (obj%{c.name}, rest, i_slot, ptr, err, why)')
      else:
        if c.allocatable:
          lines.append(f'  if (tao_res_alloc_bad(head, allocated(obj%{c.name}), err, why)) return')
          lo, hi = f'lbound(obj%{c.name},1)', f'ubound(obj%{c.name},1)'
        else:
          lo, hi = c.bounds(0)
        lines.append(f'  if (tao_res_struct_array_bad_for_slot(head, has_sub, isub, {lo}, {hi}, err, why)) return')
        lines.append(f'  call {sub} (obj%{c.name}(isub), rest, i_slot, ptr, err, why)')

    else:
      pscal, _ = PTR_COMPONENT[c.base_type]
      lines.append('  if (tao_res_slot_leaf_bad(head, rest, err, why)) return')
      if c.rank == 0:
        lines.append('  if (i_slot /= 1) then')
        lines.append('    err = .true.')
        lines.append(f"    why = 'TOO MANY VALUES FOR COMPONENT: ' // head")
        lines.append('    return')
        lines.append('  endif')
        lines.append(f'  ptr%{pscal} => obj%{c.name}')
      else:
        if c.allocatable:
          lines.append(f'  if (tao_res_alloc_bad(head, allocated(obj%{c.name}), err, why)) return')
          lo, hi = f'lbound(obj%{c.name},1)', f'ubound(obj%{c.name},1)'
        else:
          lo, hi = c.bounds(0)
        lines.append(f'  if (i_slot < 1 .or. i_slot > {hi} - {lo} + 1) then')
        lines.append('    err = .true.')
        lines.append(f"    why = 'TOO MANY VALUES FOR COMPONENT: ' // head")
        lines.append('    return')
        lines.append('  endif')
        lines.append(f'  ptr%{pscal} => obj%{c.name}({lo} + i_slot - 1)')

    lines.append('')

  lines.append('case default')
  lines.append('  err = .true.')
  lines.append(f"  why = 'NO SUCH COMPONENT: ' // head // '  (IN {name.upper()})'")
  lines.append('end select')
  lines.append('')
  lines.append(f'end subroutine {resolver_name(name)}_slot')
  lines.append('')
  return lines


def emit_resolver(struct, structs, report):
  name = struct.name
  lines = []
  lines.append('!--------------------------------------------------------------------------')
  lines.append('!--------------------------------------------------------------------------')
  lines.append('!--------------------------------------------------------------------------')
  lines.append('!+')
  lines.append(f'! Subroutine {resolver_name(name)} (obj, name, ptr, err, why)')
  lines.append('!')
  lines.append(f'! Resolve a component name of a {name} instance to a pointer.')
  lines.append(f'! Structure defined in: {struct.source_file}')
  lines.append('!-')
  lines.append('')
  lines.append(f'recursive subroutine {resolver_name(name)} (obj, name, ptr, err, why)')
  lines.append('')
  lines.append(f'type ({name}), target :: obj')
  lines.append('type (tao_ptr_struct) ptr')
  lines.append('')
  lines.append('character(*), intent(in) :: name')
  lines.append('character(*), intent(out) :: why')
  lines.append('character(:), allocatable :: head, rest')
  lines.append('')
  lines.append('integer isub')
  lines.append('')
  lines.append('logical, intent(out) :: err')
  lines.append('logical has_sub')
  lines.append('')
  lines.append('!')
  lines.append('')
  lines.append('! ptr is cleared on entry. Without this a component left associated by an')
  lines.append('! earlier call would win over the one set here, since tao_set_ptr_value takes')
  lines.append('! the first associated component it finds.')
  lines.append('')
  lines.append('ptr = tao_ptr_struct()')
  lines.append('')
  lines.append('call tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)')
  lines.append('if (err) return')
  lines.append('')
  lines.append('select case (head)')
  lines.append('')

  seen_names = set()
  unsupported = []
  for c in struct.components:
    if c.name in seen_names:
      report.append(f'WARNING: duplicate component {name}%{c.name} -- second one ignored')
      continue
    seen_names.add(c.name)
    lines.extend(emit_case(c, structs, unsupported))
    lines.append('')

  for u in unsupported:
    report.append(f'  unsupported {name}%{u}')

  lines.append('case default')
  lines.append("  err = .true.")
  lines.append(f"  why = 'NO SUCH COMPONENT: ' // head // '  (IN {name.upper()})'")
  lines.append('end select')
  lines.append('')
  lines.append(f'end subroutine {resolver_name(name)}')
  lines.append('')
  return lines


def main():
  structs = {}
  for fn, _mod in SOURCE_MODULES:
    if not os.path.exists(fn):
      sys.exit(f'Cannot find struct definition file: {fn}\nRun from the repo root.')
    fsp.parse_file(fn, structs)

  closure, missing = compute_closure(structs)
  if missing:
    sys.exit('Root structures not found: ' + ', '.join(missing))

  report = []
  body = []
  n_slot_total = 0
  for sname in closure:
    body.extend(emit_resolver(structs[sname], structs, report))
    body.extend(emit_slot_resolver(structs[sname], structs))
    flat = flatten(sname, structs)
    if flat is None:
      report.append(f'  no positional assignment for {sname} (allocatable component)')
    else:
      n_slot_total += len(flat)

  n_comp = sum(len(structs[s].components) for s in closure)

  header = [
    '!+',
    '! Module tao_attrib_resolve_mod',
    '!',
    '! *** THIS FILE IS GENERATED. DO NOT EDIT. ***',
    '!',
    '! Regenerate with, from the root of the bmad-ecosystem repo:',
    '!     python3 tao/scripts/generate_attrib_tables.py',
    '!',
    '! Each routine here resolves a structure component name such as',
    '!     "floor_plan%ele_shape(3)%ele_id"',
    '! onto a pointer to that component, replacing the use of a Fortran namelist read on a',
    '! scratch file for the purpose of setting a structure component from a string.',
    '!',
    '! IMPORTANT: the structure passed to a resolver must have the TARGET attribute.',
    '! The obj dummy argument below is declared with TARGET. If the actual argument does not',
    '! also have TARGET then the returned pointer becomes undefined as soon as the resolver',
    '! returns, since the compiler is free to pass a copy. That failure mode is nasty: it can',
    '! appear to work for one structure and segfault for another, and the crash surfaces at',
    '! whatever code runs next rather than at the call. So declare the variable as, eg:',
    '!     type (bmad_common_struct), target :: this_bmad_com',
    '! For a module variable that lacks TARGET, copy it into a local TARGET variable, resolve',
    '! and set through that, then copy back. See tao_set_opti_de_param_cmd for an example.',
    '!',
    f'! Structures covered: {len(closure)}.  Components: {n_comp}.',
    '!-',
    '',
    'module tao_attrib_resolve_mod',
    '',
    'use tao_attrib_ptr_mod',
  ]
  for _fn, mod in SOURCE_MODULES:
    header.append(f'use {mod}')
  header += [
    '',
    'implicit none',
    '',
    'contains',
    '',
  ]

  footer = ['end module tao_attrib_resolve_mod', '']

  text = '\n'.join(header + body + footer)

  # --check regenerates and compares without writing. This is what keeps the generated file
  # from drifting when someone edits a structure definition: see test_generated.py.

  if '--check' in sys.argv:
    if not os.path.exists(OUT_FILE):
      sys.exit(f'{OUT_FILE} does not exist. Run: python3 {sys.argv[0]}')
    with open(OUT_FILE) as f:
      current = f.read()
    if current != text:
      sys.exit(f'{OUT_FILE} is out of date with respect to the structure definitions.\n'
               f'Regenerate with: python3 {sys.argv[0]}')
    print(f'{OUT_FILE} is up to date.')
    return

  with open(OUT_FILE, 'w') as f:
    f.write(text)

  print(f'Wrote {OUT_FILE}: {len(closure)} structures, {n_comp} components, '
        f'{len(text.splitlines())} lines.')
  for r in report:
    print(r)


if __name__ == '__main__':
  main()
