!+
! Module tao_completion_mod
!
! Tab completion engine for Tao.
!
! The engine (tao_complete) is shared by two front ends:
!   - The interactive prompt: GNU readline calls tao_rl_complete_c (registered
!     via tao_register_completion / readline_set_completion_fn in sim_utils).
!   - The "pipe complete" command (see tao_pipe_cmd) used by PyTao and other
!     external programs.
!
! Tokens break on blanks and tabs only (matching readline_completion.c) and each
! candidate is a full replacement for the token being completed. At the prompt
! that break set also governs readline's file name fallback, so a file name glued
! to shell syntax ("spawn cat <./di") is taken as one word; quoted names work.
!-

module tao_completion_mod

use tao_struct
use tao_command_names_mod
use tao_input_struct, only: tao_plot_page_input
use attribute_mod, only: attribute_info, ele_attribute_struct, attribute_type, attribute_index, attribute_free
use bmad_routine_interface, only: pointer_to_attribute
use geodesic_lm, only: geodesic_lm_param_struct
use opti_de_mod, only: opti_de_param
use, intrinsic :: iso_c_binding

implicit none

integer, parameter, private :: max_matches$ = 200
integer, parameter, private :: max_ele_scan$ = 100   ! Elements examined for attribute intersection/value union.

contains

!------------------------------------------------------------------------------
!+
! Subroutine tao_complete (line, cursor, word_start, context, matches, common_prefix)
!
! Compute completion candidates for the whitespace-delimited token ending at
! line(cursor-1:cursor-1). Only text to the left of the cursor is considered.
!
! This routine must not do any terminal I/O (out_io, print, etc.): on the
! interactive path it runs inside readline() while the prompt is being edited.
!
! Input:
!   line        -- character(*): Command line being typed.
!   cursor      -- integer: 1-based cursor position. The token being completed
!                    ends at cursor-1. Use len_trim(line)+1 for end-of-line.
!
! Output:
!   word_start  -- integer: 1-based index in line of the start of the token
!                    being completed.
!   context     -- character(*): 'LIST' = matches(:) holds the candidates (possibly
!                    none), 'FILE' = token is a file path, 'NONE' = command not recognized.
!   matches(:)  -- character(100), allocatable: Candidates, at most max_matches$.
!   common_prefix -- character(*), optional: Text to replace the token with right away:
!                    the sole candidate when there is just one, otherwise the token as
!                    typed followed by whatever all candidates (including any beyond
!                    the max_matches$ cap) agree on. Never shorter than the token.
!-

subroutine tao_complete (line, cursor, word_start, context, matches, common_prefix)

character(*), intent(in) :: line
integer, intent(in) :: cursor
integer, intent(out) :: word_start
character(*), intent(out) :: context
character(100), allocatable, intent(out) :: matches(:)
character(*), optional, intent(out) :: common_prefix

character(100) cand(max_matches$)
character(100) token, lcp, lcp_fold, prepend
character(60) attrib_name
character(40) words(8), cmd_name, sub_name
character(2), parameter :: blanks = ' ' // achar(9)
character(1) quote
integer n_end, ix_semi, n_words, n_cand, n_lcp_hits, i, j, ix

!

context = 'NONE'
word_start = 1
n_cand = 0
n_lcp_hits = 0
lcp = ''
lcp_fold = ''
prepend = ''

if (.not. s%initialized) then
  call finish()
  return
endif

n_end = min(cursor-1, len(line))

! Commands can be chained with ";". Work on the text after the last
! semicolon that is not inside a quoted string (cf. check_for_multi_commands).

ix_semi = 0
quote = ''
do i = 1, n_end
  if (quote == '') then
    select case (line(i:i))
    case (';');       ix_semi = i
    case ("'", '"');  quote = line(i:i)
    end select
  else
    if (line(i:i) == quote) quote = ''
  endif
enddo

! The token is the trailing (possibly empty) word. The words before it, minus any
! "-switch" words, give the context.

word_start = ix_semi + 1 + scan(line(ix_semi+1:n_end), blanks, back = .true.)
token = ''
if (word_start <= n_end) token = line(word_start:n_end)

n_words = 0
i = ix_semi + 1
do
  ix = verify(line(i:word_start-1), blanks)
  if (ix == 0) exit
  i = i + ix - 1
  j = scan(line(i:word_start-1), blanks)
  if (j == 0) then
    j = word_start
  else
    j = i + j - 1
  endif
  if ((line(i:i) /= '-' .or. n_words == 0) .and. n_words < size(words)) then
    n_words = n_words + 1
    words(n_words) = line(i:j-1)
  endif
  i = j
enddo

!------------------------------------------
! First word: command names plus user defined aliases. Case sensitive.

if (n_words == 0) then
  if (token(1:1) /= '-') then
    call add_prefix_matches (tao_visible_command_names, .true.)
    do i = 1, s%com%n_alias
      call add_match_if_prefix (s%com%alias(i)%name, .true.)
    enddo
  endif
  context = 'LIST'
  call finish()
  return
endif

! Later words: context depends on the command.

call match_word (words(1), tao_command_names, ix, .true., matched_name = cmd_name)
if (ix <= 0) then
  call finish()
  return
endif

! From here on the command is known, so an unhandled position is an empty list
! rather than NONE: file names are offered only where a command takes a file.
! The add_* helpers below set context = 'LIST' themselves.

context = 'LIST'

! A token starting with "-" is a switch: complete it from the switch context table.

if (token(1:1) == '-') then
  call add_switch_matches ()
  call finish()
  return
endif

select case (cmd_name)

case ('show')
  if (n_words == 1) then
    call add_prefix_matches (tao_show_what_names, .false.)
    ! Pseudo names remapped at the top of tao_show_this.
    call add_prefix_matches ([character(20):: 'plot_page', 'bmad_com', 'ptc_com', &
                                              'space_charge_com', 'floor_plan'], .false.)
  elseif (n_words == 2) then
    call match_word (words(2), tao_show_what_names, ix, matched_name = sub_name)
    select case (sub_name)
    case ('element');            call add_element_matches ()
    case ('data');               call add_data_matches ()
    case ('variables');          call add_var_matches ()
    case ('plot', 'graph', 'curve'); call add_plot_matches (.true., .true.)
    end select
  endif

case ('set')
  sub_name = ''
  if (n_words >= 2) call match_word (words(2), tao_set_target_names, ix, .true., matched_name = sub_name)

  if (at_set_value(attrib_name)) then
    select case (sub_name)
    case ('element')
      if (n_words >= 3) call add_attribute_value_matches (words(3), attrib_name)
    case ('global', 'beam_init', 'bmad_com', 'space_charge_com', 'geodesic_lm', &
          'opti_de_param', 'plot_page', 'ptc_com')
      call add_struct_value_matches (sub_name, attrib_name)
    end select

  elseif (n_words == 1) then
    call add_prefix_matches (tao_set_target_names, .true.)
  elseif (n_words == 2) then
    select case (sub_name)
    case ('global', 'beam_init', 'bmad_com', 'space_charge_com', 'geodesic_lm', &
          'opti_de_param', 'plot_page')
      call add_set_struct_matches (sub_name)
    case ('ptc_com');  call add_prefix_matches (tao_set_ptc_com_names, .false.)
    case ('beam');     call add_prefix_matches (tao_set_beam_names, .false.)
    case ('element')
      call add_element_matches ()
    end select

  elseif (n_words == 3) then
    if (sub_name == 'element') call add_attribute_matches (words(3))
  endif

case ('pipe', 'python')
  if (n_words == 1) call add_prefix_matches (tao_pipe_cmd_names, .false.)

case ('place')
  if (n_words == 1) then
    call add_plot_matches (.true., .false.)
  elseif (n_words == 2) then
    call add_plot_matches (.false., .true.)
  endif

case ('help')
  if (n_words == 1) then
    call add_prefix_matches (tao_visible_command_names, .false.)
  elseif (n_words == 2) then
    call match_word (words(2), [character(8):: 'pipe', 'python'], ix, matched_name = sub_name)
    if (ix > 0) call add_prefix_matches (tao_pipe_cmd_names, .false.)
  endif

case ('change')
  if (n_words == 1) then
    call add_prefix_matches (tao_change_what_names, .false.)
  elseif (n_words == 2) then
    if (words(2) /= '' .and. index('element', trim(words(2))) == 1) call add_element_matches ()
  endif

case ('use', 'veto', 'restore')
  if (n_words == 1) then
    call add_prefix_matches ([character(8):: 'data', 'variable'], .true.)
  else
    call match_word (words(2), [character(8):: 'data', 'variable'], ix, .true., matched_name = sub_name)
    select case (sub_name)
    case ('data');      call add_data_matches ()
    case ('variable');  call add_var_matches ()
    end select
  endif

case ('call')
  if (n_words == 1) context = 'FILE'

case ('read')
  if (n_words == 1) then
    call add_prefix_matches (tao_read_what_names, .false.)
  else
    context = 'FILE'
  endif

case ('ls', 'spawn')
  context = 'FILE'

case ('write')
  if (n_words == 1) then
    call add_prefix_matches (tao_write_action_names, .true.)
  else
    context = 'FILE'
  endif

end select

call finish()

!------------------------------------------
contains

subroutine finish ()
integer lt
matches = cand(1:n_cand)
if (.not. present(common_prefix)) return
! Readline replaces the whole token with this, so it must never be shorter than
! what was typed, and the typed characters keep the user's case. A sole candidate
! is returned as is so readline sees it as a single match.
lt = len_trim(token)
if (n_cand == 1) then
  common_prefix = cand(1)
elseif (len_trim(lcp) > lt) then
  common_prefix = trim(prepend) // token(1:lt) // lcp(lt+1:)
else
  common_prefix = trim(prepend) // token
endif
end subroutine finish

!.................................

subroutine add_prefix_matches (names, exact_case)

character(*) names(:)
logical exact_case
integer in

context = 'LIST'
do in = 1, size(names)
  call add_match_if_prefix (names(in), exact_case)
enddo

end subroutine add_prefix_matches

!.................................

subroutine add_match_if_prefix (name, exact_case)

character(*) name
logical exact_case
character(100) n1, t1, full
integer im, lt

if (name == '') return
lt = len_trim(token)
if (lt > len_trim(name)) return

n1 = name
if (.not. exact_case) call str_upcase (n1, n1)

if (lt > 0) then
  t1 = token
  if (.not. exact_case) call str_upcase (t1, t1)
  if (n1(1:lt) /= t1(1:lt)) return
endif

! The common prefix is compared the same way the token is (case folded unless
! exact_case) but spelled as in the first hit. It counts hits beyond the candidate cap too.

if (n_lcp_hits == 0) then
  lcp = name
  lcp_fold = n1
else
  do im = 1, min(len_trim(lcp_fold), len_trim(n1))
    if (lcp_fold(im:im) /= n1(im:im)) exit
  enddo
  lcp = lcp(1:im-1)
  lcp_fold = lcp_fold(1:im-1)
endif
n_lcp_hits = n_lcp_hits + 1

if (n_cand >= max_matches$) return
full = trim(prepend) // name
do im = 1, n_cand
  if (cand(im) == full) return
enddo
n_cand = n_cand + 1
cand(n_cand) = full

end subroutine add_match_if_prefix

!.................................

subroutine add_switch_matches ()

character(tao_switch_name_len) key
character(20) what_name
character(tao_switch_name_len), allocatable :: sw(:)
integer ik

context = 'LIST'
key = cmd_name

if (cmd_name == 'place' .and. n_words > 1) return   ! -no_buffer must come first.

if (cmd_name == 'show' .and. n_words > 1) then
  call match_word (words(2), tao_show_what_names, ik, matched_name = what_name)
  if (ik <= 0) return
  key = 'show ' // trim(what_name)
endif

sw = tao_switches_for(key)
do ik = 1, size(sw)
  call add_match_if_prefix (sw(ik), .true.)
enddo

end subroutine add_switch_matches

!.................................
! Plot region and/or template names, as used by "place" and plot name arguments.

subroutine add_plot_matches (do_regions, do_templates)

logical do_regions, do_templates
integer ip

context = 'LIST'

if (do_regions .and. allocated(s%plot_page%region)) then
  do ip = 1, size(s%plot_page%region)
    if (s%plot_page%region(ip)%name == '') cycle
    call add_match_if_prefix (s%plot_page%region(ip)%name, .false.)
  enddo
endif

if (do_templates .and. allocated(s%plot_page%template)) then
  do ip = 1, size(s%plot_page%template)
    if (s%plot_page%template(ip)%phantom) cycle
    if (s%plot_page%template(ip)%name == '' .or. s%plot_page%template(ip)%name == 'scratch') cycle
    call add_match_if_prefix (s%plot_page%template(ip)%name, .false.)
  enddo
endif

end subroutine add_plot_matches

!.................................
! "set <struct> <component> = <value>" works via a namelist read (see tao_set_mod),
! so a namelist write of the struct lists exactly the component names set accepts.

subroutine add_set_struct_matches (set_word)

character(*) set_word
character(300) nml_line
character(60) comp_name
integer iu_nml, ios
logical ok

context = 'LIST'

iu_nml = struct_namelist_unit (set_word, ok)
if (.not. ok) return
do
  read (iu_nml, '(a)', iostat = ios) nml_line
  if (ios /= 0) exit
  call parse_namelist_line (nml_line, comp_name)
  if (comp_name /= '') call add_match_if_prefix (comp_name, .false.)
enddo
close (iu_nml)

! Components handled as special cases before the namelist read.

select case (set_word)
case ('global');    call add_prefix_matches ([character(16):: 'phase_units', 'quiet'], .false.)
case ('plot_page'); call add_prefix_matches ([character(16):: 'title', 'subtitle'], .false.)
end select

end subroutine add_set_struct_matches

!.................................
! Namelist dump of a "set" struct in a rewound scratch file.

function struct_namelist_unit (set_word, ok) result (iu_nml)

type (tao_global_struct) global
type (beam_init_struct) beam_init
type (bmad_common_struct) this_bmad_com
type (space_charge_common_struct) this_space_charge_com
type (geodesic_lm_param_struct) this_geodesic_lm
type (tao_plot_page_input) plot_page

character(*) set_word
integer iu_nml, ios
logical ok

namelist / nml_global / global
namelist / nml_beam_init / beam_init
namelist / nml_bmad_com / this_bmad_com
namelist / nml_space_charge_com / this_space_charge_com
namelist / nml_geodesic_lm / this_geodesic_lm
namelist / nml_opti_de_param / opti_de_param
namelist / nml_plot_page / plot_page

open (newunit = iu_nml, status = 'scratch', iostat = ios)
ok = (ios == 0)
if (.not. ok) return

select case (set_word)
case ('global');            write (iu_nml, nml = nml_global, iostat = ios)
case ('beam_init');         write (iu_nml, nml = nml_beam_init, iostat = ios)
case ('bmad_com');          write (iu_nml, nml = nml_bmad_com, iostat = ios)
case ('space_charge_com');  write (iu_nml, nml = nml_space_charge_com, iostat = ios)
case ('geodesic_lm');       write (iu_nml, nml = nml_geodesic_lm, iostat = ios)
case ('opti_de_param');     write (iu_nml, nml = nml_opti_de_param, iostat = ios)
case ('plot_page');         write (iu_nml, nml = nml_plot_page, iostat = ios)
end select

rewind (iu_nml)

end function struct_namelist_unit

!.................................
! Component name (downcased, '' for a continuation line) and value from one line
! of a namelist dump: "<struct>%<name>= value," (gfortran) or "%<name> = value," (ifort).

subroutine parse_namelist_line (nml_line, comp_name, value_str)

character(*) nml_line, comp_name
character(*), optional :: value_str
character(*), parameter :: name_chars = &
          'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_%'
character(len(nml_line)) rest
integer ix1, ix2, ieq, iend

comp_name = ''
if (present(value_str)) value_str = ''

nml_line = adjustl(nml_line)
ix1 = index(nml_line, '%')
if (ix1 == 0) return
ix2 = ix1 + 1
do while (ix2 <= len(nml_line))
  if (verify(nml_line(ix2:ix2), name_chars) /= 0) exit
  ix2 = ix2 + 1
enddo
if (ix2 == ix1 + 1) return
ieq = index(nml_line(ix2:), '=')
if (ieq == 0) return

comp_name = downcase(nml_line(ix1+1:ix2-1))

if (present(value_str)) then
  rest = adjustl(nml_line(ix2+ieq:))
  iend = scan(rest, ', ')
  if (iend == 0) iend = len_trim(rest) + 1
  value_str = rest(1:iend-1)
  if (value_str(1:1) == '"' .or. value_str(1:1) == "'") value_str = value_str(2:len_trim(value_str)-1)
endif

end subroutine parse_namelist_line

!.................................

function struct_component_value (set_word, comp) result (value_str)

character(*) set_word, comp
character(100) value_str
character(300) nml_line
character(60) comp_name
character(60) comp_lc
integer iu_nml, ios
logical ok

value_str = ''
comp_lc = downcase(comp)

iu_nml = struct_namelist_unit (set_word, ok)
if (.not. ok) return
do
  read (iu_nml, '(a)', iostat = ios) nml_line
  if (ios /= 0) exit
  call parse_namelist_line (nml_line, comp_name, value_str)
  if (comp_name == comp_lc) exit
  value_str = ''
enddo
close (iu_nml)

end function struct_component_value

!.................................
! True if the token is the value in "set <target> ... <attrib> = <value>".
! Handles "attrib = val", "attrib =val" and the glued "attrib=val". On return
! attrib holds the attribute or component name, token holds only the value
! prefix, and prepend holds any "attrib=" text that candidates must keep in front.

function at_set_value (attrib) result (is_value)

character(*) attrib
logical is_value
integer ie, lw

is_value = .false.
attrib = ''

ie = index(token, '=')
if (ie > 0) then
  if (ie > 1) then
    attrib = token(1:ie-1)
  elseif (n_words >= 1) then
    attrib = words(n_words)
  endif
  prepend = token(1:ie)
  token = token(ie+1:)
  is_value = .true.

elseif (n_words >= 3) then
  lw = len_trim(words(n_words))
  if (words(n_words) == '=') then
    attrib = words(n_words-1)
    is_value = .true.
  elseif (lw > 1 .and. words(n_words)(lw:lw) == '=') then
    attrib = words(n_words)(1:lw-1)
    is_value = .true.
  endif
endif

end function at_set_value

!.................................
! Universe named by an optional "n@" prefix, which is stripped from sel. Null if
! the index is bad.

function selector_universe (sel) result (u)

type (tao_universe_struct), pointer :: u
character(*) sel
integer ia, iu, ios

nullify (u)
iu = s%global%default_universe
ia = index(sel, '@')
if (ia > 1) then
  read (sel(1:ia-1), *, iostat = ios) iu
  if (ios /= 0) return
  sel = sel(ia+1:)
endif
if (.not. allocated(s%u)) return
if (iu < lbound(s%u, 1) .or. iu > ubound(s%u, 1)) return
u => s%u(iu)

end function selector_universe

!.................................
! Elements matching a Tao selector ("quad::*", "1:10", "2@q1", ...). lat_ele_locator
! errors are expected for a half-typed selector and must not print inside readline.

subroutine locate_elements (selector, eles, n_loc)

type (ele_pointer_struct), allocatable :: eles(:)
type (out_io_output_direct_struct) out_state
type (tao_universe_struct), pointer :: u

character(*) selector
character(len(selector)) sel
integer n_loc
logical err

n_loc = 0
sel = selector
u => selector_universe(sel)
if (.not. associated(u) .or. sel == '') return

call output_direct (get = out_state)
call output_direct (print_and_capture = .false.)
call lat_ele_locator (sel, u%model%lat, eles, n_loc, err)
call output_direct (set = out_state)
if (err) n_loc = 0

end subroutine locate_elements

!.................................
! The beginning element and control lords (overlays, groups, girders, rampers)
! have attribute tables unlike any real element. When a selector matches several
! elements they are left out of the attribute intersection and value union; a
! single match is always used.

function has_fixed_attributes (ele) result (is_fixed)
type (ele_struct) ele
logical is_fixed
select case (ele%key)
case (beginning_ele$, overlay$, group$, girder$, ramper$); is_fixed = .false.
case default;                                              is_fixed = .true.
end select
end function has_fixed_attributes

! Index of the first matched element to use as the attribute reference: for a
! single match that element itself, otherwise the first with fixed attributes.

function first_settable (eles, n_loc) result (ie)
type (ele_pointer_struct) eles(:)
integer n_loc, ie
if (n_loc == 1) then
  ie = 1
  return
endif
do ie = 1, min(n_loc, max_ele_scan$)
  if (has_fixed_attributes(eles(ie)%ele)) return
enddo
ie = 0
end function first_settable

!.................................
! Values for "set element <name> <attrib> = <value>": the values valid for this
! element's switch attributes (via tao_enum_value_names) or T/F for logicals.

subroutine add_attribute_value_matches (ele_name, attrib)

type (ele_pointer_struct), allocatable :: eles(:)
character(*) ele_name, attrib
character(40), allocatable :: names(:)
integer, allocatable :: ixs(:)
integer n_loc, ie, ie0, in

context = 'LIST'
call locate_elements (ele_name, eles, n_loc)
if (n_loc == 0) return

ie0 = first_settable(eles, n_loc)
if (ie0 == 0) return

select case (attribute_type(upcase(attrib), eles(ie0)%ele))
case (is_logical$)
  call add_prefix_matches ([character(8):: 'T', 'F'], .false.)
case (is_switch$)
  ! Union of the values valid for the matched elements. The reference element is
  ! always included: with a single match it may be a control lord.
  do ie = ie0, min(n_loc, max_ele_scan$)
    if (ie /= ie0 .and. .not. has_fixed_attributes(eles(ie)%ele)) cycle
    call tao_enum_value_names (attrib, names, ixs, eles(ie)%ele)
    if (.not. allocated(names)) cycle
    do in = 1, size(names)
      call add_match_if_prefix (downcase(names(in)), .false.)
    enddo
  enddo
end select

end subroutine add_attribute_value_matches

!.................................
! Values for "set <struct> <component> = <value>": enumerated components share
! their names with "pipe enum" (track_type, optimizer, ...); otherwise a component
! whose namelist dump value is T or F is a logical.

subroutine add_struct_value_matches (set_word, comp)

character(*) set_word, comp
character(40), allocatable :: names(:)
integer, allocatable :: ixs(:)
character(100) value_str
integer in

context = 'LIST'

call tao_enum_value_names (downcase(comp), names, ixs, switch_attribs = .false.)
if (allocated(names)) then
  do in = 1, size(names)
    call add_match_if_prefix (names(in), .false.)
  enddo
  return
endif

! ptc_common_struct has pointer components so there is no namelist dump to consult.

if (set_word == 'ptc_com') then
  if (any(tao_set_ptc_com_logical_names == downcase(comp))) call add_prefix_matches ([character(8):: 'T', 'F'], .false.)
  return
endif

value_str = struct_component_value (set_word, comp)
if (value_str == 'T' .or. value_str == 'F') call add_prefix_matches ([character(8):: 'T', 'F'], .false.)

end subroutine add_struct_value_matches

!.................................
! Attribute names for "set element <name> <attrib>", from bmad's attribute table
! for the first element matching <name> (wildcards allowed) in the default universe.

subroutine add_attribute_matches (ele_name)

type (ele_pointer_struct), allocatable :: eles(:)
type (ele_attribute_struct) attrib
type (all_pointer_struct) a_ptr

character(*) ele_name
character(40) attrib_names(num_ele_attrib_extended$)
integer n_loc, n_attr, ia, ie, ie0, i2
logical err

context = 'LIST'
call locate_elements (ele_name, eles, n_loc)
if (n_loc == 0) return

! Attributes of the first matched element (the extended range includes switches
! and logicals), intersected with every other matched element since set applies to all.

ie0 = first_settable(eles, n_loc)
if (ie0 == 0) return
ie = ie0

n_attr = 0
do ia = 1, num_ele_attrib_extended$
  attrib = attribute_info(eles(ie)%ele, ia)
  if (attrib%name == '' .or. attrib%name == null_name$ .or. attrib%name(1:1) == '!') cycle
  if (attrib%state == does_not_exist$ .or. attrib%state == private$) cycle
  n_attr = n_attr + 1
  attrib_names(n_attr) = attrib%name
enddo

do ie = ie+1, min(n_loc, max_ele_scan$)
  if (.not. has_fixed_attributes(eles(ie)%ele)) cycle
  i2 = 0
  do ia = 1, n_attr
    if (attribute_index(eles(ie)%ele, attrib_names(ia), print_error = .false.) == 0) cycle
    i2 = i2 + 1
    attrib_names(i2) = attrib_names(ia)
  enddo
  n_attr = i2
enddo

! Offer only what set_ele_attribute can set: resolvable by pointer_to_attribute and
! free (dependent attributes allowed, so b1_gradient counts even with field_master).

do ia = 1, n_attr
  call pointer_to_attribute (eles(ie0)%ele, attrib_names(ia), .false., a_ptr, err, err_print_flag = .false.)
  if (err) cycle
  if (.not. attribute_free (eles(ie0)%ele, attrib_names(ia), .false., dependent_attribs_free = .true.)) cycle
  call add_match_if_prefix (downcase(attrib_names(ia)), .false.)
enddo

end subroutine add_attribute_matches

!.................................
! Element names, plus the "n@" and "key::" selector prefixes: "quad::Q0" completes
! to quadrupoles starting with Q0, and without a "::" the element types present in
! the lattice ("quadrupole::", ...) are offered too.

subroutine add_element_matches ()

type (tao_universe_struct), pointer :: u
type (branch_struct), pointer :: branch
character(100) sel
integer ib, ie, ia, ic, ik, ix_key
logical key_present(n_key$)

context = 'LIST'

! prepend // token must always equal the text typed, so token tracks sel as
! each prefix moves to prepend.

sel = token
ia = index(sel, '@')
u => selector_universe(sel)
if (.not. associated(u)) return
if (ia > 1) then
  prepend = trim(prepend) // token(1:ia)
  token = sel
endif

ix_key = 0
ic = index(sel, '::')
if (ic > 0) then
  if (ic == 1) return
  call match_word (sel(1:ic-1), key_name, ix_key)
  if (ix_key <= 0) return
  prepend = trim(prepend) // sel(1:ic+1)
  sel = sel(ic+2:)
  token = sel
endif

! Element types go first so a big lattice cannot crowd them out of the capped list.

if (ix_key == 0) then
  key_present = .false.
  do ib = 0, ubound(u%model%lat%branch, 1)
    branch => u%model%lat%branch(ib)
    do ie = 1, branch%n_ele_max
      key_present(branch%ele(ie)%key) = .true.
    enddo
  enddo
  do ik = 1, size(key_name)
    if (.not. key_present(ik) .or. key_name(ik)(1:1) == '!') cycle
    call add_match_if_prefix (trim(downcase(key_name(ik))) // '::', .false.)
  enddo
endif

do ib = 0, ubound(u%model%lat%branch, 1)
  branch => u%model%lat%branch(ib)
  do ie = 1, branch%n_ele_max
    if (ix_key > 0 .and. branch%ele(ie)%key /= ix_key) cycle
    call add_match_if_prefix (branch%ele(ie)%name, .false.)
  enddo
enddo

end subroutine add_element_matches

!.................................

subroutine add_data_matches ()

type (tao_universe_struct), pointer :: u
type (tao_d2_data_struct), pointer :: d2
character(100) name
integer iu, id, id1

context = 'LIST'
if (.not. allocated(s%u)) return
iu = s%global%default_universe
if (iu < lbound(s%u, 1) .or. iu > ubound(s%u, 1)) return
u => s%u(iu)

do id = 1, u%n_d2_data_used
  d2 => u%d2_data(id)
  if (d2%name == '') cycle
  call add_match_if_prefix (d2%name, .true.)
  if (.not. allocated(d2%d1)) cycle
  do id1 = 1, size(d2%d1)
    name = trim(d2%name) // '.' // d2%d1(id1)%name
    call add_match_if_prefix (name, .true.)
  enddo
enddo

end subroutine add_data_matches

!.................................

subroutine add_var_matches ()

integer iv

context = 'LIST'
if (.not. allocated(s%v1_var)) return
do iv = 1, s%n_v1_var_used
  call add_match_if_prefix (s%v1_var(iv)%name, .true.)
enddo

end subroutine add_var_matches

end subroutine tao_complete

!------------------------------------------------------------------------------
!+
! Function tao_rl_complete_c (line_c, point, istart, iend, buf_c, buf_size) result (n_cand)
!
! bind(c) completion callback invoked by GNU readline via the shim in
! sim_utils/io/readline_completion.c. Registered by tao_register_completion.
!
! Input:
!   line_c    -- type(c_ptr): Null terminated readline line buffer.
!   point     -- integer(c_int): 0-based cursor offset (rl_point).
!   istart    -- integer(c_int): 0-based start of the word readline is completing.
!   iend      -- integer(c_int): 0-based end of that word (unused).
!   buf_size  -- integer(c_int): Size of the buffer at buf_c.
!
! Output:
!   buf_c     -- type(c_ptr): Null terminated buffer. Line 1 is the common prefix
!                  readline should insert; each following line is one candidate.
!   n_cand    -- integer(c_int): Number of candidates (0 = recognized position with
!                  nothing matching), or -1 = use readline's file name completion
!                  (the engine's FILE and NONE contexts).
!-

function tao_rl_complete_c (line_c, point, istart, iend, buf_c, buf_size) bind(c) result (n_cand)

type(c_ptr), value :: line_c, buf_c
integer(c_int), value :: point, istart, iend, buf_size
integer(c_int) :: n_cand

character(kind=c_char), pointer :: line_p(:), buf_p(:)
character(4000) line_f
character(100) lcp
character(8) context
character(100), allocatable :: matches(:)
integer word_start, i, k, n

!

n_cand = 0
if (.not. s%initialized) return
if (.not. c_associated(line_c) .or. .not. c_associated(buf_c)) return

! Only the text left of the cursor matters, and rl_point is exactly its length.

n = min(int(point), len(line_f))
line_f = ''
if (n > 0) then
  call c_f_pointer (line_c, line_p, [n])
  do i = 1, n
    line_f(i:i) = line_p(i)
  enddo
endif

call tao_complete (line_f, n + 1, word_start, context, matches, lcp)

if (context /= 'LIST') then
  n_cand = -1
  return
endif

! The engine and readline must agree on the token span. If they do not
! (for example a ";" inside the token), do not offer anything.

if (word_start - 1 /= istart) return

call c_f_pointer (buf_c, buf_p, [buf_size])
k = 0
if (append_line(lcp)) then
  do i = 1, size(matches)
    if (.not. append_line(matches(i))) exit
    n_cand = n_cand + 1
  enddo
endif
buf_p(k+1) = c_null_char

!------------------------------------------
contains

! Append str plus a newline, always leaving room for the final null.

function append_line (str) result (ok)
character(*) str
logical ok
integer jj, lt
lt = len_trim(str)
ok = (k + lt + 2 <= buf_size)
if (.not. ok) return
do jj = 1, lt
  k = k + 1
  buf_p(k) = str(jj:jj)
enddo
k = k + 1
buf_p(k) = c_new_line
end function append_line

end function tao_rl_complete_c

!------------------------------------------------------------------------------
!+
! Subroutine tao_register_completion ()
!
! Install tao_rl_complete_c as the readline tab completion callback and name the
! application "Tao" for inputrc "$if Tao" blocks. Idempotent.
!-

subroutine tao_register_completion ()

interface
  subroutine readline_set_completion_fn (fn) bind(c, name = 'readline_set_completion_fn')
    import :: c_funptr
    type(c_funptr), value :: fn
  end subroutine
  subroutine readline_set_app_name (name) bind(c, name = 'readline_set_app_name')
    import :: c_char
    character(kind=c_char) :: name(*)
  end subroutine
end interface

logical, save :: registered = .false.

!

if (registered) return
registered = .true.
! Lets users put Tao-only readline settings in ~/.inputrc inside "$if Tao ... $endif".
call readline_set_app_name ('Tao' // c_null_char)
call readline_set_completion_fn (c_funloc(tao_rl_complete_c))

end subroutine tao_register_completion

end module tao_completion_mod
