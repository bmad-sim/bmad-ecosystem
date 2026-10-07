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
! Completion tokens break on space/tab only. This must agree with
! rl_completer_word_break_characters set in sim_utils/io/readline_completion.c
! so constructs like "2@q1", "orbit.x" and "-universe" complete as single words.
! Candidates are always full replacements for the token being completed.
!
! Known limitation: "-switch" words are dropped when locating the word position
! but the values of value-taking switches (eg "-write file", "-universe 2") are
! not, so completion after such a switch is off by one word. Fixing this needs
! per-switch arity information in tao_switch_sets.
!-

module tao_completion_mod

use tao_struct
use tao_command_names_mod
use tao_input_struct, only: tao_plot_page_input
use attribute_mod, only: attribute_info, ele_attribute_struct, attribute_type, attribute_index
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
!   context     -- character(*): 'LIST' = matches(:) holds the candidates,
!                    'FILE' = token is a file path (caller should do file name
!                    completion), 'NONE' = nothing to offer here.
!   matches(:)  -- character(100), allocatable: Candidate token replacements.
!                    At most max_matches$ are returned.
!   common_prefix -- character(*), optional: Longest common prefix of every
!                    candidate that matched the token, computed over all matches
!                    including any beyond the max_matches$ cap. This is what an
!                    interactive front end may safely insert.
!-

subroutine tao_complete (line, cursor, word_start, context, matches, common_prefix)

character(*), intent(in) :: line
integer, intent(in) :: cursor
integer, intent(out) :: word_start
character(*), intent(out) :: context
character(100), allocatable, intent(out) :: matches(:)
character(*), optional, intent(out) :: common_prefix

character(100) cand(max_matches$)
character(100) token, lcp, prepend
character(60) attrib_name
character(40) words(8), cmd_name, sub_name
character(1), parameter :: tab_char = achar(9)
character(1) quote
integer n_end, ix_semi, n_words, n_cand, n_lcp_hits, i, j, ix

!

context = 'NONE'
word_start = 1
n_cand = 0
n_lcp_hits = 0
lcp = ''
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

! word_start = start of the trailing (possibly empty) token.

word_start = ix_semi + 1
do i = n_end, ix_semi+1, -1
  if (line(i:i) == ' ' .or. line(i:i) == tab_char) then
    word_start = i + 1
    exit
  endif
enddo

token = ''
if (word_start <= n_end) token = line(word_start:n_end)

! Gather the complete words before the token, dropping "-switch" words.

n_words = 0
i = ix_semi + 1
do while (i < word_start)
  do while (i < word_start)
    if (line(i:i) /= ' ' .and. line(i:i) /= tab_char) exit
    i = i + 1
  enddo
  if (i >= word_start) exit
  j = i
  do while (j < word_start)
    if (line(j:j) == ' ' .or. line(j:j) == tab_char) exit
    j = j + 1
  enddo
  if (line(i:i) /= '-' .or. n_words == 0) then
    if (n_words < size(words)) then
      n_words = n_words + 1
      words(n_words) = line(i:j-1)
    endif
  endif
  i = j
enddo

!------------------------------------------
! First word: command names plus user defined aliases. Case sensitive.

if (n_words == 0) then
  if (token(1:1) /= '-') then
    call add_command_name_matches (.true.)
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
  if (at_set_value(attrib_name)) then
    call match_word (words(2), tao_set_target_names, ix, .true., matched_name = sub_name)
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
    call match_word (words(2), tao_set_target_names, ix, .true., matched_name = sub_name)
    select case (sub_name)

    ! These set targets work via a namelist read so a namelist write gives the component names.
    case ('global', 'beam_init', 'bmad_com', 'space_charge_com', 'geodesic_lm', &
          'opti_de_param', 'plot_page')
      call add_set_struct_matches (sub_name)

    ! Set via select case in tao_set_ptc_com_cmd. Keep in sync.
    case ('ptc_com')
      call add_prefix_matches ([character(24):: 'vertical_kick', 'cut_factor', &
            'max_fringe_order', 'old_integrator', 'exact_model', 'exact_misalign', &
            'use_orientation_patches', 'print_info_messages', 'pancake_symplectic', &
            'pancake_canonical'], .false.)

    ! Set via select case in tao_set_beam_cmd. Keep in sync (deprecated aliases omitted).
    case ('beam')
      call add_prefix_matches ([character(24):: 'beginning', 'comb_ds_save', &
            'always_reinit', 'track_start', 'track_end', 'beam_init_position_file', &
            'dump_file', 'dump_at', 'saved_at', 'add_saved_at', 'subtract_saved_at'], .false.)

    case ('element')
      call add_element_matches ()
    end select

  elseif (n_words == 3) then
    call match_word (words(2), tao_set_target_names, ix, .true., matched_name = sub_name)
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
    call add_command_name_matches (.false.)
  elseif (n_words == 2) then
    call match_word (words(2), [character(8):: 'pipe', 'python'], ix, matched_name = sub_name)
    if (ix > 0) call add_prefix_matches (tao_pipe_cmd_names, .false.)
  endif

case ('change')
  if (n_words == 1) then
    call add_prefix_matches ([character(20):: 'element', 'variable', 'tune', 'z_tune', &
                                              'particle_start'], .false.)
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

case ('call', 'read')
  if (n_words == 1) context = 'FILE'

end select

call finish()

!------------------------------------------
contains

subroutine finish ()
matches = cand(1:n_cand)
if (present(common_prefix)) common_prefix = trim(prepend) // lcp
end subroutine finish

!.................................

! The candidate-adding helpers below all mark the position as recognized by
! setting context = 'LIST'. Thus a recognized position with no matching
! candidate still reports 'LIST' (empty), distinct from 'NONE' (not a
! completion position). Callers therefore never set context themselves.

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

if (lt > 0) then
  n1 = name
  t1 = token
  if (.not. exact_case) then
    call str_upcase (n1, n1)
    call str_upcase (t1, t1)
  endif
  if (n1(1:lt) /= t1(1:lt)) return
endif

! Track the common prefix over every hit, including those beyond the candidate
! cap, so the interactive front end never inserts text a truncated list would
! not justify.

if (n_lcp_hits == 0) then
  lcp = name
else
  do im = 1, min(len_trim(lcp), len_trim(name))
    if (lcp(im:im) /= name(im:im)) exit
  enddo
  lcp = lcp(1:im-1)
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
! Top level command names, omitting internal commands that are not documented.

subroutine add_command_name_matches (exact_case)

logical exact_case
integer in

context = 'LIST'
do in = 1, size(tao_command_names)
  if (any(tao_command_names(in) == tao_hidden_command_names)) cycle
  call add_match_if_prefix (tao_command_names(in), exact_case)
enddo

end subroutine add_command_name_matches

!.................................
! Switch ("-flag") candidates from the context table in tao_command_names_mod.
! The context key is the command name, refined by the resolved subcommand for "show".

subroutine add_switch_matches ()

character(tao_switch_name_len) key
character(20) what_name
character(tao_switch_name_len), allocatable :: sw(:)
integer ik

context = 'LIST'

select case (cmd_name)
case ('show')
  if (n_words == 1) then
    key = 'show'
  else
    call match_word (words(2), tao_show_what_names, ik, matched_name = what_name)
    if (ik <= 0) return
    key = 'show ' // trim(what_name)
  endif
case ('place')
  ! The -no_buffer switch must come first.
  if (n_words > 1) return
  key = cmd_name
case default
  key = cmd_name
end select

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
! The "set" commands for these structs work by writing "<struct>%<component> = <value>"
! to a scratch file and doing a namelist read (see tao_set_mod). A namelist write of the
! same struct therefore enumerates exactly the component names that "set" accepts.

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
! Write the namelist for a "set" struct to a scratch file and return its unit
! (rewound, ready to read). ok is False if the file could not be opened. Scratch
! file I/O is safe inside the readline callback: only terminal output would
! corrupt the display.

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
! Split one namelist output line into its component name (downcased, '' if the
! line is not a component line) and its value text. A component line looks like
! "<struct>%<name>= value," (gfortran) or "<struct>%<name> = value," (ifort).
! Continuation lines of array values have no "<name>=" and give comp_name = ''.

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
! Current value text of one component from a "set" struct's namelist dump.

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
! Elements matching a Tao element selector ("Q01W", "quad::*", "1:10", "b>>q*",
! "2@q1", ...) via lat_ele_locator. Its error messages are suppressed: a selector
! that is still being typed is expected here and nothing may be printed from
! inside the readline callback.

subroutine locate_elements (selector, eles, n_loc)

type (ele_pointer_struct), allocatable :: eles(:)
type (out_io_output_direct_struct) out_state
type (tao_universe_struct), pointer :: u

character(*) selector
character(len(selector)) sel
integer n_loc, iuni, ia, ios
logical err

n_loc = 0
if (.not. allocated(s%u)) return

sel = selector
iuni = s%global%default_universe
ia = index(sel, '@')
if (ia > 1) then
  read (sel(1:ia-1), *, iostat = ios) iuni
  if (ios /= 0) return
  sel = sel(ia+1:)
endif
if (iuni < lbound(s%u, 1) .or. iuni > ubound(s%u, 1)) return
if (sel == '') return
u => s%u(iuni)

call output_direct (get = out_state)
call output_direct (print_and_capture = .false.)
call lat_ele_locator (sel, u%model%lat, eles, n_loc, err)
call output_direct (set = out_state)
if (err) n_loc = 0

end subroutine locate_elements

!.................................
! The beginning element and control lords (overlays, groups, girders, rampers)
! have attribute tables unlike any real element. When a selector matches several
! elements they are left out of the attribute intersection and value union.

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
integer n_loc, ie, in

context = 'LIST'
call locate_elements (ele_name, eles, n_loc)
if (n_loc == 0) return

ie = first_settable(eles, n_loc)
if (ie == 0) return

select case (attribute_type(upcase(attrib), eles(ie)%ele))
case (is_logical$)
  call add_prefix_matches ([character(8):: 'T', 'F'], .false.)
case (is_switch$)
  ! Union of the values valid for the matched elements.
  do ie = ie, min(n_loc, max_ele_scan$)
    if (.not. has_fixed_attributes(eles(ie)%ele)) cycle
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

! Attributes of the first matched element, then keep only those that every
! other matched element also has, since "set" applies to all of them.

ie0 = first_settable(eles, n_loc)
if (ie0 == 0) return
ie = ie0

! The extended range holds the non-value attributes: method switches such as
! tracking_method and space_charge_method, apertures, and logicals like field_master.

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

! Offer only what "set element" can actually set, which it resolves with
! pointer_to_attribute (see tao_set_mod). This drops lattice-file-only constructs
! such as superimpose or wall, and components not currently allocated.

do ia = 1, n_attr
  call pointer_to_attribute (eles(ie0)%ele, attrib_names(ia), .false., a_ptr, err, err_print_flag = .false.)
  if (err) cycle
  call add_match_if_prefix (downcase(attrib_names(ia)), .false.)
enddo

end subroutine add_attribute_matches

!.................................

! Element selector completion. Besides element names this understands the "n@"
! universe prefix and the "key::" element-type prefix of Tao's selector syntax:
! "quad::Q0" completes to the quadrupoles starting with Q0, and without a "::"
! the element types present in the lattice ("quadrupole::", ...) are offered too.

subroutine add_element_matches ()

type (tao_universe_struct), pointer :: u
type (branch_struct), pointer :: branch
character(100) sel
integer iu, ib, ie, ia, ic, ik, ix_key, ios
logical key_present(n_key$)

context = 'LIST'
if (.not. allocated(s%u)) return

sel = token
iu = s%global%default_universe
ia = index(sel, '@')
if (ia > 1) then
  read (sel(1:ia-1), *, iostat = ios) iu
  if (ios /= 0) return
  prepend = trim(prepend) // sel(1:ia)
  sel = sel(ia+1:)
endif
if (iu < lbound(s%u, 1) .or. iu > ubound(s%u, 1)) return
u => s%u(iu)

ix_key = 0
ic = index(sel, '::')
if (ic > 0) then
  if (ic == 1) return
  call match_word (sel(1:ic-1), key_name, ix_key)
  if (ix_key <= 0) return
  prepend = trim(prepend) // sel(1:ic+1)
  sel = sel(ic+2:)
endif
token = sel

! Element types present in the lattice go first so that a big lattice filling the
! candidate cap with element names cannot crowd them out.

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
    ! Element names are stored upcased so match case insensitively.
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
!   n_cand    -- integer(c_int): Number of candidates (0 = recognized position
!                  with nothing matching), or -1 meaning "not Tao's to complete:
!                  use readline's default file name completion". The engine's
!                  FILE and NONE contexts both map to -1 so the behavior before
!                  Tao completion existed is preserved there.
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
! Install tao_rl_complete_c as the readline tab completion callback. Idempotent.
!
! Called from tao_get_user_input immediately before a terminal line is read, so
! it runs only when Tao itself owns the prompt. Embedders such as PyTao (which
! drive Tao through tao_c_command) never reach that path and so never have the
! process-wide readline completion state changed under them.
!-

subroutine tao_register_completion ()

interface
  subroutine readline_set_completion_fn (fn) bind(c, name = 'readline_set_completion_fn')
    import :: c_funptr
    type(c_funptr), value :: fn
  end subroutine
end interface

logical, save :: registered = .false.

!

if (registered) return
registered = .true.
call readline_set_completion_fn (c_funloc(tao_rl_complete_c))

end subroutine tao_register_completion

end module tao_completion_mod
