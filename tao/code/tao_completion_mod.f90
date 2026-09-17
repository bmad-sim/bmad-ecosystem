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
!-

module tao_completion_mod

use tao_struct
use tao_command_names_mod
use tao_input_struct, only: tao_plot_page_input
use attribute_mod, only: attribute_info, ele_attribute_struct
use geodesic_lm, only: geodesic_lm_param_struct
use opti_de_mod, only: opti_de_param
use, intrinsic :: iso_c_binding

implicit none

integer, parameter, private :: max_matches$ = 200

contains

!------------------------------------------------------------------------------
!+
! Subroutine tao_complete (line, cursor, word_start, context, matches)
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
!-

subroutine tao_complete (line, cursor, word_start, context, matches)

character(*), intent(in) :: line
integer, intent(in) :: cursor
integer, intent(out) :: word_start
character(*), intent(out) :: context
character(100), allocatable, intent(out) :: matches(:)

character(100) cand(max_matches$)
character(100) token
character(40) words(8), cmd_name, sub_name
character(1), parameter :: tab_char = achar(9)
character(1) quote
integer n_end, ix_semi, n_words, n_cand, i, j, ix

!

context = 'NONE'
word_start = 1
n_cand = 0

if (.not. s%initialized) then
  allocate (matches(0))
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

! Switch completion is not supported.

if (token(1:1) == '-') then
  allocate (matches(0))
  return
endif

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
  call add_prefix_matches (tao_command_names, .true.)
  do i = 1, s%com%n_alias
    call add_match_if_prefix (s%com%alias(i)%name, .true.)
  enddo
  context = 'LIST'
  matches = cand(1:n_cand)
  return
endif

! Later words: context depends on the command.

call match_word (words(1), tao_command_names, ix, .true., matched_name = cmd_name)
if (ix <= 0) then
  allocate (matches(0))
  return
endif

select case (cmd_name)

case ('show')
  if (n_words == 1) then
    call add_prefix_matches (tao_show_what_names, .false.)
    ! Pseudo names remapped at the top of tao_show_this.
    call add_prefix_matches ([character(20):: 'plot_page', 'bmad_com', 'ptc_com', &
                                              'space_charge_com', 'floor_plan'], .false.)
    context = 'LIST'
  elseif (n_words == 2) then
    call match_word (words(2), tao_show_what_names, ix, matched_name = sub_name)
    if (sub_name == 'element') then
      call add_element_matches ()
      context = 'LIST'
    endif
  endif

case ('set')
  if (n_words == 1) then
    call add_prefix_matches (tao_set_target_names, .true.)
    context = 'LIST'
  else
    call match_word (words(2), tao_set_target_names, ix, .true., matched_name = sub_name)
    select case (sub_name)

    ! These set targets work via a namelist read so a namelist write gives the component names.
    case ('global', 'beam_init', 'bmad_com', 'space_charge_com', 'geodesic_lm', &
          'opti_de_param', 'plot_page')
      if (n_words == 2) then
        call add_set_struct_matches (sub_name)
        context = 'LIST'
      endif

    ! Set via select case in tao_set_ptc_com_cmd. Keep in sync.
    case ('ptc_com')
      if (n_words == 2) then
        call add_prefix_matches ([character(24):: 'vertical_kick', 'cut_factor', &
              'max_fringe_order', 'old_integrator', 'exact_model', 'exact_misalign', &
              'use_orientation_patches', 'print_info_messages', 'pancake_symplectic', &
              'pancake_canonical'], .false.)
        context = 'LIST'
      endif

    ! Set via select case in tao_set_beam_cmd. Keep in sync (deprecated aliases omitted).
    case ('beam')
      if (n_words == 2) then
        call add_prefix_matches ([character(24):: 'beginning', 'comb_ds_save', &
              'always_reinit', 'track_start', 'track_end', 'beam_init_position_file', &
              'dump_file', 'dump_at', 'saved_at', 'add_saved_at', 'subtract_saved_at'], .false.)
        context = 'LIST'
      endif

    case ('element')
      if (n_words == 2) then
        call add_element_matches ()
        context = 'LIST'
      elseif (n_words == 3) then
        call add_attribute_matches (words(3))
        context = 'LIST'
      endif
    end select
  endif

case ('pipe', 'python')
  if (n_words == 1) then
    call add_prefix_matches (tao_pipe_cmd_names, .false.)
    context = 'LIST'
  endif

case ('help')
  if (n_words == 1) then
    call add_prefix_matches (tao_command_names, .false.)
    context = 'LIST'
  elseif (n_words == 2) then
    call match_word (words(2), [character(8):: 'pipe', 'python'], ix, matched_name = sub_name)
    if (ix > 0) then
      call add_prefix_matches (tao_pipe_cmd_names, .false.)
      context = 'LIST'
    endif
  endif

case ('change')
  if (n_words == 1) then
    call add_prefix_matches ([character(20):: 'element', 'variable', 'tune', 'z_tune', &
                                              'particle_start'], .false.)
    context = 'LIST'
  elseif (n_words == 2) then
    if (words(2) /= '' .and. index('element', trim(words(2))) == 1) then
      call add_element_matches ()
      context = 'LIST'
    endif
  endif

case ('use', 'veto', 'restore')
  if (n_words == 1) then
    call add_prefix_matches ([character(8):: 'data', 'variable'], .true.)
    context = 'LIST'
  else
    call match_word (words(2), [character(8):: 'data', 'variable'], ix, .true., matched_name = sub_name)
    select case (sub_name)
    case ('data')
      call add_data_matches ()
      context = 'LIST'
    case ('variable')
      call add_var_matches ()
      context = 'LIST'
    end select
  endif

case ('call', 'read')
  if (n_words == 1) context = 'FILE'

end select

matches = cand(1:n_cand)

!------------------------------------------
contains

subroutine add_prefix_matches (names, exact_case)

character(*) names(:)
logical exact_case
integer in

do in = 1, size(names)
  call add_match_if_prefix (names(in), exact_case)
enddo

end subroutine add_prefix_matches

!.................................

subroutine add_match_if_prefix (name, exact_case)

character(*) name
logical exact_case
character(100) n1, t1
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

if (n_cand >= max_matches$) return
do im = 1, n_cand
  if (cand(im) == name) return
enddo
n_cand = n_cand + 1
cand(n_cand) = name

end subroutine add_match_if_prefix

!.................................
! The "set" commands for these structs work by writing "<struct>%<component> = <value>"
! to a scratch file and doing a namelist read (see tao_set_mod). A namelist write of the
! same struct therefore enumerates exactly the component names that "set" accepts.

subroutine add_set_struct_matches (set_word)

type (tao_global_struct) global
type (beam_init_struct) beam_init
type (bmad_common_struct) this_bmad_com
type (space_charge_common_struct) this_space_charge_com
type (geodesic_lm_param_struct) this_geodesic_lm
type (tao_plot_page_input) plot_page

character(*) set_word
character(300) nml_line
integer iu_nml, ios, ix1, ix2

namelist / nml_global / global
namelist / nml_beam_init / beam_init
namelist / nml_bmad_com / this_bmad_com
namelist / nml_space_charge_com / this_space_charge_com
namelist / nml_geodesic_lm / this_geodesic_lm
namelist / nml_opti_de_param / opti_de_param
namelist / nml_plot_page / plot_page

! Scratch file I/O is safe here: only terminal output would corrupt the readline display.

open (newunit = iu_nml, status = 'scratch', iostat = ios)
if (ios /= 0) return

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
do
  read (iu_nml, '(a)', iostat = ios) nml_line
  if (ios /= 0) exit
  nml_line = adjustl(nml_line)
  ix1 = index(nml_line, '%')
  ix2 = index(nml_line, '=')
  if (ix1 == 0 .or. ix2 <= ix1 + 1) cycle
  if (index(nml_line(1:ix2), ' ') /= 0 .or. index(nml_line(1:ix2), '"') /= 0) cycle
  call add_match_if_prefix (downcase(nml_line(ix1+1:ix2-1)), .false.)
enddo
close (iu_nml)

! Components handled as special cases before the namelist read.

select case (set_word)
case ('global');    call add_prefix_matches ([character(16):: 'phase_units', 'quiet'], .false.)
case ('plot_page'); call add_prefix_matches ([character(16):: 'title', 'subtitle'], .false.)
end select

end subroutine add_set_struct_matches

!.................................
! Attribute names for "set element <name> <attrib>", from bmad's attribute table
! for the first element matching <name> (wildcards allowed) in the default universe.

subroutine add_attribute_matches (ele_name)

type (tao_universe_struct), pointer :: u
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (ele_attribute_struct) attrib

character(*) ele_name
character(60) name_up
integer iuni, ib, ie, ia

if (ele_name == '') return
if (.not. allocated(s%u)) return
iuni = s%global%default_universe
if (iuni < lbound(s%u, 1) .or. iuni > ubound(s%u, 1)) return
u => s%u(iuni)

nullify(ele)
name_up = upcase(ele_name)
branch_loop: do ib = 0, ubound(u%model%lat%branch, 1)
  branch => u%model%lat%branch(ib)
  do ie = 1, branch%n_ele_max
    if (.not. match_wild(branch%ele(ie)%name, trim(name_up))) cycle
    ele => branch%ele(ie)
    exit branch_loop
  enddo
enddo branch_loop
if (.not. associated(ele)) return

do ia = 1, num_ele_attrib$
  attrib = attribute_info(ele, ia)
  if (attrib%name == null_name$) cycle
  if (attrib%state == private$) cycle
  call add_match_if_prefix (downcase(attrib%name), .false.)
enddo

end subroutine add_attribute_matches

!.................................

subroutine add_element_matches ()

type (tao_universe_struct), pointer :: u
type (branch_struct), pointer :: branch
integer iu, ib, ie

if (.not. allocated(s%u)) return
iu = s%global%default_universe
if (iu < lbound(s%u, 1) .or. iu > ubound(s%u, 1)) return
u => s%u(iu)

do ib = 0, ubound(u%model%lat%branch, 1)
  branch => u%model%lat%branch(ib)
  do ie = 1, branch%n_ele_max
    ! Element names are stored upcased so match case insensitively.
    call add_match_if_prefix (branch%ele(ie)%name, .false.)
    if (n_cand >= max_matches$) return
  enddo
enddo

end subroutine add_element_matches

!.................................

subroutine add_data_matches ()

type (tao_universe_struct), pointer :: u
type (tao_d2_data_struct), pointer :: d2
character(100) name
integer iu, id, id1

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
!   buf_c     -- type(c_ptr): Buffer filled with newline separated candidates,
!                  null terminated.
!   n_cand    -- integer(c_int): Number of candidates, or -1 for "file path:
!                  use readline's default file name completion".
!-

function tao_rl_complete_c (line_c, point, istart, iend, buf_c, buf_size) bind(c) result (n_cand)

type(c_ptr), value :: line_c, buf_c
integer(c_int), value :: point, istart, iend, buf_size
integer(c_int) :: n_cand

character(kind=c_char), pointer :: line_p(:), buf_p(:)
character(4000) line_f
character(8) context
character(100), allocatable :: matches(:)
integer word_start, i, j, k, lt

!

n_cand = 0
if (.not. s%initialized) return
if (.not. c_associated(line_c) .or. .not. c_associated(buf_c)) return

call c_f_pointer (line_c, line_p, [len(line_f)])
line_f = ''
do i = 1, len(line_f)
  if (line_p(i) == c_null_char) exit
  line_f(i:i) = line_p(i)
enddo

call tao_complete (line_f, min(point, len(line_f)) + 1, word_start, context, matches)

if (context == 'FILE') then
  n_cand = -1
  return
endif

! The engine and readline must agree on the token span. If they do not
! (for example a ";" inside the token), do not offer anything.

if (word_start - 1 /= istart) return

call c_f_pointer (buf_c, buf_p, [buf_size])
k = 0
do i = 1, size(matches)
  lt = len_trim(matches(i))
  if (k + lt + 2 > buf_size) exit
  do j = 1, lt
    k = k + 1
    buf_p(k) = matches(i)(j:j)
  enddo
  k = k + 1
  buf_p(k) = c_new_line
  n_cand = n_cand + 1
enddo
k = k + 1
buf_p(k) = c_null_char

end function tao_rl_complete_c

!------------------------------------------------------------------------------
!+
! Subroutine tao_register_completion ()
!
! Install tao_rl_complete_c as the readline tab completion callback.
! Called once at startup (see tao_top_level). Only affects interactive
! terminal input; command files and the pipe interface never enter readline.
!-

subroutine tao_register_completion ()

interface
  subroutine readline_set_completion_fn (fn) bind(c, name = 'readline_set_completion_fn')
    import :: c_funptr
    type(c_funptr), value :: fn
  end subroutine
end interface

!

call readline_set_completion_fn (c_funloc(tao_rl_complete_c))

end subroutine tao_register_completion

end module tao_completion_mod
