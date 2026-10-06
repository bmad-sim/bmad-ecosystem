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
use attribute_mod, only: attribute_info, ele_attribute_struct
use geodesic_lm, only: geodesic_lm_param_struct
use opti_de_mod, only: opti_de_param
use, intrinsic :: iso_c_binding

implicit none

integer, parameter, private :: max_matches$ = 200

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
character(100) token, lcp
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
  if (n_words == 1) then
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
if (present(common_prefix)) common_prefix = lcp
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
do im = 1, n_cand
  if (cand(im) == name) return
enddo
n_cand = n_cand + 1
cand(n_cand) = name

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

type (tao_global_struct) global
type (beam_init_struct) beam_init
type (bmad_common_struct) this_bmad_com
type (space_charge_common_struct) this_space_charge_com
type (geodesic_lm_param_struct) this_geodesic_lm
type (tao_plot_page_input) plot_page

character(*) set_word
character(300) nml_line
character(*), parameter :: name_chars = &
          'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_%'
integer iu_nml, ios, ix1, ix2

namelist / nml_global / global
namelist / nml_beam_init / beam_init
namelist / nml_bmad_com / this_bmad_com
namelist / nml_space_charge_com / this_space_charge_com
namelist / nml_geodesic_lm / this_geodesic_lm
namelist / nml_opti_de_param / opti_de_param
namelist / nml_plot_page / plot_page

context = 'LIST'

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
  ! A component line looks like "<struct>%<name>= value," (gfortran) or
  ! "<struct>%<name> = value," (ifort). Take the name as the run of name
  ! characters after the "%" and require an "=" after it, which also skips
  ! continuation lines of array values.
  nml_line = adjustl(nml_line)
  ix1 = index(nml_line, '%')
  if (ix1 == 0) cycle
  ix2 = ix1 + 1
  do while (ix2 <= len(nml_line))
    if (verify(nml_line(ix2:ix2), name_chars) /= 0) exit
    ix2 = ix2 + 1
  enddo
  if (ix2 == ix1 + 1) cycle
  if (index(nml_line(ix2:), '=') == 0) cycle
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

context = 'LIST'
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

context = 'LIST'
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
