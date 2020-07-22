!+
! Module ptc_interface_mod
!
! Module of basic PTC interface routines.
! Also see: ptc_layout_mod
!-

module ptc_interface_mod

use bmad_interface

interface assignment (=)
  module procedure real_8_equal_bmad_taylor
  module procedure ptc_taylor_equal_bmad_taylor
  module procedure bmad_taylor_equal_real_8
  module procedure universal_equal_universal
  module procedure bmad_taylor_equal_ptc_taylor
  module procedure bmad_taylors_equal_ptc_taylors
  module procedure complex_taylor_equal_c_taylor
  module procedure complex_taylors_equal_c_taylors
end interface

interface operator (+)
  module procedure taylor_plus_taylor
end interface

interface operator (-)
  module procedure taylor_minus_taylor
end interface

type ptc_parameter_struct
  logical exact_model
  logical exact_misalign
  real(dp) phase0
end type

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Function taylor_plus_taylor (taylor1, taylor2) result (taylor3)
!
! Routine to add two taylor maps.
!
! Input:
!   taylor1 -- Taylor_struct:
!   taylor2 -- Taylor_struct:
!
! Output:
!   taylor3 -- Taylor_struct:
!-

function taylor_plus_taylor (taylor1, taylor2) result (taylor3)

use polymorphic_taylor, only: alloc, kill, operator(+), real_8

implicit none

type (taylor_struct), intent(in) :: taylor1(:), taylor2(:)
type (taylor_struct) taylor3(size(taylor1))
type (real_8) y1(size(taylor1)), y2(size(taylor1)), y3(size(taylor1))

integer i

! set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

!

call alloc(y1)
call alloc(y2)
call alloc(y3)

y1 = taylor1
y2 = taylor2

do i = 1, size(taylor1)
  y3(i) = y1(i) + y2(i)
enddo

taylor3 = y3

call kill (y1)
call kill (y2)
call kill (y3)

end function taylor_plus_taylor

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Function taylor_minus_taylor (taylor1, taylor2) result (taylor3)
!
! Routine to add two taylor maps.
!
! Input:
!   taylor1 -- Taylor_struct:
!   taylor2 -- Taylor_struct:
!
! Output:
!   taylor3 -- Taylor_struct:
!-

function taylor_minus_taylor (taylor1, taylor2) result (taylor3)

use polymorphic_taylor, only: kill, operator(-), real_8, alloc

implicit none

type (taylor_struct), intent(in) :: taylor1(:), taylor2(:)
type (taylor_struct) taylor3(size(taylor1))
type (real_8) y1(size(taylor1)), y2(size(taylor1)), y3(size(taylor1))

integer i

! set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

!

call alloc(y1)
call alloc(y2)
call alloc(y3)

y1 = taylor1
y2 = taylor2

do i = 1, size(taylor1)
  y3(i) = y1(i) - y2(i)
enddo

taylor3 = y3

call kill (y1)
call kill (y2)
call kill (y3)

end function taylor_minus_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Function map_coef(y, i, j, k, l)
!
! Function to return the coefficient of the map y(:) up to 3rd order.
! Example: 
!   In 4-dimensional space with
!      X = (x_1, x_2, x_3, x_4)
!   And if y(:) is the map of X(in) to X(out) then
!      map_coef(y, 2, 1, 4)
!   Gives the coefficient for
!     x_2(out) = ... + coef * x_1(in) * x_4(in) + ...
!
! Notice that map_coef(y, i, j) just gives the linear (matrix) part of the map.
!
! Input:
!   y(:)  -- Real_8: Taylor Map.
!   i     -- Integer: output index.
!   j     -- Integer, optional: 1st input index, needed for 1st order and above.
!   k     -- Integer, optional: 2nd input index, needed for 2nd order and above.
!   l     -- Integer, optional: 3rd input index, needed for 3rd order.
!
! Output:
!   map_coef -- Real(8): Coefficient.
!-

function map_coef (y, i, j, k, l) result (m_coef)

use polymorphic_taylor, only: operator (.sub.), operator(*), real_8

implicit none

type (real_8) y(:)


real(8) m_coef

integer i
integer, optional :: j, k, l
integer arr(40), n_max, sgn, ii

character str*40
character, parameter :: str1(4) = ['1', '2', '3', '4' ]

logical use_bmad

!

arr = 0
sgn = 1
str = '0000000000000000000000000000000000000000'
n_max = 1

call map_index(j)
call map_index(k)
call map_index(l)

call map_index(i, ii)
m_coef = sgn * y(ii)%t .sub. str(1:n_max)

!--------------------------------------------------------------
contains

subroutine map_index (iz, i_in)

integer, optional :: iz, i_in
integer n0

!

if (.not. present(iz)) return

  n0 = iz

  if (use_bmad) then
    if (iz == 5) then
      n0 = 6
      sgn = -sgn
    elseif (iz == 6) then
      n0 = 5
    endif
  endif

! i_in is present only with the input index.

  if (present(i_in)) then
    i_in = n0
    return
  endif

  arr(n0) = arr(n0) + 1
  str(n0:n0) = str1(arr(n0))

  n_max = max(n_max, n0)

end subroutine map_index

end function map_coef

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine type_map1 (y, type0, n_dim)
!
! Subroutine to type the transfer map up to first order.
!
! Input:
!   y(:)  -- Real_8: 
!   type0 -- Logical: Type zeroth order map
!   n_dim -- Integer: Number of dimensions to type: 4 or 6
!-

subroutine type_map1 (y, type0, n_dim)

use definition, only: real_8

implicit none

type (real_8), intent(in) :: y(:)

integer, intent(in) :: n_dim
integer :: i, j

logical, intent(in) :: type0

character(80) line
character(16), parameter :: r_name = 'type_map1'
!

if (type0) then
  call out_io (s_blank$, r_name, '0th Order Map:')
  do i = 1, n_dim
    write (line, '(6f11.5)') map_coef(y(:), i) 
    call out_io (s_blank$, r_name, line)
  enddo
  call out_io (s_blank$, r_name, '')
endif

call out_io (s_blank$, r_name, '1st Order Map:')
do i = 1, n_dim
  write (line, '(6f11.5)') (map_coef(y(:), i, j), j = 1, n_dim)
  call out_io (s_blank$, r_name, line)
enddo

end subroutine type_map1

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine type_ptc_internal_state (intern_state, line, n_lines)
!
! Routine to print information on a PTC internal state.
!
! Input:
!   intern_state -- Internal_state, optional: PTC state. If not present then the PTC 
!                     state DEFAULT is used.
!
! Output:
!   lines(:)  -- character(100), optional, allocatable: Character array to hold the output.
!   n_lines   -- integer, optional: Number of lines used in lines(:)
!-

subroutine type_ptc_internal_state (intern_state, lines, n_lines)

use s_status, only: internal_state, default

implicit none

type (internal_state), optional, target :: intern_state
type (internal_state), pointer :: state_ptr

integer, optional :: n_lines
integer i, nl

character(*), allocatable, optional :: lines(:)
character(100), allocatable :: li(:)

!

state_ptr => default
if (present(intern_state)) state_ptr => intern_state

allocate (li(20))

nl = 0
nl=nl+1; li(nl) = 'Internal_state:'
nl=nl+1; write (li(nl), '(a, t16, i0)') '  %totalpath:    ', state_ptr%totalpath
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %time:         ', state_ptr%time
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %radiation:    ', state_ptr%radiation
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %nocavity:     ', state_ptr%nocavity
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %fringe:       ', state_ptr%fringe
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %stochastic:   ', state_ptr%stochastic
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %envelope:     ', state_ptr%envelope
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %para_in:      ', state_ptr%para_in
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %only_4d:      ', state_ptr%only_4d
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %delta:        ', state_ptr%delta
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %spin:         ', state_ptr%spin
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %modulation:   ', state_ptr%modulation
nl=nl+1; write (li(nl), '(a, t16, l1)') '  %only_2d:      ', state_ptr%only_2d

call type_end_stuff(li, nl, lines, n_lines)

end subroutine type_ptc_internal_state

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine type_ptc_fibre (ptc_fibre, print_coords, lines, n_lines)
!
! Routine to put information on a PTC fibre element into a string array.
! If "lines" is not present, the information will be printed to the screen.
!
! Input:
!   ptc_fibre    -- fibre, pointer: Fibre to type info of.
!   print_coords -- logical, optional: If True then print coordinate and patch information. 
!                     Default is True.
!
! Output:
!   lines(:)  -- character(120), optional, allocatable: Character array to hold the output.
!   n_lines   -- integer, optional: Number of lines used in lines(:)
!-

subroutine type_ptc_fibre (ptc_fibre, print_coords, lines, n_lines)

use definition, only: patch, element, elementp, magnet_chart, lp, chart, magnet_frame

implicit none

type (fibre), pointer :: ptc_fibre
type (patch), pointer :: ptch
type (element), pointer :: mag
type (elementp), pointer :: magp
type (magnet_chart), pointer :: p
type (chart), pointer :: chrt
type (magnet_frame), pointer :: frame

integer, optional :: n_lines
integer i, nl

character(*), allocatable, optional :: lines(:)
character(160), allocatable :: li(:)
character(160) str
character(16), parameter :: r_name = 'type_ptc_fibre'

logical, optional :: print_coords
logical printit, all_zero

!

nl = 0
call re_allocate(li, 100, .false.)

!

if (.not. associated(ptc_fibre)) then
  nl=nl+1; li(nl) = 'Fibre pointer not associated!'
  call type_end_stuff(li, nl, lines, n_lines)
  return
endif

! Fibre

nl=nl+1; li(nl) = 'Fibre:'
nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%pos    [Index in layout]:    ', int_of(ptc_fibre%pos, 'i0')
nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%loc    [Index in universe]:  ', int_of(ptc_fibre%loc, 'i0')
nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%dir    [Direction]:          ', int_of(ptc_fibre%dir, 'i0')
nl=nl+1; write (li(nl), '(2x, a, t39, a)') '%beta0  [Beta velocity]:      ', real_of(ptc_fibre%beta0, 'es13.5')
nl=nl+1; write (li(nl), '(2x, a, t39, a)') '%mass   [Mass]:               ', real_of(ptc_fibre%mass, 'es13.5', 1e9_rp)
nl=nl+1; write (li(nl), '(2x, a, t39, a)') '%charge [Charge]:             ', real_of(ptc_fibre%charge, 'es13.5')

!

if (logic_option(.true., print_coords)) then

  ! ptc_fibre%chart

  nullify(frame)
  chrt => ptc_fibre%chart

  nl=nl+1; li(nl) = ''
  if (associated(chrt)) then
    nl=nl+1; li(nl) = 'fibre%chart (Type: chart):'
    nl=nl+1; li(nl) = '  [Fibre to Magnet patch at Entrance End:]'
    call write_reals ('%d_in   [displacement]  ', chrt%d_in,   '3f11.7')
    call write_reals ('%ang_in [angle]         ', chrt%ang_in, '3f11.7')
    nl=nl+1; li(nl) = '  [Magnet to Fibre patch at Exit End:]'
    call write_reals ('%d_out   [displacement]  ', chrt%d_out,   '3f11.7')
    call write_reals ('%ang_out [angle]         ', chrt%ang_out, '3f11.7')
    frame => chrt%f
  else
    nl=nl+1; li(nl) = 'fibre%chart (Type: chart): Not Associated.'
  endif

  ! ptc_fibre%chart%f

  nl=nl+1; li(nl) = ''
  if (associated(frame)) then
    nl=nl+1; li(nl) = 'fibre%chart%f (Type: magnet_frame):'
    nl=nl+1; li(nl) = '  [Global orientation of the Entrance end of the Fibre]:'
    call write_reals  ('%a   [position]:   ', frame%a,   '3f11.7')
    call write_real2d ('%ent [Orientation]:', frame%ent, '3f11.7')
    nl=nl+1; li(nl) = '  [Global orientation of the Center of the Fibre]:'
    call write_reals  ('%o   [position]:   ', frame%o,   '3f11.7')
    call write_real2d ('%mid [Orientation]:', frame%mid, '3f11.7')
    nl=nl+1; li(nl) = '  [Global orientation of the Exit end of the Fibre]:'
    call write_reals  ('%b   [position]:   ', frame%b,   '3f11.7')
    call write_real2d ('%exi [Orientation]:', frame%exi, '3f11.7')

  else
    nl=nl+1; li(nl) = 'fibre%chart (Type: chart): Not Associated'
  endif

endif

! ptc_fibre%mag

nullify(p)

mag => ptc_fibre%mag
nl=nl+1; li(nl) = ''
if (associated(mag)) then
  nl=nl+1; li(nl) = 'fibre%mag (Type: element):'
  nl=nl+1; write (li(nl), '(2x, a, t40, a)')  '%name   [Name]:     ', name_of(mag%name)
  nl=nl+1; write (li(nl), '(2x, a, t40, a)')  '%kind   [Kind]:     ', kind_name(mag%kind)
  call is_associated (associated(mag%d0),     '%d0     [Drift]:')
  call is_associated (associated(mag%k2),     '%k2     [Integrator]:')
  call is_associated (associated(mag%k3),     '%k3     [Thin Kick]:')
  call is_associated (associated(mag%c4),     '%c4     [Cavity]:')
  call is_associated (associated(mag%s5),     '%s5     [Solenoid]:')
  call is_associated (associated(mag%t6),     '%t6     [Integrator thick,slow]:')
  call is_associated (associated(mag%t7),     '%t7     [Integrator thick,fast]:')
  call is_associated (associated(mag%s8),     '%s8     [Normal SMI]:')
  call is_associated (associated(mag%s9),     '%s9     [Skew SMI]:')
  call is_associated (associated(mag%tp10),   '%tp10   [Sector Teapot]:')
  call is_associated (associated(mag%mon14),  '%mon14  [Monitor/Instrument]:')
  call is_associated (associated(mag%sep15),  '%sep15  [Monitor/Instrument]:')
  call is_associated (associated(mag%k16),    '%k16    [Exact Straight Integrator]:')
  call is_associated (associated(mag%enge17), '%enge17 [Solenoid Sixtrack style]:')
  call is_associated (associated(mag%rcol18), '%rcol18 [Rcollimator]:')
  call is_associated (associated(mag%ecol19), '%ecol19 [Ecollimator]:')
  call is_associated (associated(mag%cav21),  '%cav21  [Cavity. Traveling Wave]:')
  call is_associated (associated(mag%wi),     '%wi     [Wiggler]:')
  call is_associated (associated(mag%pa),     '%pa     [General B]:')
  call is_associated (associated(mag%he22),   '%he22   [Helical Dipole]:')

  call write_real ('%l [Integration Len] ',  mag%l,       'es13.5', writeit = .true.)
  call write_real ('%volt                ',  mag%volt,    'es13.5')
  call write_real ('%freq                ',  mag%freq,    'es13.5')
  call write_real ('%phas                ',  mag%phas,    'es13.5')
  call write_real ('%lag                 ',  mag%lag,     'es13.5')
  call write_real ('%delta_e             ',  mag%delta_e, 'es13.5')
  call write_reals('%fint                ',  mag%fint,    '2es13.5')
  call write_reals('%hgap                ',  mag%hgap,    '2es13.5')
  call write_real ('%h1                  ',  mag%h1,      'es13.5')
  call write_real ('%h2                  ',  mag%h2,      'es13.5')
  call write_real ('%b_sol               ',  mag%b_sol,   'es13.5')
  call write_real ('%va [sad f1]         ',  mag%va,      'es13.5')
  call write_real ('%vs [sad f2]         ',  mag%vs,      'es13.5')

  if (associated(mag%an)) then
    all_zero = .true.
    do i = lbound(mag%an, 1), ubound(mag%an, 1)
      if (mag%an(i) == 0 .and. mag%bn(i) == 0) cycle
      nl=nl+1; write (li(nl), '(2x, 2(a, i0), a, t38, 2es13.5)') &
                          '%an(', i, '), %bn(', i, '):           ', mag%an(i), mag%bn(i)
      all_zero = .false.
    enddo
    if (all_zero) then
      nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%an(:), %bn(:):', 'All Zero'
    endif
  else
    nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%an(:), %bn(:):', 'Not Associated'
  endif

  nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%thin:        ', logical_of(mag%thin)

  p =>  ptc_fibre%mag%p

else
  nl=nl+1; li(nl) = 'fibre%mag (Type: element): Not Associated.'
endif

! ptc_fibre%mag%p

nullify(frame)

nl=nl+1; li(nl) = ''
if (associated(p)) then
  nl=nl+1; li(nl) = 'fibre%mag%p (Type: magnet_chart):'
  call write_real  ('%p0c   [P0C]     ',  p%p0c,   'es13.5', 1e9_rp)
  call write_real  ('%ld    [L]:      ',  p%ld,    'es13.5')
  call write_real  ('%lc    [L_chord]:',  p%lc,    'es13.5')
  call write_real  ('%tiltd [Tilt]:   ',  p%tiltd, 'es13.5')
  call write_real  ('%b0    [Rho]:    ',  p%b0,    'es13.5')
  call write_reals ('%edge  [E1, E2]: ',  p%edge,  '10es13.5')
  nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%method [Integration order]: ', int_of(p%method, 'i0')
  nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%nst    [Integration steps]: ', int_of(p%nst, 'i0')
  nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%exact  [Exact integration]: ', logical_of(p%exact)
  frame => p%f

else
  nl=nl+1; li(nl) = 'fibre%mag%p (Type: magnet_chart): Not Associated.'
endif

if (logic_option(.true., print_coords)) then

  ! ptc_fibre%mag%p%f

  nl=nl+1; li(nl) = ''
  if (associated(frame)) then
    nl=nl+1; li(nl) = 'fibre%mag%p%f (Type: magnet_frame):'
    nl=nl+1; li(nl) = '  [Global orientation of the Entrance end of the Magnet]:'
    call write_reals  ('%a   [position]:   ', frame%a,   '3f11.7')
    call write_real2d ('%ent [Orientation]:', frame%ent, '3f11.7')
    nl=nl+1; li(nl) = '  [Global orientation of the Exit end of the Magnet]:'
    call write_reals  ('%b   [position]:   ', frame%b,   '3f11.7')
    call write_real2d ('%exi [Orientation]:', frame%exi, '3f11.7')

  else
    nl=nl+1; li(nl) = 'fibre%chart (Type: chart): Not Associated'
  endif

endif

! ptc_fibre%magp

magp => ptc_fibre%magp
nl=nl+1; li(nl) = ''
if (associated(magp)) then
  nl=nl+1; li(nl) = 'fibre%magp (Type: elementp):'
  nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%knob:        ', logical_of(magp%knob)

else
  nl=nl+1; li(nl) = 'fibre%magp (Type: elementp): Not Associated.'
endif

!

if (logic_option(.true., print_coords)) then

  ptch => ptc_fibre%patch
  nl=nl+1; li(nl) = ''
  if (associated(ptch)) then
    nl=nl+1; li(nl) = 'fibre%patch (Type: patch):'
    nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%patch  [Patch at]:            ', patch_name(ptch%patch)

    printit = .false.
    if (integer2_value(0, ptch%patch) == 1 .or. integer2_value(0, ptch%patch) == 3) printit = .true.

    if (printit) then
      nl=nl+1
      call write_reals ('%a_d    [Entrance translation]:', ptch%a_d,   '10f13.9', writeit = .true.)
      call write_reals ('%a_ang  [Entrance angle]:      ', ptch%a_ang, '10f13.9', writeit = .true.)
    endif

    printit = .false.
    if (integer2_value(0, ptch%patch) == 2 .or. integer2_value(0, ptch%patch) == 3) printit = .true.

    if (printit) then
      call write_reals ('%b_d    [Exit translation]:', ptch%b_d,   '10f13.9', writeit = .true.)
      call write_reals ('%b_ang  [Exit angle]:      ', ptch%b_ang, '10f13.9', writeit = .true.)
    endif

    nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%energy [Energy Patch at]:     ', patch_name(ptch%energy)
    nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%time   [Time Patch at]:       ', patch_name(ptch%time)

    printit = .false.
    if (integer2_value(0, ptch%time) == 1 .or. integer2_value(0, ptch%time) == 3) printit = .true.

    if (printit) then
      call write_real ('%a_t    [Patch Entrance dTime]:', ptch%a_t,   'es13.5', writeit = .true.)
    endif

    printit = .false.
    if (integer2_value(0, ptch%time) == 2 .or. integer2_value(0, ptch%time) == 3) printit = .true.

    if (printit) then
      call write_real ('%b_t    [Patch Exit dTime]:', ptch%b_t,   'es13.5', writeit = .true.)
    endif



  else
    nl=nl+1; li(nl) = 'fibre%patch (Type: patch): Not Associated.'
  endif

endif

!

call type_end_stuff(li, nl, lines, n_lines)

!-----------------------------------------------------------------------------
contains

function integer2_value (default, int_ptr) result (int_val)

integer default, int_val
integer(2), pointer :: int_ptr

!

if (associated(int_ptr)) then
  int_val = int_ptr
else
  int_val = default
endif

end function integer2_value

!-----------------------------------------------------------------------------
! contains

function logical_of (logical_in) result (str)

logical(lp), pointer :: logical_in
character(20) str

!

if (associated(logical_in)) then
  write (str, '(l1)') logical_in
else
  str = 'Not Associated'
endif

end function logical_of

!-----------------------------------------------------------------------------
! contains

subroutine is_associated (is_assoc, str)

logical is_assoc
character(*) str

!

if (.not. is_assoc) return
nl= nl+1; write (li(nl), '(2x, a, t40, a)') str, 'Associated'

end subroutine is_associated

!-----------------------------------------------------------------------------
! contains

function patch_name (patch) result (patch_str)

integer(2), pointer :: patch
character(40) patch_str

!

if (.not. associated(patch)) then
  patch_str = 'Not Associated.'
  return
endif

! Internal patches are used with energy offsets where the information about the before and
! after reference energies is contained within the fibre and not contained in the 
! previous or next fibres.

select case (patch)
case (0); patch_str = 'No patches (0)'
case (1); patch_str = 'Entrance patch (1)'
case (2); patch_str = 'Exit patch (2)'
case (3); patch_str = 'Patch both ends (3)'
case (4); patch_str = 'Internal entrance patch (4)'
case (5); patch_str = 'Internal exit patch (5)'
case default; write (patch_str, '(a, i0, a)') 'UNKNOWN! [', patch, ']'
end select

end function patch_name

!-----------------------------------------------------------------------------
! contains

function name_of (name_in) result (name_out)

character(*), pointer :: name_in
character(len(name_in)) :: name_out

if (associated(name_in)) then
  name_out = name_in
else
  name_out = 'Not Associated'
endif

end function name_of

!-----------------------------------------------------------------------------
! contains

function real_of (real_in, fmt, mult) result (str)

real(dp), pointer :: real_in
real(rp), optional :: mult

character(*) fmt
character(20) fmt2
character(100) str

if (associated(real_in)) then
  fmt2 = '(' // trim(fmt) // ')'
  if (present(mult)) then
    write (str, fmt2) real_in * mult
  else
    write (str, fmt2) real_in
  endif
  str = adjustl(str)
  if (str(1:1) /= '-') str = ' ' // str
else
  str = 'Not Associated'
endif

end function real_of

!-----------------------------------------------------------------------------
! contains

subroutine write_real (pre_str, real_in, real_fmt, mult, writeit)

real(dp), pointer :: real_in
real(rp), optional :: mult

logical, optional :: writeit
character(*) real_fmt, pre_str
character(20) fmt2
character(100) str

! Default is to no write if real_in is not present

if (.not. associated(real_in)) then
  if (.not. logic_option(.false., writeit)) return
  str = ' Not Associated'
else
  fmt2 = '(' // trim(real_fmt) // ')'
  if (present(mult)) then
    write (str, fmt2) real_in * mult
  else
    write (str, fmt2) real_in
  endif
  str = adjustl(str)
  if (str(1:1) /= '-') str = ' ' // str
endif

nl = nl + 1
write (li(nl), '(2x, a, t39, a)') pre_str, str

end subroutine write_real

!-----------------------------------------------------------------------------
! contains

subroutine write_reals (pre_str, real_in, real_fmt, mult, writeit)

real(dp), pointer :: real_in(:)
real(rp), optional :: mult

logical, optional :: writeit
character(*) real_fmt, pre_str
character(20) fmt2
character(100) str

! Default is to no write if real_in is not present

if (.not. associated(real_in)) then
  if (.not. logic_option(.false., writeit)) return
  str = ' Not Associated'
else
  fmt2 = '(' // trim(real_fmt) // ')'
  if (present(mult)) then
    write (str, fmt2) real_in * mult
  else
    write (str, fmt2) real_in
  endif
  str = adjustl(str)
  if (str(1:1) /= '-') str = ' ' // str
endif

nl = nl + 1
write (li(nl), '(2x, a, t39, a)') pre_str, str

end subroutine write_reals

!-----------------------------------------------------------------------------
! contains

subroutine write_real2d (pre_str, real2d, real_fmt, mult, writeit)

real(dp), pointer :: real2d(:, :)
real(dp), pointer :: real1d(:)
real(rp), optional :: mult

integer i, n
logical, optional :: writeit
character(*) real_fmt, pre_str
character(100) str

! Default is to no write if real_in is not present

if (.not. associated(real2d)) then
  nullify(real1d)
  call write_reals(pre_str, real1d, real_fmt, mult, writeit)
else
  real1d  => real2d(1,:)
  call write_reals(pre_str, real1d, real_fmt, mult, writeit)

  str = ''
  do i = 2, size(real2d, 1)
    real1d  => real2d(i,:)
    call write_reals(str(1:len(pre_str)), real1d, real_fmt, mult, writeit)
  enddo

endif

end subroutine write_real2d

!-----------------------------------------------------------------------------
! contains

function int_of (int_in, fmt) result (str)

integer, pointer :: int_in
character(*) fmt
character(20) fmt2
character(100) str

if (associated(int_in)) then
  fmt2 = '(' // trim(fmt) // ')'
  write (str, fmt2) int_in
else
  str = 'Not Associated'
endif

end function int_of

end subroutine type_ptc_fibre

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!
!-

subroutine type_end_stuff(li, nl, lines, n_lines)

integer, optional :: n_lines
integer i, nl

character(*), allocatable, optional :: lines(:)
character(*), allocatable :: li(:)

if (present(lines)) then
  call re_allocate(lines, nl, .false.)
  n_lines = nl
  lines(1:nl) = li(1:nl)
else
  do i = 1, nl
    print *, trim(li(i))
  enddo
endif

end subroutine type_end_stuff

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Function kind_name (this_kind) result (kind_str)
!
! function to return the name of a PTC kind.
!
! Input:
!   this_kind -- Integer, pointer: PTC kind 
!
! Output:
!   kind_str -- Character(60): String representation
!-

function kind_name (this_kind) result (kind_str)

use s_status, only: kind0, kind1, kind2, kind3, kind4, kind5, kind6, kind7, kind8, kind9, &
                    kind10, kind11, kind12, kind13, kind14, kind15, kind16, kind17, kind18, &
                    kind19, kind20, kind21, kind22, kind23, kindwiggler, kindpa

implicit none

integer, pointer :: this_kind
character(60) kind_str

!

if (.not. associated(this_kind)) then
  kind_str = 'Not Associated.'
  return
endif

! Generally the integer value of kindX is X+30

select case (this_kind)
case (KIND0);        write (kind_str, '(i0, a)') kind0, ': KIND0 [Marker, Patch]'
case (KIND1);        write (kind_str, '(i0, a)') kind1, ': KIND1 [Drift]'
case (KIND2);        write (kind_str, '(i0, a)') kind2, ': KIND2 [Drift-Kick-Drift EXACT_MODEL = F]'
case (KIND3);        write (kind_str, '(i0, a)') kind3, ': KIND3 [Thin Element, L == 0]'
case (KIND4);        write (kind_str, '(i0, a)') kind4, ': KIND4 [RFcavity]'
case (KIND5);        write (kind_str, '(i0, a)') kind5, ': KIND5 [Solenoid]'
case (KIND6);        write (kind_str, '(i0, a)') kind6, ': KIND6 [Kick-SixTrack-Kick]'
case (KIND7);        write (kind_str, '(i0, a)') kind7, ': KIND7 [Matrix-Kick-Matrix]'
case (KIND8);        write (kind_str, '(i0, a)') kind8, ': KIND8 [Normal SMI]'
case (KIND9);        write (kind_str, '(i0, a)') kind9, ': KIND9 [Skew SMI]'
case (KIND10);       write (kind_str, '(i0, a)') kind10, ': KIND10 [Sector Bend, Exact_Model]'
case (KIND11);       write (kind_str, '(i0, a)') kind11, ': KIND11 [Monitor]'
case (KIND12);       write (kind_str, '(i0, a)') kind12, ': KIND12 [HMonitor]'
case (KIND13);       write (kind_str, '(i0, a)') kind13, ': KIND13 [VMonitor]'
case (KIND14);       write (kind_str, '(i0, a)') kind14, ': KIND14 [Instrument]'
case (KIND15);       write (kind_str, '(i0, a)') kind15, ': KIND15 [ElSeparator]'
case (KIND16);       write (kind_str, '(i0, a)') kind16, ': KIND16 [True RBend, EXACT_MODEL = T]'
case (KIND17);       write (kind_str, '(i0, a)') kind17, ': KIND17 [SixTrack Solenoid]'
case (KIND18);       write (kind_str, '(i0, a)') kind18, ': KIND18 [Rcollimator]'
case (KIND19);       write (kind_str, '(i0, a)') kind19, ': KIND19 [Ecollimator]'
case (KIND20);       write (kind_str, '(i0, a)') kind20, ': KIND20 [Straight Geometry MAD RBend]'
case (KIND21);       write (kind_str, '(i0, a)') kind21, ': KIND21 [Traveling Wave Cavity]'
case (KIND22);       write (kind_str, '(i0, a)') kind22, ': KIND22'
case (KIND23);       write (kind_str, '(i0, a)') kind23, ': KIND23'
case (KINDWIGGLER);  write (kind_str, '(i0, a)') kindwiggler, ': KINDWIGGLER [Wiggler]'
case (KINDPA);       write (kind_str, '(i0, a)') kindpa, ': KINDPA'
case default;        write (kind_str, '(a, i0, a)') 'UNKNOWN! [', this_kind, ']'
end select

end function kind_name

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine set_ptc (e_tot, particle, taylor_order, integ_order, n_step, &
!                      no_cavity, exact_modeling, exact_misalign, init_complex, force_init)
!
! Subroutine to initialize PTC.
!
! Note: At some point before you use PTC to compute Taylor maps etc.
!   you have to call set_ptc with both e_tot and particle args present. 
!   Always supply both of these args together or not at all. 
!
! Note: If you just want to use FPP without PTC then call the FPP routine init directly.
!
! Note: This subroutine cannot be used if you want to have "knobs" (in the PTC sense).
!
! Note: Use the routine get_ptc_params to get the state of PTC parameters.
!
! Note: Call this routine to transfer the value of the electric dipole moment from 
!   bmad_com%electric_dipole_moment to PTC.
!
! Input:
!   e_tot  -- Real(rp), optional: Energy in eV.
!   particle     -- Integer, optional: Type of particle:
!                     electron$, proton$, etc.
!   taylor_order -- Integer, optional: Maximum order of the taylor polynomials.
!                     0 => Use default.
!   integ_order  -- Integer, optional: Default Order for the drift-kick-drift 
!                     sympletic integrator. Possibilities are: 2, 4, or 6
!                     Default = 2
!   n_step       -- Integer, optional: Default Number of integration steps.
!                     Default = 1
!   no_cavity    -- Logical, optional: No RF Cavity exists? 
!                     Default = False.
!                     Corresponds to the nocavity option of the PTC init routine.
!                     no_cavity = .true. will turn any cavity into a drift.
!   exact_modeling -- logical, optional: Sets the PTC EXACT_MODEL variable.
!                       Default = False.
!                       See the PTC guide for more details.
!   exact_misalign -- logical, optional: Sets the PTC ALWAYS_EXACTMIS variable.
!                       Default = true.
!                       See the PTC guide for more details.
!   init_complex   -- logical, optional: If present and True then init complex PTC.
!                       Note: Complex PTC will also be initialized with bmad_com%spin_tracking_on = T.
!   force_init     -- logical, optional: If present and True then force a PTC init.
!-

subroutine set_ptc (e_tot, particle, taylor_order, integ_order, n_step, &
                        no_cavity, exact_modeling, exact_misalign, init_complex, force_init) 

use mad_like, only: make_states, exact_model, always_exactmis, pmaMUON, pmaE, &
              assignment(=), nocavity0, default, operator(+), in_bmad_units, &
              berz, init, set_madx, lp, superkill, TIME0, PHASE0, HIGHEST_FRINGE, init_all, SPIN0
use madx_ptc_module, only: ptc_ini_no_append, append_empty_layout, m_u, bmadl, use_info, &
              use_info_m, check_longitudinal
use c_tpsa, only: c_verbose, E_MUON, USE_QUATERNION

implicit none

integer, optional :: integ_order, particle, n_step, taylor_order
integer this_method, this_steps, t_order

real(rp), optional :: e_tot
real(rp), save :: old_e_tot = 0
real(dp) this_energy

logical, optional :: no_cavity, exact_modeling, exact_misalign, init_complex, force_init
logical, save :: init_ptc_needed = .true., init_init_needed = .true., init_spin_needed = .true.
logical params_present, c_verbose_save

character(16) :: r_name = 'set_ptc'

! ptc cannot be used with photons

if (logic_option(.false., force_init)) init_ptc_needed = .true.

if (present(particle)) then
  if (particle == photon$) return
endif

! Some init

USE_QUATERNION = .TRUE.
E_MUON = bmad_com%electric_dipole_moment
CHECK_LONGITUDINAL = .false. ! MAD-X uses the True setting.
call in_bmad_units

if (init_init_needed) then
  EXACT_MODEL = .false.
  ALWAYS_EXACTMIS = .true.
  init_init_needed = .false.
endif

! More init

HIGHEST_FRINGE = bmad_com%ptc_max_fringe_order

! do not call set_mad

params_present = present(e_tot) .and. present(particle)

if (init_ptc_needed .and. params_present) then
  if (particle == muon$ .or. particle == antimuon$) then
    call make_states (pmaMUON/pmaE)
  elseif (particle == positron$ .or. particle == electron$) then
    call make_states(.true._lp)
  elseif (particle == proton$ .or. particle == antiproton$) then
    if (particle /= proton$ .and. particle /= antiproton$) then
      call out_io (s_error$, r_name, 'PTC IS NOT ABLE TO HANDLE PARTICLES OF TYPE: ' // species_name(particle), 'USING PROTON/ANTIPROTON')
    endif
    call make_states(.false._lp)
  else
    call make_states (mass_of(particle)/mass_of(electron$), anomalous_moment_of(particle), real(charge_of(particle), dp))
    call out_io (s_warn$, r_name, 'Note: Radiation calculation in PTC not correct for particles of type: ' // species_name(particle))
  endif

  ! Use PTC time tracking
  DEFAULT = DEFAULT + TIME0
  PHASE0 = 0
endif

if (present (exact_modeling))     EXACT_MODEL = exact_modeling
if (present (exact_misalign))     ALWAYS_EXACTMIS = exact_misalign
if (present(no_cavity))           DEFAULT = DEFAULT + NOCAVITY0
if (bmad_com%spin_tracking_on)    DEFAULT = DEFAULT + SPIN0

if (present (integ_order)) then
  this_method = integ_order
  bmad_com%default_integ_order = integ_order
else
  this_method = bmad_com%default_integ_order
endif

if (present (n_step)) then
  this_steps = n_step
else
  this_steps = 10
endif

if (present(taylor_order)) then
  t_order = taylor_order
  if (t_order == 0) t_order = ptc_com%taylor_order_saved
  ptc_com%taylor_order_saved = t_order
endif

if (params_present) then
  if (init_ptc_needed .or. old_e_tot /= e_tot .or. present(integ_order) .or. present(n_step)) then
    this_energy = 1d-9 * e_tot
    if (this_energy == 0) then
      call out_io (s_fatal$, r_name, 'E_TOT IS 0.')
      if (global_com%exit_on_error) call err_exit
    endif
    call set_madx (energy = this_energy, method = this_method, step = this_steps)
    old_e_tot  = e_tot
    ! Only do this once
    if (init_ptc_needed) call ptc_ini_no_append 
    init_ptc_needed = .false.
  endif
endif

! Do not call init before the call to make_states
! Note: Once complex_ptc is set to True it remains True forever.

ptc_com%complex_ptc_used = ptc_com%complex_ptc_used .or. logic_option(.false., init_complex) .or. &
                                                                            bmad_com%spin_tracking_on 

if (.not. init_ptc_needed .or. logic_option(.false., force_init)) then  ! If make_states has been called
  t_order = 0
  if (present(taylor_order)) t_order = taylor_order
  if (t_order == 0) t_order = bmad_com%taylor_order
  if (t_order == 0) t_order = ptc_com%taylor_order_saved
  if (ptc_com%taylor_order_ptc /= t_order) then
    ! Due to Bmad vs PTC units bug, call init with nocavity
    call init (DEFAULT+NOCAVITY0, t_order, 0)
    init_spin_needed = .true.
    c_verbose_save = c_verbose
    c_verbose = .false.
    c_verbose = c_verbose_save
    ptc_com%taylor_order_ptc = t_order
  endif

  if (ptc_com%complex_ptc_used .and. init_spin_needed) then
    call init_all (DEFAULT, t_order, 0)
    init_spin_needed = .false.
  endif
endif

! Superkill tells PTC to do a through cleanup when killing a fibre.

SUPERKILL = .false.

!

use_info   = .true.    ! Enable storage space in fibre%i
use_info_m = .true.    ! Enable matrix storage in fibre%i%m

end subroutine set_ptc

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine get_ptc_params (ptc_param)
!
! Routine to return ptc parameters.
!
! Output:
!   ptc_param -- ptc_parameter_struct: PTC parameters.
!-

subroutine get_ptc_params (ptc_param)

use mad_like, only: EXACT_MODEL, ALWAYS_EXACTMIS, DEFAULT, PHASE0

implicit none

type (ptc_parameter_struct) ptc_param

!

ptc_param%exact_model     = EXACT_MODEL
ptc_param%exact_misalign  = ALWAYS_EXACTMIS
ptc_param%phase0          = PHASE0

end subroutine get_ptc_params

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine bmad_taylor_equal_real_8 (bmad_taylor, y8)
!
! Subroutine to convert from a real_8 taylor map in Etienne's PTC 
! to a taylor map in Bmad. This does not do any
! conversion between Bmad units (z, dp/p0) and PTC units (dE/p0, c*t).
!
! Subroutine overloads "=" in expressions
!       bmad_taylor = y8
!
! Input:
!   y8(:) -- real_8: PTC Taylor map.
!
! Output:
!   bmad_taylor(:) -- Taylor_struct: Input taylor map.
!-

subroutine bmad_taylor_equal_real_8 (bmad_taylor, y8)

use polymorphic_taylor, only: assignment (=), universal_taylor, real_8
use definition, only: real_8

implicit none

type (real_8), intent(in) :: y8(:)
type (taylor_struct), intent(inout) :: bmad_taylor(:)
type (universal_taylor) :: u_t(size(y8))

integer i, n_taylor

!

do i = 1, size(y8)
  u_t(i) = 0  ! nullify
  u_t(i) = y8(i)%t
  call universal_to_bmad_taylor (u_t(i), bmad_taylor(i))
  u_t(i) = -1  ! deallocate
enddo

end subroutine bmad_taylor_equal_real_8

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine real_8_equal_bmad_taylor (y8, bmad_taylor)
!
! Subroutine to convert from a taylor map in Bmad to a real_8 taylor map in Etienne's PTC. 
!
! Subroutine overloads "=" in expressions
!       y8 = bmad_taylor
!
! Input:
!   bmad_taylor(:) -- Taylor_struct: Input taylor map.
!
! Output:
!   y8(:) -- real_8: PTC Taylor map.
!-

subroutine real_8_equal_bmad_taylor (y8, bmad_taylor)

use polymorphic_taylor, only: kill, assignment(=), real_8, universal_taylor, alloc

implicit none

type (real_8), intent(inout) :: y8(:)
type (taylor_struct), intent(in) :: bmad_taylor(:)
type (universal_taylor) :: u_t

integer i, j, n, n_taylor

! init

call kill (y8)
call alloc(y8)

!

n_taylor = size(bmad_taylor)  ! Generally n_taylor = 6
do i = 1, n_taylor

  n = size(bmad_taylor(i)%term)
  allocate (u_t%n, u_t%nv, u_t%c(n), u_t%j(n,n_taylor))
  u_t%n = n
  u_t%nv = n_taylor

  do j = 1, n
    u_t%j(j,:) = bmad_taylor(i)%term(j)%expn(:)
    u_t%c(j) = bmad_taylor(i)%term(j)%coef
  enddo

  y8(i) = u_t
  u_t = -1   ! deallocate
enddo

end subroutine real_8_equal_bmad_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ptc_taylor_equal_bmad_taylor (ptc_taylor, bmad_taylor)
!
! Subroutine to convert from a taylor map in Bmad to a taylor map in Etienne's PTC. 
!
! Subroutine overloads "=" in expressions
!       ptc_taylor = bmad_taylor
!
! Input:
!   bmad_taylor(:) -- taylor_struct: Input taylor map.
!
! Output:
!   ptc_taylor(:)   -- taylor: PTC Taylor map.
!-

subroutine ptc_taylor_equal_bmad_taylor (ptc_taylor, bmad_taylor)

use polymorphic_taylor, only: alloc, kill, assignment(=), taylor, universal_taylor

implicit none

type (taylor), intent(inout) :: ptc_taylor
type (taylor_struct), intent(in) :: bmad_taylor
type (universal_taylor) :: u_t

integer j, n
character(*), parameter :: r_name = 'ptc_taylor_equal_bmad_taylor'

!

if (.not. associated(bmad_taylor%term)) then
  call out_io (s_error$, r_name, 'TAYLOR SERIES NOT DEFINED!')
  return
endif

n = size(bmad_taylor%term)
allocate (u_t%n, u_t%nv, u_t%c(n), u_t%j(n, 6))
u_t%n = n
u_t%nv = 6

do j = 1, n
  u_t%j(j,:) = bmad_taylor%term(j)%expn(:)
  u_t%c(j) = bmad_taylor%term(j)%coef
enddo

ptc_taylor = u_t
u_t = -1   ! deallocate

end subroutine ptc_taylor_equal_bmad_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine sigma_mat_ptc_to_bmad (sigma_mat_ptc, beta0, sigma_mat_bmad)
!
! Routine to convert a PTC sigma matrix to a Bmad sigma matrix.
! The conversion includes the conversion between Bmad and PTC time coordinate systems.
!
! Since PTC uses delta_E/P0c and Bmad uses delta_P/P0c coordinates, and since 
! the relationship between delta_E and delta_P is nonlinear, this routine 
! simplifies the calculation and assumes that the particle beta is constant
! over the range of particle energies.
!
! Input:
!   sigma_mat_ptc(6,6)  -- real(rp): PTC sigma matrix.
!   beta0               -- real(rp): Reference particle velocity
!
! Output:
!   sigma_mat_bmad(6,6) -- real(rp): Bmad sigma matrix.
!-

subroutine sigma_mat_ptc_to_bmad (sigma_mat_ptc, beta0, sigma_mat_bmad)

implicit none

real(rp) sigma_mat_bmad(6,6), sigma_mat_ptc(6,6)
real(rp) beta0, temp_mat(6,6)

!

temp_mat = sigma_mat_ptc
sigma_mat_bmad = temp_mat

sigma_mat_bmad(5,1:4) = temp_mat(6,1:4) * beta0
sigma_mat_bmad(6,1:4) = temp_mat(5,1:4) / beta0

sigma_mat_bmad(1:4,5) = temp_mat(1:4,6) * beta0
sigma_mat_bmad(1:4,6) = temp_mat(1:4,5) / beta0

sigma_mat_bmad(5,5) = temp_mat(6,6) * beta0**2
sigma_mat_bmad(6,6) = temp_mat(5,5) / beta0**2

end subroutine sigma_mat_ptc_to_bmad 

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Elemental subroutine universal_equal_universal (ut1, ut2)
!
! Subroutine to transfer the values in one universal taylor variable to
! another. Note: ut1 needs to have been initialized.
!
! Input:
!   ut2 -- Universal_taylor:
!
! Output:
!   ut1 -- Universal_taylor:
!-

elemental subroutine universal_equal_universal (ut1, ut2)

use definition, only: universal_taylor

implicit none

type (universal_taylor), intent(inout) :: ut1
type (universal_taylor), intent(in)    :: ut2

!

if (associated (ut1%n)) deallocate (ut1%n, ut1%nv, ut1%c, ut1%j)
allocate (ut1%n, ut1%nv, ut1%c(ut2%n), ut1%j(ut2%n, ut2%nv))

ut1%n  = ut2%n
ut1%nv = ut2%nv
ut1%c  = ut2%c
ut1%j  = ut2%j

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Elemental subroutine universal_to_bmad_taylor (u_taylor, bmad_taylor)
!
! Subroutine to convert from a universal_taylor map in Etienne's PTC 
! to a taylor map in Bmad.
!
! Input:
!   u_taylor(:) -- Universal_taylor: Universal_taylor map.
!
! Output:
!   bmad_taylor(:)   -- Taylor_struct:
!-

subroutine universal_to_bmad_taylor (u_taylor, bmad_taylor)

use definition, only: universal_taylor

implicit none

type (universal_taylor), intent(in) :: u_taylor
type (taylor_struct), intent(inout) :: bmad_taylor

integer :: j, k, n

! Remember to suppress any terms that have a zero coef.  

n = count(u_taylor%c(:) /= 0)

if (associated(bmad_taylor%term)) then
  if (size(bmad_taylor%term) /= n) deallocate(bmad_taylor%term)
endif
if (.not. associated(bmad_taylor%term)) allocate(bmad_taylor%term(n))

k = 0
do j = 1, u_taylor%n
  if (u_taylor%c(j) == 0) cycle
  k = k + 1
  bmad_taylor%term(k)%expn = u_taylor%j(j,:)
  bmad_taylor%term(k)%coef = u_taylor%c(j)
enddo

end subroutine universal_to_bmad_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine complex_taylor_equal_c_taylor (bmad_complex_taylor, ptc_c_taylor)
!
! Subroutine to convert a PTC c_taylor to a bmad_complex_taylor
!
! Subroutine overloads "=" in expressions
!       bmad_complex_taylor = ptc_c_taylor
! 
!
!
! Input:
!   ptc_c_taylor -- c_taylor: PTC complex Taylor series
!
! Output:
!   bmad_complex_taylor -- complex_taylor_struct: Bmad complex Taylor series
!-

subroutine complex_taylor_equal_c_taylor (bmad_complex_taylor, ptc_c_taylor)

use polymorphic_taylor, only: assignment (=), alloc, universal_taylor, real_8, &
                              complextaylor, c_taylor, kill
use c_TPSA, only: assignment(=)

implicit none

type (c_taylor), intent(in) :: ptc_c_taylor
type (complextaylor) :: ct
type (complex_taylor_struct), intent(inout) :: bmad_complex_taylor
type (universal_taylor) :: re_u, im_u
type (taylor_struct) :: re, im

integer :: i, n_taylor

! Convert to complextaylor
call alloc(ct)
ct = ptc_c_taylor

! Convert to two universal_taylor
re_u = 0
im_u = 0
re_u = ct%r
im_u = ct%i
call universal_to_bmad_taylor (re_u, re)
call universal_to_bmad_taylor (im_u, im)

! Form 
call form_complex_taylor(re, im, bmad_complex_taylor)

! Cleanup
call kill(ct)
re_u = -1
im_u = -1

end subroutine complex_taylor_equal_c_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!  Subroutine complex_taylors_equal_c_taylors (bmad_complex_taylors(:), ptc_c_taylors(:))
!  Vector version of complex_taylor_equal_c_taylor
!
!  Subroutine overloads '=' for bmad_complex_taylor(:) = ptc_c_taylor(:)
!
!-
subroutine complex_taylors_equal_c_taylors (bmad_complex_taylors, ptc_c_taylors)
use polymorphic_taylor, only: c_taylor
implicit none
type (c_taylor), intent(in) :: ptc_c_taylors(:)
type (complex_taylor_struct), intent(inout) :: bmad_complex_taylors(:)
integer :: i
!
do i = 1, size(ptc_c_taylors)
  bmad_complex_taylors(i) = ptc_c_taylors(i)
enddo

end subroutine complex_taylors_equal_c_taylors



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine bmad_taylor_equal_ptc_taylor (bmad_taylor, ptc_taylor)
!
! Subroutine to convert a PTC taylor to a bmad_taylor
!
! Subroutine overloads "=" in expressions
!       bmad_taylor = ptc_taylor
!
! Input:
!   ptc_taylor -- taylor: PTC complex Taylor series
!
! Output:
!   bmad_taylor -- taylor_struct: Bmad Taylor series
!-

subroutine bmad_taylor_equal_ptc_taylor (bmad_taylor, ptc_taylor)

use polymorphic_taylor, only: assignment (=), universal_taylor,  taylor

implicit none

type (taylor), intent(in) :: ptc_taylor
type (taylor_struct), intent(inout) :: bmad_taylor
type (universal_taylor) :: utaylor

! Convert to universal_taylor, then to bmad_taylor
utaylor = 0
utaylor = ptc_taylor
call universal_to_bmad_taylor (utaylor, bmad_taylor)
! Cleanup
utaylor = -1

end subroutine bmad_taylor_equal_ptc_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!  Subroutine complex_taylors_equal_c_taylors (bmad_complex_taylors(:), ptc_c_taylors(:))
!  Vector version of bmad_taylor_equal_ptc_taylor
!
!  Subroutine overloads '=' for bmad_taylor(:) = ptc_taylor(:)
!
!-
subroutine bmad_taylors_equal_ptc_taylors (bmad_taylors, ptc_taylors)
use polymorphic_taylor, only: taylor
implicit none
type (taylor), intent(in) :: ptc_taylors(:)
type (taylor_struct), intent(inout) :: bmad_taylors(:)
integer :: i
!
do i = 1, size(ptc_taylors)
  bmad_taylors(i) = ptc_taylors(i)
enddo

end subroutine bmad_taylors_equal_ptc_taylors


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine form_complex_taylor(re_taylor, im_taylor, complex_taylor)
!
! Subroutine to form a complex taylor from two taylor series representing 
!   the real and imaginary parts
!
! Input:
!   re_taylor       -- taylor_struct: Real part
!   im_taylor       -- taylor_struct: Imaginary part
!
! Output:
!   complex_taylor  -- complex_taylor_struct: combined complex taylor
!-
subroutine form_complex_taylor(re_taylor, im_taylor, complex_taylor)

implicit none

type (complex_taylor_struct) :: complex_taylor
type (taylor_struct) :: re_taylor, im_taylor
type (taylor_struct) :: taylor1, taylor2
real(rp) :: re, im
integer :: n, n1, n2, n_tot, t1, t2, ix1, ix2, expn(6)
!

n1 = size(re_taylor%term)
n2 = size(im_taylor%term)

! Easy cases
if (n1 == 0 .or. n2 ==0 ) then 
  
  if (n1 > 0) then
  ! purely real
     call init_complex_taylor_series (complex_taylor, n1)
     complex_taylor%ref = cmplx(re_taylor%ref, 0.0_rp, rp)
     do n = 1, n1
       complex_taylor%term(n)%coef = cmplx(re_taylor%term(n)%coef, 0.0_rp, rp)
       complex_taylor%term(n)%expn = re_taylor%term(n)%expn
     enddo
 
  else if (n2 > 0) then
  ! purely imaginary
    call init_complex_taylor_series (complex_taylor, n2)
    complex_taylor%ref = cmplx(0.0_rp, im_taylor%ref, rp)
    do n = 1, n2
      complex_taylor%term(n)%coef = cmplx(0.0_rp, im_taylor%term(n)%coef, rp)
      complex_taylor%term(n)%expn = im_taylor%term(n)%expn
    enddo
  
  else if (n1 == 0 .and. n2 == 0 ) then
    ! No terms!
    call init_complex_taylor_series (complex_taylor, 0)
  endif
  
  return
endif

! Init scratch series for sort
call init_taylor_series (taylor1, n1)
call init_taylor_series (taylor2, n2)

! Sort
call sort_taylor_terms (re_taylor, taylor1)
call sort_taylor_terms (im_taylor, taylor2)

! First pass to count unique terms
t1 = 1
t2 = 1
n_tot = 1
ix1 = taylor_exponent_index(taylor1%term(t1)%expn)
ix2 = taylor_exponent_index(taylor2%term(t2)%expn)
do
  if (t1 == n1 .and. t2 == n2) exit
  n_tot = n_tot + 1
  if (ix1 == ix2) then
    ! re and im parts      
    call iter1()
    call iter2()
  else if (ix1 < ix2) then
    ! re term only
    call iter1()
  else
    ! im term only
    call iter2()
  endif
end do

! Initialize output taylor
call init_complex_taylor_series(complex_taylor, n_tot)
complex_taylor%ref = cmplx(taylor1%ref, taylor2%ref, rp)

! Second pass to assign values
n = 1
t1 = 1
t2 = 1
ix1 = taylor_exponent_index(taylor1%term(t1)%expn)
ix2 = taylor_exponent_index(taylor2%term(t2)%expn)
do
  if (ix1 == ix2) then
    ! re and im parts  
    expn = taylor1%term(t1)%expn
    re = taylor1%term(t1)%coef
    im = taylor2%term(t2)%coef
    call iter1()
    call iter2()
  
  else if (ix1 < ix2) then
    ! re term only
    expn = taylor1%term(t1)%expn
    re = taylor1%term(t1)%coef
    im = 0.0_rp
    call iter1()
  
  else
    ! im term only
    expn = taylor2%term(t2)%expn
    re = 0.0_rp
    im = taylor2%term(t2)%coef
    call iter2()
  endif
  
  ! Assign complex coef
  complex_taylor%term(n)%coef = cmplx(re, im, rp)
  complex_taylor%term(n)%expn = expn
  if (n == n_tot) exit
  n = n + 1
end do

! Cleanup
deallocate(taylor1%term)
deallocate(taylor2%term)

contains
  ! Stepping helper routines
  subroutine iter1()
  implicit none
    if (t1 < n1) then
      t1 = t1 + 1 
      ix1 = taylor_exponent_index(taylor1%term(t1)%expn)
    else
      ix1 = huge(0) ! Set to largest integer
    endif
  end subroutine

  subroutine iter2()
  implicit none
    if (t2 < n2) then
      t2 = t2 + 1 
      ix2 = taylor_exponent_index(taylor2%term(t2)%expn)
    else
      ix2 = huge(0)
    endif
  end subroutine

end subroutine



!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine concat_real_8 (y1, y2, y3, r2_ref, keep_y1_const_terms)
!
! Subroutine to concatinate two real_8 taylor maps.
!       y3 = y2(y1)
! This subroutine assumes that y1, y2, and y3 have been allocated.
!
! For any map "M", Bmad treats the argument of the map as M(r-r_ref) where
! r are the coordinates and r_ref are the reference coordinates at the beginning of M.
! The reference coordinates at the end of M is thus M(0).
!
! If r2_ref is not present, it is assumed that the referece orbit at the end of 
! y1 (which is equal to y1(0)), is equal to the referene orbit at the beginning of y2.
! If r2_ref is present, the maps are "aligned" so that
!   y3(r - r1_ref) = y2(y1(r - r1_ref) - r2_ref)
! What this means in practice is that, if r2_ref is not present, it is assumed that 
! r2_ref = y1(0) so that the constant terms y1(0) are thrown away. 
!
! Input:
!   y1(:)     -- real_8: First Input map.
!   y2(:)     -- real_8: Second Input map.
!   r2_ref(:) -- real(dp), optional: Reference orbit at beginning of y2. See above.
!                 Cannot be present if keep_y1_const_terms is present
!   keep_y1_const_terms
!             -- logical, optional: If present and True, just concatenate y1 and y2 retaining any constant terms in y1.
!                 That is, ignore the reference orbit. Cannot be present if r2_ref is present.
!
! Output
!   y3(:) -- real_8: Concatinated map.
!-

subroutine concat_real_8 (y1, y2, y3, r2_ref, keep_y1_const_terms)

use s_fitting, only: alloc, assignment(=), kill, damap, operator(*), operator(+), real_8, operator(.o.), operator(.sub.), print

implicit none

type (real_8) :: y1(:), y2(:), y3(:)
type (damap) da1, da2, da3, id
real(dp), optional :: r2_ref(:)

integer i
logical, optional :: keep_y1_const_terms

! set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

! Allocate temp vars

call alloc(id)
call alloc(da1)
call alloc(da2)
call alloc(da3)

! concat

da1 = y1
da2 = y2

if (logic_option(.false., keep_y1_const_terms)) then
  da3 = da2 .o. da1

elseif (present(r2_ref)) then
  id = 1
  y3 = id + (-r2_ref) ! Identity map - r2_ref
  da3 = y3
  da1 = da3 .o. da1
  da3 = da2 .o. da1  ! Includes constant term of da1.

else
  da3 = da2 * da1
endif

y3 = da3

! kill temp vars

call kill (id)
call kill (da1)
call kill (da2)
call kill (da3)

end subroutine concat_real_8

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_to_genfield (bmad_taylor, ptc_genfield, c0)
!
! Subroutine to construct a genfield (partially inverted map) from a taylor
! map.
! Note: The constant terms of the taylor map are removed in the process.
! Note: The genfield uses PTC coordinates.
!
! Moudules needed:
!   use ptc_interface_mod
!
! Input:
!   bmad_taylor(:) -- Taylor_struct: Input taylor map.
!
! Output:
!   ptc_genfield   -- Genfield: Output partially inverted map.
!   c0(:)          -- Real(rp): The constant part of the bmad_taylor map
!-

subroutine taylor_to_genfield (bmad_taylor, ptc_genfield, c0)

use s_fitting, only: alloc, kill, assignment(=), damap, real_8

implicit none

type (taylor_struct), intent(in) :: bmad_taylor(:)
type (genfield), intent(inout) :: ptc_genfield
type (taylor_struct) taylor(size(bmad_taylor))
type (damap) da_map
type (real_8) y(size(bmad_taylor))

real(rp), intent(out) :: c0(:)

! set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

! Remove constant terms from the taylor map first. This is probably
! not needed but we do it to make sure everything is alright.
! Also remove terms that have higher order then bmad_com%taylor_order

call remove_constant_taylor (bmad_taylor, taylor, c0, .true.)

! allocate pointers

call alloc (ptc_genfield)
call alloc (da_map)
call alloc (y)

! calculate the ptc_genfield

y = taylor
da_map = y
ptc_genfield = da_map

! cleanup

call kill (da_map)
call kill (y)
call kill_taylor (taylor)

end subroutine taylor_to_genfield

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine remove_constant_taylor (taylor_in, taylor_out, c0, remove_higher_order_terms)
!
! Subroutine to remove the constant part of a taylor map.
! Optionally terms that are higher order than bmad_com%taylor_order can
! be removed.
!
! Note: It is assumed that taylor_out has been deallocated before the call to
! this routine. Calling this routine with the first two actual arguments the
! same is prohibited.
!
! Moudules needed:
!   use ptc_interface_mod
!
! Input:
!  taylor_in(:)              -- Taylor_struct: Input taylor map.
!  remove_higher_order_terms -- Logical: If True then terms that are higher
!                               order than bmad_com%taylor_order are removed.
!
! Output:
!   taylor_out(:)  -- Taylor_struct: Taylor with constant terms removed.
!   c0(:)          -- Real(rp): The constant part of the taylor map
!-

subroutine remove_constant_taylor (taylor_in, taylor_out, c0, remove_higher_order_terms)

implicit none

type (taylor_struct), intent(in) :: taylor_in(:)
type (taylor_struct) taylor_out(:)

real(rp), intent(out) :: c0(:)

integer i, j, n, nn, ss

logical, intent(in) :: remove_higher_order_terms

!

c0 = 0

do i = 1, size(taylor_in)

  if (.not. associated(taylor_in(i)%term)) cycle
  n = size(taylor_in(i)%term)

  do j = 1, size(taylor_in(i)%term)
    if (all(taylor_in(i)%term(j)%expn == 0)) then
      n = n - 1
      c0(i) = taylor_in(i)%term(j)%coef
    endif
    if (remove_higher_order_terms) then
      if (sum(taylor_in(i)%term(j)%expn) > ptc_com%taylor_order_ptc) n = n - 1
    endif
  enddo

  if (associated(taylor_out(i)%term)) then
    if (size(taylor_out(i)%term) /= n) deallocate(taylor_out(i)%term)
  endif
  if (.not. associated(taylor_out(i)%term)) allocate (taylor_out(i)%term(n))

  nn = 0
  do j = 1, size(taylor_in(i)%term)
    ss = sum(taylor_in(i)%term(j)%expn)
    if (ss == 0 .or. (remove_higher_order_terms .and. ss > ptc_com%taylor_order_ptc)) cycle
    nn = nn + 1
    taylor_out(i)%term(nn) = taylor_in(i)%term(j)
  enddo

enddo

end subroutine remove_constant_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_inverse (taylor_in, taylor_inv, err)
!
! Subroutine to invert a taylor map. Since the inverse map is truncated, it is not exact.
!
! Moudules needed:
!   use ptc_interface_mod
!
! Input:
!   taylor_in(:)  -- Taylor_struct: Input taylor map.
!
! Output:
!   taylor_inv(:) -- Taylor_struct: Inverted taylor map.
!   err           -- Logical, optional: Set True if there is no inverse.
!                     If not present then print an error message.
!-

subroutine taylor_inverse (taylor_in, taylor_inv, err)

use s_fitting, only: assignment(=), alloc, kill, operator(**), damap, real_8

implicit none

type (taylor_struct) :: taylor_in(:)
type (taylor_struct) :: taylor_inv(:)
type (taylor_struct) tlr(size(taylor_in)), tlr2(size(taylor_in))
type (real_8) y(size(taylor_in))
type (damap) da

real(rp) c0(size(taylor_in))

integer i, n_taylor
integer :: expn0(6) = 0

real(dp) c8(size(taylor_in)), c_ref(size(taylor_in))

logical, optional :: err

character(16) :: r_name = 'taylor_inverse'

! Set the taylor order in PTC if not already done so.

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

call alloc(da)
call alloc(y)

! The inverse operation of PTC ignores constant terms so we have 
! to take them out and then put them back in.

call remove_constant_taylor (taylor_in, tlr, c0, .true.)

! Each taylor_in(i) must have at least one term for an inverse to exist.

n_taylor = size(taylor_in)
do i = 1, n_taylor
  if (.not. associated(tlr(i)%term) .or. size(tlr(i)%term) == 0) then
    if (present(err)) then
      err = .true.
    else
      call out_io (s_error$, r_name, 'Taylor map does not have an inverse.')
    endif
    return
  endif
enddo

! Compute inverse.

y = tlr
da = y
da = da**(-1)
y = da

! Transfer inverse to taylor_inv.
! If the Map is written as:
!   R_out = T' (R_in - R_ref) + C
! Where C = constant part of map and T' = T - C 
! Then:
!   R_in = T'^-1 (R_out - C) + R_ref
! Thus:
!   T_inv = T'^-1 + R_ref
!   R_ref_inv = C

taylor_inv = y
taylor_inv%ref = c0

do i = 1, size(taylor_in)
  call add_taylor_term(taylor_inv(i), taylor_in(i)%ref, expn0)
enddo

! Clean up

call kill (da)
call kill (y)
call kill_taylor (tlr)

if (present(err)) err = .false.

end subroutine taylor_inverse

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine concat_taylor (taylor1, taylor2, taylor3)
! 
! Subroutine to concatinate two taylor maps:
!   taylor3[x] = taylor2(taylor1[x])
!
! Note: In general, if taylor2 is a component of an ele_struct, use
! concat_ele_taylor instead.
!
! Input:
!   taylor1(:) -- Taylor_struct: Taylor map.
!   taylor2(:) -- Taylor_struct: Taylor map.
!
! Output
!   taylor3(:) -- Taylor_struct: Concatinated map
!-

subroutine concat_taylor (taylor1, taylor2, taylor3)

use s_fitting, only: assignment(=), alloc, kill, real_8

implicit none

type (taylor_struct) :: taylor1(:), taylor2(:)
type (taylor_struct) :: taylor3(:)
type (real_8) y1(size(taylor1)), y2(size(taylor2)), y3(size(taylor3))

! Set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

! Allocate temp vars

call alloc (y1)
call alloc (y2)
call alloc (y3)

! Concat

y1 = taylor1
y2 = taylor2

call concat_real_8 (y1, y2, y3, taylor2%ref)

taylor3 = y3
taylor3(:)%ref = taylor1(:)%ref

! Cleanup

call kill (y1)
call kill (y2)
call kill (y3)

end subroutine concat_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine concat_ele_taylor (taylor1, ele, taylor3)
!
! Routine to concatinate two taylor maps:
!   taylor3[x] = ele_taylor(taylor1[x])
! If ele%taylor_map_includes_offsets = True:  ele_taylor == ele%taylor 
! If ele%taylor_map_includes_offsets = False: ele_taylor == ele%taylor + offset corrections. 
!
! Also see: concat_taylor
!
! Input:
!   taylor1(6) -- Taylor_struct: Taylor map.
!   ele        -- ele_struct: Element containing a Taylor map.
!
! Output
!   taylor3(6) -- Taylor_struct: Concatinated map
!-

Subroutine concat_ele_taylor (taylor1, ele, taylor3)

use s_tracking, only: mis_fib, alloc, kill, dtiltd, assignment(=), default, real_8, fibre

implicit none

type (ele_struct) ele
type (taylor_struct) taylor1(:), taylor3(:)
type (lat_param_struct) param
type (real_8) x_ele(6), x_body(6), x1(6), x3(6)
type (fibre), pointer :: fib

real(rp) beta0, beta1, tilt
real(8) x_dp(6)

! Match elements are not implemented in PTC so just use the matrix.

if (ele%key == match$) then
  call mat6_to_taylor (ele%vec0, ele%mat6, ele%taylor)
  call concat_taylor (taylor1, ele%taylor, taylor3)
  return
endif

! ele%taylor_map_includes_offsets = T means that misalignment effects are already included 
! in ele%taylor.

if (ele%taylor_map_includes_offsets .or. ele%key == taylor$) then
  call concat_taylor (taylor1, ele%taylor, taylor3)
  return
endif

! Here when we need to include the misalignment effects.
! First set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

! Create a PTC fibre that holds the misalignment info
! and create map corresponding to ele%taylor.

param%particle = positron$  ! Actually this does not matter to the calculation
call ele_to_fibre (ele, fib, param, use_offsets = .true.)

! Init

call alloc(x_ele)
call alloc(x_body)
call alloc(x1)
call alloc(x3)

beta0 = ele%value(p0c_start$)/ele%value(e_tot_start$)
beta1 = ele%value(p0c$)/ele%value(e_tot$)

x_dp = 0
x_ele = x_dp  ! x_ele = Identity map 

if (ele%key == sbend$) then
  tilt = ele%value(ref_tilt_tot$)
else
  tilt = ele%value(tilt_tot$)
endif

call dtiltd (tilt, 1, x_ele)
call mis_fib (fib, x_ele, DEFAULT, .true., entering = .true.)

x_body = ele%taylor

call concat_real_8 (x_ele, x_body, x_ele)
call mis_fib (fib, x_ele, DEFAULT, .true., entering = .false.)
call dtiltd (tilt, 2, x_ele)

! Concat with taylor1

x1 = taylor1
call concat_real_8 (x1, x_ele, x3, ele%taylor%ref)

! convert x3 to final result taylor3

taylor3(:)%ref = taylor1(:)%ref
taylor3 = x3

! Cleanup

call kill(x_ele)
call kill(x_body)
call kill(x1)
call kill(x3)

end subroutine concat_ele_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_to_real_8 (bmad_taylor, beta0, beta1, ptc_re8, ref_orb_ptc, exi_orb_ptc)
!
! Routine to convert a Bmad Taylor map to PTC real_8 map.
! The constant term of ptc_re8 will be removed.
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Input taylor map.
!   beta0          -- real(rp): Reference particle velocity at beginning of map
!   beta1          -- real(rp): Reference particle velocity at end of map
!   ptc_re8(6)     -- real_8: PTC Taylor polymorph. Note: it is the duty of the calling routine
!                       to call alloc beforehand and kill afterwards.
!
! Output:
!   ptc_re8(6)      -- real_8: PTC Taylor polymorph. Note: it is the duty of the calling routine
!                       to call alloc beforehand and kill afterwards.
!   ref_orb_ptc(6)  -- real(rp): PTC starting reference orbit.
!   exi_orb_ptc(6)  -- real(rp): Constant part of the map = orbit at the exit end.
!-

subroutine taylor_to_real_8 (bmad_taylor, beta0, beta1, ptc_re8, ref_orb_ptc, exi_orb_ptc)

use ptc_spin

implicit none

type (taylor_struct) :: bmad_taylor(6)

type(real_8) :: ptc_re8(6) ! It is the duty of the calling routine to alloc and kill this.

type(real_8), allocatable :: expn(:, :)
type(real_8) diff_orb(6), start_orb(6)
type(damap) id

real(rp) beta0, beta1, ref_ptc(6)
real(rp), optional :: exi_orb_ptc(6), ref_orb_ptc(6)
integer i, j, k, ie, e_max

!

ref_ptc = bmad_taylor%ref

call alloc(diff_orb)
call alloc(start_orb)
call alloc(id)

id = 1
start_orb = id + ref_ptc
ref_orb_ptc = ref_ptc

do i = 1, 6
  diff_orb(i) = start_orb(i) - bmad_taylor(i)%ref
enddo

! size cache matrix

e_max = 0 

do i = 1, 6
  do j = 1, size(bmad_taylor(i)%term)
    e_max = max (e_max, maxval(bmad_taylor(i)%term(j)%expn)) 
  enddo
enddo

allocate (expn(0:e_max, 6))
do i = 0, e_max
do j = 1, 6
  call alloc (expn(i, j))
enddo
enddo

! Fill in cache matrix

do i = 1, 6
  expn(0, i) = 1.0d0
enddo

do j = 1, e_max
do i = 1, 6
  expn(j, i) = expn(j-1, i) * diff_orb(i)
enddo
enddo

! Compute taylor map

call alloc (ptc_re8)

do i = 1, 6
  ptc_re8(i) = 0
  do j = 1, size(bmad_taylor(i)%term)
    ptc_re8(i) = ptc_re8(i) + bmad_taylor(i)%term(j)%coef * &
                       expn(bmad_taylor(i)%term(j)%expn(1), 1) * &
                       expn(bmad_taylor(i)%term(j)%expn(2), 2) * &
                       expn(bmad_taylor(i)%term(j)%expn(3), 3) * &
                       expn(bmad_taylor(i)%term(j)%expn(4), 4) * &
                       expn(bmad_taylor(i)%term(j)%expn(5), 5) * &
                       expn(bmad_taylor(i)%term(j)%expn(6), 6)
  enddo
enddo

! Remove constant.

exi_orb_ptc = ptc_re8
do i = 1, 6
  ptc_re8(i) = ptc_re8(i) - exi_orb_ptc(i)
enddo

! Cleanup.

call kill(diff_orb)
call kill(start_orb)
call kill(id)

do i = 0, e_max
do j = 1, 6
  call kill(expn(i, j))
enddo
enddo

end subroutine taylor_to_real_8

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_propagate1 (bmad_taylor, ele, param, track_particle)
!
! Subroutine to track (symplectic integration) a taylor map through an element.
! The alternative routine, if ele has a taylor map, is concat_taylor.
!
! This routine will fail if there is no corresponding ptc fibre for this
! element. In general, the transfer_map_calc routine should be used instead.
!
! Input:
!   bmad_taylor(6)   -- Taylor_struct: Map to be tracked
!   ele              -- Ele_struct: Element to track through
!   param            -- lat_param_struct: 
!   track_particle   -- coord_struct, optional: Particle to be tracked. ref_particle$, electron$, etc.
!                         Must be present if the particle to be tracked is not the reference particle or
!                         if the direction of propagation is backwards.
!
! Output:
!   bmad_taylor(6)  -- Taylor_struct: Map through element
!-

subroutine taylor_propagate1 (bmad_taylor, ele, param, track_particle)

use s_tracking
use mad_like, only: real_8, fibre, ptc_track => track
use ptc_spin, only: track_probe_x
use madx_ptc_module, only: bmadl

implicit none

type (taylor_struct) bmad_taylor(:)
type (real_8), save :: ptc_tlr(6)
type (ele_struct) ele, drift_ele
type (lat_param_struct) param
type (fibre), pointer :: ptc_fibre
type (coord_struct), optional :: track_particle

real(rp) beta0, beta1, m2_rel

! If the element is a taylor then just concat since this is faster.

if (ele%key == taylor$) then
  call concat_ele_taylor (bmad_taylor, ele, bmad_taylor)
  return
endif

! set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

! Init ptc map with bmad map

beta0 = ele%value(p0c_start$) / ele%value(e_tot_start$)
beta1 = ele%value(p0c$) / ele%value(e_tot$)

call alloc (ptc_tlr)
ptc_tlr = bmad_taylor

! Track entrance drift if PTC is using a hard edge model

if (tracking_uses_end_drifts(ele)) then
  call create_hard_edge_drift (ele, upstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, ptc_fibre, param, .true.)
  call ptc_track (ptc_fibre, ptc_tlr, default)  ! "track" in PTC
endif

! track the map

call ele_to_fibre (ele, ptc_fibre, param, .true., track_particle = track_particle)
call track_probe_x (ptc_tlr, DEFAULT, fibre1 = bmadl%start)

! Track exit side drift if PTC is using a hard edge model

if (tracking_uses_end_drifts(ele)) then
  call create_hard_edge_drift (ele, downstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, ptc_fibre, param, .true.)
  call ptc_track (ptc_fibre, ptc_tlr, default)  ! "track" in PTC
endif

! transfer ptc map back to bmad map

bmad_taylor = ptc_tlr

! cleanup

call kill (ptc_tlr)

end subroutine taylor_propagate1

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ele_to_taylor (ele, param, orb0, taylor_map_includes_offsets, orbital_taylor, spin_taylor)
!
! Subroutine to make orbital and spin (if spin tracking is on) taylor maps for an element. 
! The order of the map is set by set_ptc
!
! Input:
!   ele     -- Element_struct: 
!     %value(integrator_order$)    -- Order for the symplectic integrator: 2, 4, or 6.
!     %value(ds_step$)             -- Integrater step size.
!     %taylor_map_includes_offsets -- Make Taylor map with element offsets, pitches, and tilt?
!   orb0    -- Coord_struct, optional: Starting coords around which the Taylor map is evaluated.
!   param   -- lat_param_struct: 
!     %e_tot   -- Needed for wigglers.
!   taylor_map_includes_offsets -- Logical, optional: If present then overrides 
!                         ele%taylor_map_includes_offsets.
!
! Output:
!   orbital_taylor(6) -- taylor_struct, optional: Orbital taylor map.
!                         If not present then the map is put in ele%taylor.
!   spin_taylor(0:3)  -- taylor_struct, optional: Spin taylor map. 
!                         If not present then the map is put in ele%spin_taylor.
!-

subroutine ele_to_taylor (ele, param, orb0, taylor_map_includes_offsets, orbital_taylor, spin_taylor)

use s_tracking
use mad_like, only: real_8, fibre, ring_l, survey, make_node_layout, CONVERSION_XPRIME_IN_ABELL
use ptc_spin, only: track_probe_x, track_probe
use s_family, only: survey
use madx_ptc_module, only: bmadl

implicit none

type (ele_struct), target :: ele, drift_ele
type (lat_param_struct) :: param
type (coord_struct), optional, intent(in) :: orb0
type (coord_struct) c0
type (taylor_struct), optional, target :: orbital_taylor(6), spin_taylor(0:3)
type (taylor_struct), pointer :: orb_tylr(:), spin_tylr(:)
type (probe) ptc_probe
type (probe_8) ptc_probe8
type (fibre), pointer :: ptc_fibre
type (real_8) y2(6)
type (c_damap) ptc_cdamap

real(dp) x(6), beta
integer i, print12

logical, optional :: taylor_map_includes_offsets
logical use_offsets, err_flag

character(16) :: r_name = 'ele_to_taylor'

!

CONVERSION_XPRIME_IN_ABELL = (.not. bmad_com%convert_to_kinetic_momentum) ! Only affects cylindrical map eles

if (present(orbital_taylor)) then
  orb_tylr => orbital_taylor
else
  orb_tylr => ele%taylor
endif

if (present(spin_taylor)) then
  spin_tylr => spin_taylor
else
  spin_tylr => ele%spin_taylor
endif

! Match elements and helical wiggler/undulators without a map are not implemented in PTC so just use the matrix.

if (ele%key == match$ .or. &
            (.not. associated(ele%cylindrical_map) .and. .not. associated(ele%cartesian_map) .and. &
            .not. associated(ele%taylor_field) .and. (ele%key == wiggler$ .or. ele%key == undulator$) .and. &
            ele%field_calc == helical_model$)) then
  call mat6_to_taylor (ele%vec0, ele%mat6, orb_tylr)
  return
endif

! Init. Note that the fibre must be made before any map manipulation in case ele contains a taylor_field.

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

use_offsets = logic_option(ele%taylor_map_includes_offsets, taylor_map_includes_offsets)

if (tracking_uses_end_drifts(ele)) then
  call create_hard_edge_drift (ele, upstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, ptc_fibre, param, .true.)  ! Do not kill ptc_fibre2
  call ele_to_fibre (ele, ptc_fibre, param, use_offsets, track_particle = orb0, kill_layout = .false.)
  call create_hard_edge_drift (ele, downstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, ptc_fibre, param, .true., kill_layout = .false.)
  bmadl%closed = .true.
  call ring_l (bmadl, bmadl%closed)
  call survey (bmadl)
  call set_ptc_quiet (12, set$, print12)
  call make_node_layout(bmadl)
  call set_ptc_quiet (12, unset$, print12)
else
  call ele_to_fibre (ele, ptc_fibre, param, use_offsets, track_particle = orb0)
endif

call alloc(ptc_cdamap)
call alloc(ptc_probe8)

call attribute_bookkeeper (ele, .true.)

! Initial map

if (present(orb0)) then
  orb_tylr(:)%ref = orb0%vec
  x = orb0%vec
else
  orb_tylr(:)%ref = 0
  x = 0
endif

ptc_probe = 0
ptc_probe = x

ptc_cdamap = 1
ptc_probe8 = ptc_cdamap + ptc_probe ! = IdentityMap + const

! 

if (bmad_com%spin_tracking_on) then
  call track_probe (ptc_probe8, DEFAULT+SPIN0, fibre1 = bmadl%start)
else
  call track_probe (ptc_probe8, DEFAULT-SPIN0, fibre1 = bmadl%start)
endif

do i = 0, 3
  spin_tylr(i) = ptc_probe8%q%x(i)%t
enddo

! take out the offset

!if (any(x /= 0)) then
!  call alloc(y2)
!  y2 = -x  ! y2 = IdentityMap - x
!  call concat_real_8 (y2, y0, y0)
!  call kill(y2)
!endif

! convert to orb_tylr_struct

orb_tylr = ptc_probe8%x

call kill (ptc_probe8)
call kill (ptc_cdamap)

if (associated (ele%ptc_genfield%field)) call kill_ptc_genfield (ele%ptc_genfield%field)

call set_ele_status_stale (ele, mat6_group$)

CONVERSION_XPRIME_IN_ABELL = .true. ! Reset to normal.

end subroutine ele_to_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine type_real_8_taylors (y)
!
! Subroutine to type out the taylor map from a real_8 array.
!
! Input
!   y(6)     -- Real_8: 6 taylor map: 
!-

subroutine type_real_8_taylors (y)

use definition, only: real_8

implicit none

type (real_8) y(:)
type (taylor_struct) b_taylor(6)

!

b_taylor = y
call type_taylors (b_taylor)
call kill_taylor (b_taylor)

end subroutine type_real_8_taylors

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine sort_universal_terms (ut_in, ut_sorted)
!
! Subroutine to sort the taylor terms from "lowest" to "highest".
! This subroutine is needed since what comes out of PTC is not sorted.
!
! The number associated with a taylor_term that is used for the sort is:
!     number = sum(exp(i))*10^6 + exp(6)*10^5 + ... + exp(1)*10^0
! Where exp(1) is the exponent for x, exp(2) is the exponent for P_x, etc.
!
! Note: ut_sorted needs to have been initialized.
! Note: ut_sorted cannot be ut_in. That is it is not legal to write:
!           call sort_universal_terms (this_ut, this_ut)
!
! Input:
!   ut_in     -- Universal_taylor: Unsorted taylor series.
!
! Output:
!   ut_sorted -- Universal_taylor: Sorted taylor series.
!-

subroutine sort_universal_terms (ut_in, ut_sorted)

use nr, only: indexx
use definition, only: universal_taylor

implicit none

type (universal_taylor), target, intent(in)  :: ut_in
type (universal_taylor) :: ut_sorted

integer, allocatable :: ix(:), ord(:)
integer i, j, n, nv
integer, pointer :: expn(:)

! init

n = ut_in%n
nv = ut_in%nv

if (associated(ut_sorted%n)) deallocate(ut_sorted%n, ut_sorted%nv, ut_sorted%c, ut_sorted%j)
allocate(ut_sorted%n, ut_sorted%nv, ut_sorted%c(n), ut_sorted%j(n,nv), ix(n), ord(n))

ut_sorted%n = n
ut_sorted%nv = nv

!

do i = 1, n
  expn => ut_in%j(i,:)
  ord(i) = sum(expn)*10**nv
  do j = 1, nv-1
    ord(i) = ord(i) + expn(j) * 10**(j-1)
  enddo
enddo

call indexx (ord, ix)

do i = 1, n
  ut_sorted%c(i)= ut_in%c(ix(i))
  ut_sorted%j(i,:)= ut_in%j(ix(i),:)
enddo

deallocate(ord, ix)

end subroutine sort_universal_terms

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine type_map (y)
!
! Subroutine to type the transfer maps of a real_8 array.
!
! Input:
!   y(:)  -- Real_8: 
!-

subroutine type_map (y)

use s_fitting, only: assignment(=), real_8, universal_taylor

implicit none

type (real_8), intent(in) :: y(:)
type (universal_taylor) ut

integer :: i, j, k

!

do i = 1, size(y)
  ut = 0
  ut = y(i)%t
  print *, '!-----------------'
  print *, '! Term:', i
  print *, 'Order            Coef    Exponents'
  do j = 1, ut%n
    print '(i6, f18.14, 20i3)', sum(ut%j(j,:)), ut%c(j), (ut%j(j,k), k = 1, ut%nv)
  enddo
  ut = -1
enddo

end subroutine type_map

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+                                
! Subroutine ele_to_fibre (ele, ptc_fibre, param, use_offsets, integ_order, steps, 
!                                  for_layout, track_particle, use_hard_edge_drifts, kill_layout)
!
! Routine to convert a Bmad element to a PTC fibre element.
!
! Note: You need to call set_ptc before using this routine.
! Note: If ele contains a taylor_field, this routine may not be called in between calls to 
!   FPP alloc and kill since the setting up of the PTC pancake uses FPP.
!
! Input:
!   ele                  -- Ele_struct: Bmad element.
!   param                -- lat_param_struct: 
!   use_offsets          -- Logical: Does ptc_fibre include element offsets, pitches and tilt?
!   integ_order          -- Integer, optional: Order for the 
!                             sympletic integrator. Possibilities are: 2, 4, or 6
!                             Overrides ele%value(integrator_order$).
!                             default = 2 (if not set with set_ptc).
!   steps                -- Integer, optional: Number of integration steps.
!                             Overrides ele%value(ds_step$).
!   for_layout           -- Logical, optional: If True then fibre will be put in the PTC layout.
!                             Default is False.
!   track_particle       -- coord_struct, optional: Particle to be tracked. ref_particle$, electron$, etc.
!                             This argument should only be present when the fibre is not to be put in a layout.
!   use_hard_edge_drifts -- logical, optional: If False then hard edge drifts are not used.
!                                              If True then this argument has no effect.
!                              Default is set by bmad_com%use_hard_edge_drifts.
!                              Default bmad_com%use_hard_edge_drifts is True.
!   kill_layout          -- logical, optional: If False, and for_layout = False, do not kill the special layout
!                             containing any previous "no layout" fibres. Default is True.
!
! Output:
!   ptc_fibre -- Fibre: PTC fibre element.
!+

subroutine ele_to_fibre (ele, ptc_fibre, param, use_offsets, integ_order, steps, &
                                       for_layout, track_particle, use_hard_edge_drifts, kill_layout)

use madx_ptc_module

implicit none
 
type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct), optional :: track_particle
type (ele_struct), pointer :: field_ele, ele2
type (cartesian_map_term1_struct), pointer :: wt
type (cartesian_map_struct), pointer :: cm
type (cylindrical_map_struct), pointer :: cy
type (fibre), pointer :: ptc_fibre
type (keywords) ptc_key
type (work) energy_work
type (tree_element), pointer :: arbre(:)
type (c_damap) ptc_c_damap
type (real_8) ptc_re8(6)
type (taylor_field_struct), pointer :: tf
type (taylor_field_plane1_struct), pointer :: plane
type (em_taylor_term_struct), pointer :: tm
type(taylor), allocatable :: pancake_field(:,:)
type (taylor) ptc_taylor

real(rp) leng, hk, vk, s_rel, z_patch, phi_tot, norm, rel_charge, k1l, t1
real(rp) dx, dy, cos_t, sin_t, coef, coef_e, coef_b, kick_magnitude, ap_lim(2), ap_dxy(2), e1, e2
real(rp) beta0, beta1, ref0(6), ref1(6), fh, dz_offset
real(rp), pointer :: val(:)
real(rp), target :: value0(num_ele_attrib$)
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) ld, hd, lc, hc, angc, xc, dc

complex(rp) k_0

integer, optional :: integ_order, steps
integer i, ii, j, k, m, n, key, n_term, exception, ix, met, net, ap_type, ap_pos, ns, n_map
integer np, max_order, ix_pole_max, exact_model_saved, nn

logical use_offsets, kill_spin_fringe, onemap, found, is_planar_wiggler, use_taylor
logical, optional :: for_layout, use_hard_edge_drifts, kill_layout

character(16) :: r_name = 'ele_to_fibre'

!

val => ele%value
key = ele%key

if (.not. ele%is_on) then
  val => value0
  val = ele%value

  select case (ele%key)
  case (sbend$)
    val = 0
    val(l$)         = ele%value(l$)
    val(g$)         = ele%value(g$)
    val(g_err$)     = -ele%value(g$)
    val(angle$)     = ele%value(angle$)
    val(rho$)       = ele%value(rho$)
    val(ref_tilt_tot$) = val(ref_tilt_tot$)

  case (lcavity$)
    val(voltage$) = 0
    val(gradient$) = 0

  case (taylor$)
    key = drift$
    val(l$) = 0

  case default
    key = drift$
    val     = 0
    val(l$) = ele%value(l$) 
  end select
endif

!

n_map = 0
if (associated(ele%cylindrical_map)) n_map = n_map + 1
if (associated(ele%cartesian_map)) n_map = n_map + 1
if (associated(ele%taylor_field)) n_map = n_map + 1

if (n_map > 1) then
  call out_io (s_fatal$, r_name, 'PTC TRACKING IS ONLY ABLE TO HANDLE A SINGLE FIELD MAP IN AN ELEMENT.', &
                                 'ELEMENT HAS MULTIPLE FIELD MAPS: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

use_taylor = (n_map == 0 .and. (key == wiggler$ .or. key == undulator$) .and. ele%field_calc == helical_model$)
if (use_taylor) key = match$

! 

call zero_key(ptc_key)  ! init key

select case (ele%ptc_integration_type)
case (drift_kick$);  ptc_key%model = 'DRIFT_KICK'
case (matrix_kick$); ptc_key%model = 'MATRIX_KICK'
case (ripken_kick$); ptc_key%model = 'DELTA_MATRIX_KICK'
end select

if (key == sbend$ .and. ele%value(angle$) == 0 .and. ptc_key%model /= 'DRIFT_KICK') then
  ptc_key%model = 'DRIFT_KICK'
  ! Only need to issue a warning if K1 is nonzero.
  !if (ele%value(k1$) /= 0) call out_io (s_warn$, r_name, &
  !          'BEND WITH ZERO BENDING ANGLE WILL USE PTC_INTEGRATION_TYPE OF DRIFT_KICK: ' // ele%name)
endif

if (key == sbend$ .and. ele%value(g$) + ele%value(g_err$) == 0 .and. ptc_key%model /= 'DRIFT_KICK') then
  ptc_key%model = 'DRIFT_KICK'
  ! Only need to issue a warning if K1 is nonzero.
  !if (ele%value(k1$) /= 0) call out_io (s_warn$, r_name, &
  !          'BEND WITH ZERO NET BENDING FIELD WILL USE PTC_INTEGRATION_TYPE OF DRIFT_KICK: ' // ele%name)
endif

if (present(track_particle)) then
  rel_charge = charge_of(track_particle%species) / charge_of(param%particle)
else
  rel_charge = charge_of(default_tracking_species(param)) / charge_of(param%particle)
endif

leng = ele%value(l$)

if (use_offsets) then
  if (key == sbend$) then
    ptc_key%tiltd = ele%value(ref_tilt_tot$)
  else
    ptc_key%tiltd = ele%value(tilt_tot$)
  endif
else
  ptc_key%tiltd = 0
endif

ptc_key%method = nint(ele%value(integrator_order$))
if (ptc_key%method == 0) ptc_key%method = bmad_com%default_integ_order 
if (present(integ_order)) ptc_key%method = integ_order

if (present(steps)) then
  ptc_key%nstep = steps
elseif (leng == 0) then
  ptc_key%nstep = 1
elseif (key == taylor$ .or. key == match$ .or. key == multipole$ .or. &
                                key == ab_multipole$ .or. key == patch$) then
  ptc_key%nstep = 1
  leng = 0  ! Problem is that PTC will not ignore the length in tracking which is different from the Bmad convention.
else
  if (ele%value(ds_step$) == 0) then
    call out_io (s_fatal$, r_name, 'DS_STEP IS ZERO FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
  endif

  ptc_key%nstep = nint(abs(leng) / ele%value(ds_step$))
  if (ptc_key%nstep == 0) ptc_key%nstep = 1
endif

ptc_key%list%name = ele%name
ptc_key%list%l    = leng

! Fringes

ix = both_ends$
if (attribute_index(ele, 'FRINGE_AT') > 0) ix = nint(ele%value(fringe_at$))
kill_spin_fringe = is_false(ele%value(spin_fringe_on$))

ptc_key%list%kill_ent_fringe = (ix == exit_end$ .or. ix == no_end$)
ptc_key%list%kill_exi_fringe = (ix == entrance_end$ .or. ix == no_end$)
ptc_key%list%kill_ent_spin = (ix == exit_end$ .or. ix == no_end$ .or. kill_spin_fringe)
ptc_key%list%kill_exi_spin = (ix == entrance_end$ .or. ix == no_end$ .or. kill_spin_fringe)

!

if (key == sbend$ .and. ele%value(l$) == 0) key = kicker$
ele2 => pointer_to_field_ele(ele, 1, s_rel)
if (ele2%field_calc == fieldmap$ .and. ele2%tracking_method /= bmad_standard$) key = wiggler$

select case (key)

case (crab_cavity$)
  if (ele%value(rf_frequency$) == 0) then
    call out_io (s_fatal$, r_name, 'RF FREQUENCY IS ZERO FOR: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  ptc_key%magnet = 'rfcavity'
  ptc_key%list%n_bessel = 0
  !!ptc_key%list%volt = 1d-6 * e_accel_field(ele, voltage$)
  ptc_key%list%permfringe = 0
  ptc_key%list%cavity_totalpath = 0
  ptc_key%list%freq0 = ele%value(rf_frequency$)
  phi_tot = ele%value(phi0$) + ele%value(phi0_multipass$) + 0.25_rp 
  ptc_key%list%delta_e = 0     ! For radiation calc.
  ptc_key%list%lag = twopi * phi_tot

case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 
  if (ele%value(hkick$) == 0 .and. ele%value(vkick$) == 0) then
    ptc_key%magnet = 'drift'
  else
    ptc_key%magnet = 'quadrupole'
  endif

case (quadrupole$) 
  ptc_key%magnet = 'quadrupole'
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

case (sad_mult$)
  if (ele%value(l$) == 0) then
    ptc_key%magnet = 'multipole'  ! No Bz field if zero length.
  else
    ptc_key%magnet = 'solenoid'
    ptc_key%list%bsol = val(ks$)
    ! PTC tracking uses Matrix-Kick where the Matrix step uses the solenoid component and the Kick step
    ! uses the quadrupole (and higher order) components. SAD and Bmad combine the quadrupole component
    ! in the Matrix step. Thus PTC may need a smaller step size.
    if (.not. present(integ_order) .and. .not. present(steps)) then
      call multipole1_ab_to_kt (ele%a_pole(1), ele%b_pole(1), 1, k1l, t1)
      ptc_key%nstep = max(ptc_key%nstep, nint(ele%value(l$) * abs(k1l) / (ele%value(eps_step_scale$) * bmad_com%ptc_cut_factor)))
      ptc_key%method = 2
      if (ptc_key%nstep > 18) then
        ptc_key%nstep = nint(ptc_key%nstep / 7.0)
        ptc_key%method = 6
      elseif (ptc_key%nstep > 4) then
        ptc_key%nstep = nint(ptc_key%nstep / 3.0)
        ptc_key%method = 4
      endif
    endif
  endif

  if (ele%value(x_pitch_mult$) /= 0 .or. ele%value(y_pitch_mult$) /= 0) then
    call out_io (s_error$, r_name, &
          'NON-ZERO X OR Y_PITCH_MULT NOT IMPLEMENTED IN PTC FOR SAD_MULT ELEMENT: ' // ele%name)
  endif

case (sbend$) 
  ! PTC does not consider a finite e1/e2 part of the fringe so must zero e1/e2 if needed.
  ix = nint(ele%value(fringe_type$))

  e1 = ele%value(e1$)
  if (ptc_key%list%kill_ent_fringe .or. ix == none$) e1 = 0

  e2 = ele%value(e2$)
  if (ptc_key%list%kill_exi_fringe .or. ix == none$) e2 = 0

  ptc_key%magnet = 'sbend'
  ptc_key%list%b0   = ele%value(g$) * leng
  ptc_key%list%t1   = e1
  ptc_key%list%t2   = e2

  ptc_key%list%hgap = ele%value(hgap$)
  ptc_key%list%fint = ele%value(fint$)
  ptc_key%list%hgap2 = ele%value(hgapx$)
  ptc_key%list%fint2 = ele%value(fintx$)

  ix = nint(ele%value(ptc_field_geometry$))
  if (ix == straight$) then
    ptc_key%magnet = 'wedgrbend'
    ptc_key%list%t1   = e1 - ele%value(angle$)/2
    ptc_key%list%t2   = e2 - ele%value(angle$)/2
  elseif (ix == true_rbend$) then
    ptc_key%magnet = 'truerbend'
    ptc_key%list%t1   = e1 - ele%value(angle$)/2
    !! ptc_key%list%t2   = ele%value(e2$) - ele%value(angle$)/2 ! Determined by %t1 in this case.
  endif

case (sextupole$)
  ptc_key%magnet = 'sextupole'
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

case (taylor$, match$)
  ptc_key%magnet = 'marker'

case (octupole$)
  ptc_key%magnet = 'octupole'
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

case (solenoid$)
  ptc_key%magnet = 'solenoid'
  ptc_key%list%bsol = val(ks$)
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

case (sol_quad$)
  ptc_key%magnet = 'solenoid'
  ptc_key%list%bsol = val(ks$)
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

case (marker$, detector$, fork$, photon_fork$, beginning_ele$, patch$, floor_shift$, fiducial$)
  ptc_key%magnet = 'marker'
  ptc_key%nstep = 1

case (kicker$, hkicker$, vkicker$)
  ptc_key%magnet = 'kicker'

case (rfcavity$, lcavity$)
  if (ele%value(rf_frequency$) == 0) then
    call out_io (s_fatal$, r_name, 'RF FREQUENCY IS ZERO FOR: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  select case (nint(ele%value(cavity_type$)))
  case (traveling_wave$)
    ptc_key%magnet = 'twcavity'
    ptc_key%list%volt = 1d-6 * e_accel_field(ele, voltage$)
    ptc_key%list%cavity_totalpath = 1  ! 
  case (standing_wave$)
    ptc_key%magnet = 'rfcavity'
    ptc_key%list%cavity_totalpath = 1  ! 
    if (ptc_key%nstep == 1) ptc_key%nstep = 5  ! Avoid bug with nstep = 1.
    ! If an element has a length that is less than the pillbox length then avoid drifts with negative
    ! length and use a "fake" RF cavity model. 
    if (tracking_uses_end_drifts(ele)) then  ! Has end drifts means end drifts have positive length.
      ptc_key%list%n_bessel = -1   ! pillbox cavity.
      ptc_key%list%volt = 2d-6 * e_accel_field(ele, voltage$)
    else
      ptc_key%list%n_bessel = 1 ! Fake model
      ptc_key%list%volt = 1d-6 * e_accel_field(ele, voltage$)
    endif
  case (ptc_standard$)
    ptc_key%magnet = 'rfcavity'
    ptc_key%list%volt = 1d-6 * e_accel_field(ele, voltage$)
    ptc_key%list%n_bessel = 0
    ptc_key%list%permfringe = 0
    ptc_key%list%cavity_totalpath = 0
  end select

  ptc_key%list%freq0 = ele%value(rf_frequency$)
  phi_tot = ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + ele%value(phi0_autoscale$)
  if (tracking_uses_end_drifts(ele, use_hard_edge_drifts)) ptc_key%list%l = hard_edge_model_length(ele)

  if (key == lcavity$) then
    ptc_key%list%lag = pi / 2 - twopi * phi_tot
  else
    ptc_key%list%lag = twopi * phi_tot
  endif

  ptc_key%list%delta_e = 0     ! For radiation calc.

! ptc elsep cannot do spin tracking so use general electrostatic element instead.
case (elseparator$)
!  ptc_key%magnet = 'elseparator'
!  hk = val(hkick$) / leng
!  vk = val(vkick$) / leng
!  if (hk == 0 .and. vk == 0) then
!    ptc_key%tiltd = 0
!  else
!    if (param%particle < 0) then
!      hk = -hk
!      vk = -vk
!    endif
!    ptc_key%tiltd = -atan2 (hk, vk) + ele%value(tilt_tot$)
!  endif
!  ptc_key%list%volt = 1d-6 * ele%value(e_tot$) * sqrt(hk**2 + vk**2)

case (ab_multipole$, multipole$)
  ptc_key%magnet = 'multipole'

! beambeam element in PTC is a special drift that must be setup after the integration 
! node array of the fibre is created.

case (beambeam$)
  ptc_key%magnet = 'drift'

case (wiggler$, undulator$)
  ptc_key%magnet = 'wiggler'

case default
  call out_io (s_fatal$, r_name, 'CONVERSION TO PTC NOT IMPLEMENTED FOR ELEMENTS OF TYPE ' // trim(key_name(ele%key)), &
                                 'FOR ELEMENT: ' // trim(ele%name))
  if (global_com%exit_on_error) call err_exit
end select

! Fringe

if (key == sbend$ .and. ele%value(l$) /= 0) then

  ix = nint(ele%value(ptc_fringe_geometry$))
  ptc_key%list%bend_fringe = (ix == x_invariant$)

  ix = nint(ele%value(fringe_type$))
  select case (ix)
  case (none$)
    ptc_key%list%permfringe = 0
  case (basic_bend$, linear_edge$)
    ptc_key%list%permfringe = 0
  case (full$)
    ptc_key%list%permfringe = 1
  case (hard_edge_only$)
    ptc_key%list%permfringe = 1
  case (soft_edge_only$)
    ptc_key%list%permfringe = 2
  case (sad_full$)
    ptc_key%list%permfringe = 3
  end select

elseif (attribute_index(ele, 'FRINGE_TYPE') > 0) then  ! If fringe_type is a valid attribute

  ptc_key%list%bend_fringe = .false.

  ix = nint(ele%value(fringe_type$))
  select case (ix)
  case (none$)
    ptc_key%list%permfringe = 0
  case (hard_edge_only$)
    ptc_key%list%permfringe = 1
  case (soft_edge_only$)
    ptc_key%list%permfringe = 2
  case (full$)
    ptc_key%list%permfringe = 3
  end select

  if (key == sad_mult$ .and. ele%value(l$) == 0) ptc_key%list%permfringe = 0
endif

! Electric fields present? Everything is an sbend!

if (key /= multipole$ .and. (associated(ele%a_pole_elec) .or. key == elseparator$)) then
  ptc_key%magnet = 'sbend'
  ptc_key%model = 'DRIFT_KICK'   ! PTC demands this.
  ptc_key%exact = .true.  ! PTC does not implement a non-exact model when there are electric fields.
  SOLVE_ELECTRIC = .true.
  if (leng == 0) then
    call out_io (s_fatal$, r_name, 'ZERO LENGTH ELEMENT WITH AN ELECTRIC FIELD NOT ALLOWED IN PTC: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

! Magnetic multipole components

call ele_to_ptc_magnetic_an_bn (ele, param, ptc_key%list%k, ptc_key%list%ks, ptc_key%list%nmul)

! cylindrical field

if (associated(ele2%cylindrical_map) .and. ele2%field_calc == fieldmap$) then
  PUT_A_ABELL = 1
  ptc_key%magnet = 'abell_dragt'
  n_abell = 0
  do i = 1, size(ele%cylindrical_map)
    n_abell = max(2, n_abell, size(ele%cylindrical_map(i)%ptr%term))
  enddo
  m_abell = maxval(ele%cylindrical_map%m)
endif

! taylor_field

if (associated(ele2%taylor_field) .and. ele2%field_calc == fieldmap$) then
  tf => ele2%taylor_field(1)

  np = size(tf%ptr%plane)

  if (key == sbend$ .and. ele%value(g$) /= 0) then
    ld = ele%value(l$)
    hd = ele%value(g$)
    lc = (np-1) * tf%dz   ! Integration length

    if (tf%curved_ref_frame) then
      hc = ele%value(g$) ! pancake curvature = ele curvature
      angc = (ele%value(angle$) - (np-1) * tf%dz * ele%value(g$)) / 2
      xc = tf%r0(1) * cos(ele%value(angle$)/2) - ele%value(rho$) * (1 - cos(angc))
      dc = tf%r0(1) * sin(ele%value(angle$)/2) + ele%value(rho$) * sin(angc)
    else
      hc = 0.d0     ! pancake curvature
      angc = ele%value(angle$) / 2
      xc = tf%r0(1) + ele%value(rho$) * (1 - cos(ele%value(angle$)/2))
      dc = (ele%value(l$) - (np-1) * tf%dz) / 2
    endif

  else
    ld = ele%value(l$)
    lc = (np-1) * tf%dz   ! Integration length

    hd = 0
    hc = 0                ! pancake curvature

    angc = 0
    xc = tf%r0(1)
    dc = (ele%value(l$) - (np-1) * tf%dz) / 2
  endif

  if (mod(np, 2) == 0 .or. np < 5) then
    call out_io (s_fatal$, r_name, 'NUMBER OF PLANES FOR A TAYLOR_FIELD MUST BE ODD AND AT LEAST 5 IF USED WITH PTC.', &
                                   'TAYLOR_FIELD USED IN ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  if (tf%field_type /= magnetic$) then
    call out_io (s_fatal$, r_name, 'FIELD TYPE MUST BE MAGNETIC FOR A TAYLOR_FIELD IF USED WITH PTC.', &
                                   'TAYLOR_FIELD USED IN ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  n = len_trim(tf%ptr%file)
  call set_pancake_constants(np, angc, xc, dc, tf%r0(2), hc, lc, hd, ld, .not. tf%canonical_tracking, tf%ptr%file(max(1,n-23):n))

  max_order = 0
  do i = lbound(tf%ptr%plane, 1), ubound(tf%ptr%plane, 1)
    plane => tf%ptr%plane(i)
    do j = 1, 3
      do k = 1, size(plane%field(j)%term)
        max_order = max(max_order, sum(plane%field(j)%term(k)%expn))
      enddo
    enddo
  enddo

  call init (max_order+1, 1, 0, 0)
  call allocate_for_pancake (pancake_field)

  do i = lbound(tf%ptr%plane, 1), ubound(tf%ptr%plane, 1)
    ii = i + 1 - lbound(tf%ptr%plane, 1)
    plane => tf%ptr%plane(i)
    do j = 1, 3
      do k = 1, size(plane%field(j)%term)
        tm => plane%field(j)%term(k)
        coef = tm%coef * tf%field_scale * 1d9 * master_parameter_value(tf%master_parameter, ele2) ! * c_light / ele2%value(p0c$)
        pancake_field(j,ii)=pancake_field(j,ii)+(coef .mono. tm%expn)
      enddo
    enddo
  enddo

  ptc_key%magnet = 'INTERNALPANCAKE'

endif

!-----------------------------
! Create ptc_fibre
! The EXCEPTION argument is an error_flag. Set to 1 if there is an error. Never reset.

! The integration node array of the fibre is not setup if for_layout = True.
! [In this case the setup is done on the entire layout later.]
! Therefore only do the beambeam setup (which needs this array to exist) if for_layout = False.

call set_ptc_quiet(12, set$, n)

if (logic_option(.false., for_layout)) then
  call create_fibre_append (.true., m_u%end, ptc_key, EXCEPTION, br = pancake_field)   ! ptc routine
  ptc_fibre => m_u%end%end
  if (present (track_particle)) then
    call out_io (s_fatal$, r_name, 'TRACK_PARTICLE ARGUMENT SHOULD NOT BE PRESENT WHEN FOR_LAYOUT IS TRUE!')
    if (global_com%exit_on_error) call err_exit
    return
  endif

else
  call set_madx (energy = ele%value(e_tot$), method = ptc_key%method , step = ptc_key%nstep)

  if (logic_option(.true., kill_layout)) then
    if (associated(bmadl%start)) then
      call kill (bmadl)
      call set_up(bmadl)
    endif
    call create_fibre_append (.false., bmadl, ptc_key, EXCEPTION, br = pancake_field)   ! ptc routine
  else
    call create_fibre_append (.true., bmadl, ptc_key, EXCEPTION, br = pancake_field)   ! ptc routine
  endif

  ptc_fibre => bmadl%end

  ! NB: Set of pointers only needed if doing stuff other than tracking (like calculating misalignments).

  ptc_fibre%mag%p%dir     => ptc_fibre%dir
  ptc_fibre%mag%p%beta0   => ptc_fibre%beta0
  ptc_fibre%mag%p%gamma0i => ptc_fibre%gamma0i
  ptc_fibre%mag%p%gambet  => ptc_fibre%gambet
  ptc_fibre%mag%p%mass    => ptc_fibre%mass
  ptc_fibre%mag%p%charge  => ptc_fibre%charge

  ptc_fibre%magp%p%dir     => ptc_fibre%dir
  ptc_fibre%magp%p%beta0   => ptc_fibre%beta0
  ptc_fibre%magp%p%gamma0i => ptc_fibre%gamma0i
  ptc_fibre%magp%p%gambet  => ptc_fibre%gambet
  ptc_fibre%magp%p%mass    => ptc_fibre%mass
  ptc_fibre%magp%p%charge  => ptc_fibre%charge

  if (key == beambeam$) call beambeam_fibre_setup(ele, ptc_fibre, param)
endif

if (associated(ele2%taylor_field) .and. ele2%field_calc == fieldmap$) then
  call kill_for_pancake(pancake_field)
  call init_all (DEFAULT, ptc_com%taylor_order_saved, 0)
endif

ptc_fibre%dir = ele%orientation
if (present(track_particle)) ptc_fibre%dir = ptc_fibre%dir * track_particle%direction

call set_ptc_quiet(12, unset$, n)

! Cylindrical field

if (associated(ele2%cylindrical_map) .and. ele2%field_calc == fieldmap$) then
  do m = 0, m_abell
    found = .false.
    do j = 1, size(ele%cylindrical_map)
      cy => ele%cylindrical_map(j)
      if (m /= cy%m) cycle
      found = .true.
      exit
    enddo

    if (found) then
      if (cy%r0(1) /= 0 .or. cy%r0(2) /= 0) then
        call out_io (s_error$, r_name, 'TRANSVERSE R0 OFFSETS NOT YET IMPLEMENTED FOR PTC TRACKING FOR: ' // ele%name)
        if (global_com%exit_on_error) call err_exit
        return
      endif
      coef_e = cy%field_scale * master_parameter_value(cy%master_parameter, ele)
      if (key == lcavity$ .or. key == rfcavity$) coef_e = coef_e * ele%value(field_autoscale$)
      coef_b = coef_e * c_light / ele%value(p0c$)
      coef_e = -coef_e * 1d-6  ! Notice negative sign.

      ptc_fibre%mag%ab%t(m)  = cy%theta0_azimuth  ! Magnetic theta0
      ptc_fibre%mag%ab%te(m) = cy%theta0_azimuth  ! Electric theta0
      ptc_fibre%mag%ab%dz(m) = cy%dz

      n = size(cy%ptr%term%b_coef)
      select case (cy%ele_anchor_pt)
      case (anchor_beginning$)
        k_0 = -I_imaginary * twopi * cy%r0(3) / (n * cy%dz)
      case (anchor_center$)
        k_0 = -I_imaginary * twopi * (cy%r0(3) + ele%value(l$)/2) / (n * cy%dz)
      case (anchor_end$)
        k_0 = -I_imaginary * twopi * (cy%r0(3) + ele%value(l$)) / (n * cy%dz)
      end select

      nn = max(n, 2)  ! n = 1 is a singular case so treat it as n = 2
      do ii = 1, nn
        if (ii == 2 .and. n == 1) then
          k = ii - 1 - nn  ! = -1
          ptc_fibre%mag%ab%b(m,k) = 0
          ptc_fibre%mag%ab%e(m,k) = 0
        elseif (ii <= nn/2) then
          k = ii - 1 
          ptc_fibre%mag%ab%b(m,k) = coef_b * exp(k_0*k) * cy%ptr%term(ii)%b_coef
          ptc_fibre%mag%ab%e(m,k) = coef_e * exp(k_0*k) * cy%ptr%term(ii)%e_coef
        else
          k = ii - 1 - nn
          ptc_fibre%mag%ab%b(m,k) = coef_b * exp(k_0*k) * cy%ptr%term(ii)%b_coef
          ptc_fibre%mag%ab%e(m,k) = coef_e * exp(k_0*k) * cy%ptr%term(ii)%e_coef
        endif
      enddo
      !! ptc_fibre%mag%ab%b(m,:) = [cy%ptr%term(n/2+1:n)%b_coef, cy%ptr%term(1:n/2)%b_coef] * coef
    else
      ptc_fibre%mag%ab%t(m)   = 0
      ptc_fibre%mag%ab%dz(m)  = 1
      ptc_fibre%mag%ab%b(m,:) = 0
      ptc_fibre%mag%ab%e(m,:) = 0
    endif

    ptc_fibre%magp%ab%t(m)   = ptc_fibre%mag%ab%t(m)
    ptc_fibre%magp%ab%te(m)  = ptc_fibre%mag%ab%te(m)
    ptc_fibre%magp%ab%dz(m)  = ptc_fibre%mag%ab%dz(m)
    do ii = lbound(ptc_fibre%magp%ab%b, 2), ubound(ptc_fibre%magp%ab%b, 2)
      ptc_fibre%magp%ab%b(m,ii) = ptc_fibre%mag%ab%b(m,ii)
      ptc_fibre%magp%ab%e(m,ii) = ptc_fibre%mag%ab%e(m,ii)
    enddo

    ptc_fibre%mag%ab%xprime = xprime_abell   ! is_false(ele%value(ptc_canonical_coords$))
    ptc_fibre%magp%ab%xprime = ptc_fibre%mag%ab%xprime
  enddo

endif

!-----------------------------
! The E-field units that PTC wants on input are MV/m (MAD convention). 
! Note: PTC convert MV/m to GV/m internally.

if (key /= multipole$ .and. (associated(ele%a_pole_elec) .or. key == elseparator$)) then
  if (key /= sbend$) then
    ptc_fibre%mag%p%bend_fringe = .false.
    ptc_fibre%magp%p%bend_fringe = .false.
  endif
  fh = 1d-9 * sign_of(charge_of(param%particle)) / VOLT_C
  if (key == sbend$ .and. nint(ele%value(exact_multipoles$)) == vertically_pure$ .and. ele%value(g$) /= 0) then
    call multipole_ele_to_ab(ele, .false., ix_pole_max, a_pole, b_pole, electric$, include_kicks$)
    call convert_bend_exact_multipole(ele%value(g$), horizontally_pure$, a_pole, b_pole)
  else
    call multipole_ele_to_ab(ele, .false., ix_pole_max, a_pole, b_pole, electric$)
  endif

  if (key == elseparator$) then
    a_pole(0) = a_pole(0) + val(vkick$) * ele%value(p0c$) / leng
    b_pole(0) = b_pole(0) + val(hkick$) * ele%value(p0c$) / leng
  endif

  do i = 0, n_pole_maxx
    if (a_pole(i) /= 0) call add (ptc_fibre, -(i+1), 0, fh*a_pole(i), electric = .true.)
    if (b_pole(i) /= 0) call add (ptc_fibre,  (i+1), 0, fh*b_pole(i), electric = .true.)
  enddo
  SOLVE_ELECTRIC = .false.
endif

! Aperture

if (ele%aperture_at /= no_aperture$ .and. (ele%value(x1_limit$) /= 0 .or. ele%value(x2_limit$) /= 0) .and. &
            (ele%value(y1_limit$) /= 0 .or. ele%value(y2_limit$) /= 0) .and. &
            (ele%aperture_type == rectangular$ .or. ele%aperture_type == elliptical$)) then

  select case (ele%aperture_type)
  case (rectangular$);    ap_type = 2
  case (elliptical$);     ap_type = 1
  end select

  select case (ele%aperture_at)
  case (both_ends$);      ap_pos = 0
  case (entrance_end$);   ap_pos = -1
  case (exit_end$);       ap_pos = 1
  end select

  ap_lim = [(ele%value(x1_limit$) + ele%value(x2_limit$)) / 2, (ele%value(y1_limit$) + ele%value(y2_limit$)) / 2]
  ap_dxy = [(ele%value(x1_limit$) - ele%value(x2_limit$)) / 2, (ele%value(y1_limit$) - ele%value(y2_limit$)) / 2]
  call assign_aperture (ptc_fibre, ap_type, ap_lim, ap_lim(1), ap_lim(2), ap_dxy(1), ap_dxy(2), pos = ap_pos)
endif

! sad_mult & quadrupole
! Following the SAD convention: A zero length sad_mult has no fringe.

if ((key == sad_mult$ .and. ele%value(l$) /= 0) .or. key == quadrupole$) then
  ptc_fibre%mag%va  = -sign(sqrt(24 * abs(ele%value(fq1$))), ele%value(fq1$))
  ptc_fibre%magp%va = -sign(sqrt(24 * abs(ele%value(fq1$))), ele%value(fq1$))
  ptc_fibre%mag%vs  = ele%value(fq2$)
  ptc_fibre%magp%vs = ele%value(fq2$)
endif

if (key == sad_mult$) then
  cos_t = cos(ele%value(tilt_tot$))
  sin_t = sin(ele%value(tilt_tot$))
  dx =  ele%value(x_offset_mult$) * cos_t + ele%value(y_offset_mult$) * sin_t 
  dy = -ele%value(x_offset_mult$) * sin_t + ele%value(y_offset_mult$) * cos_t

  if (ptc_fibre%mag%kind == kind5) then
    ptc_fibre%mag%s5%dx = dx
    ptc_fibre%mag%s5%dy = dy

    ptc_fibre%magp%s5%dx = dx
    ptc_fibre%magp%s5%dy = dy
 elseif (ptc_fibre%mag%kind == kind3) then
    ptc_fibre%mag%k3%dx = dx
    ptc_fibre%mag%k3%dy = dy

    ptc_fibre%magp%k3%dx = dx
    ptc_fibre%magp%k3%dy = dy
  else
    call out_io (s_fatal$, r_name, 'INTERNAL ERROR SETTING MULT OFFSET. PLEASE CONTACT DAVID SAGAN.')
  endif
endif

! Set reference energy to the exit reference energy.

energy_work = 0
call find_energy (energy_work, p0c =  1d-9 * ele%value(p0c$))
ptc_fibre = energy_work

! FieldMap cartesian_map element. 
! Include all wiggler elements even planar_model with field_calc = bmad_standard$

if ((associated(ele2%cartesian_map) .and. ele2%field_calc == fieldmap$) .or. key == wiggler$ .or. key == undulator$) then

  is_planar_wiggler = ((key == wiggler$ .or. key == undulator$) .and. ele2%field_calc == planar_model$) 

  if (associated(ele2%grid_field)) then
    call out_io (s_fatal$, r_name, 'PTC TRACKING IS NOT ABLE TO USE GRID_FIELDS. FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  if (n_map == 0) then
    if (is_planar_wiggler) then
      allocate(cm)
      call create_wiggler_cartesian_map(ele2, cm)
    else
      call out_io (s_fatal$, r_name, 'NOT ABLE TO DO PTC TRACKING FOR NON-PLANAR WIGGLER WITHOUT A CARTESIAN (OR OTHER TYPE OF) MAP.', &
                                     'FOR ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
      return
    endif
  else
    cm => ele2%cartesian_map(1)
  endif

  if (cm%field_type == electric$) then
    call out_io (s_fatal$, r_name, 'PTC IS NOT ABLE TO HANDLE CARTESIAN_MAP WITH ELECTRIC FIELDS. FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif  

  select case (cm%ele_anchor_pt)
  case (anchor_beginning$); s_rel = s_rel - cm%r0(3)
  case (anchor_center$);    s_rel = s_rel - cm%r0(3) + ele%value(l$) / 2
  case (anchor_end$);       s_rel = s_rel - cm%r0(3) + ele%value(l$)
  end select

  n_term = size(cm%ptr%term)
  call POINTERS_W (ptc_fibre%mag%wi%w, n_term, 0)  ! n_term_electric needed   

  ptc_fibre%mag%wi%w%k(1,1:n_term)   = cm%ptr%term%kx
  ptc_fibre%mag%wi%w%k(2,1:n_term)   = cm%ptr%term%ky
  ptc_fibre%mag%wi%w%k(3,1:n_term)   = cm%ptr%term%kz
  ptc_fibre%mag%wi%w%f(1:n_term)     = cm%ptr%term%phi_z + s_rel * cm%ptr%term%kz
  ptc_fibre%mag%wi%w%x0(1:n_term)    = cm%ptr%term%x0
  ptc_fibre%mag%wi%w%y0(1:n_term)    = cm%ptr%term%y0
  ptc_fibre%mag%wi%w%form(1:n_term)  = 3*(cm%ptr%term%family - 1) + cm%ptr%term%form

  if (ele%is_on) then
    do i = 1, size(ptc_fibre%mag%wi%w%a(1:n_term))
      wt => cm%ptr%term(i)
      ptc_fibre%mag%wi%w%a(i) = c_light * ele2%value(polarity$) * wt%coef / ele%value(p0c$)
    enddo
  else
    ptc_fibre%mag%wi%w%a(1:n_term) = 0
  endif

  ! Correct z-position 

  z_patch = ele%value(delta_ref_time$) * c_light * ele%value(p0c$) / ele%value(e_tot$) - ele%value(l$)
  ptc_fibre%mag%wi%internal(6) = z_patch

  call copy (ptc_fibre%mag, ptc_fibre%magp)

  if (n_map == 0 .and. is_planar_wiggler) then
    deallocate(cm%ptr)
    deallocate(cm)
  endif
endif

! Misalignments and patches...

call misalign_ele_to_fibre (ele, use_offsets, ptc_fibre)

! Set charge

ptc_fibre%charge = rel_charge

! Taylor maps
! In theory can put in a taylor map for any element but for now only setup Bmad taylor and match elements.

if (key == taylor$ .or. key == match$) then
  ! The map can be split into pieces by taking the log of the map.
  ! onemap = T means do not split.
  ! At this point in time there is no splitting allowed.

  if (.not. associated(ele%taylor(1)%term)) call ele_to_taylor (ele, param)

  onemap = .true.

  call alloc(ptc_c_damap)

  ! 

  beta0 = ele%value(p0c_start$)/ele%value(e_tot_start$)
  beta1 = ele%value(p0c$)/ele%value(e_tot$)  

  call alloc (ptc_re8)
  call alloc (ptc_taylor)

  call taylor_to_real_8 (ele%taylor, beta0, beta1, ptc_re8, ref0, ref1)
  ptc_c_damap = ptc_re8

  if (associated(ele%spin_taylor(1)%term)) then
    do j = 0, 3
      ptc_taylor = ele%spin_taylor(j)
      ptc_c_damap%q%x(j) = ptc_taylor
    enddo
  endif

  call kill (ptc_taylor)
  call kill (ptc_re8);

  ! 

  if (.not. onemap) then
    call nth_root(ptc_c_damap, ptc_c_damap, ptc_fibre%mag%p%nst)
  endif

  if (ptc_fibre%dir == 1) then
    if (.not.associated(ptc_fibre%mag%forward)) then 
      allocate(ptc_fibre%mag%forward(3))
    else
      call KILL(ptc_fibre%mag%forward)  ! The kill only zeros the Berz-part, it stays associated.
    endif

    call SET_TREE_G_complex(ptc_fibre%mag%forward, ptc_c_damap)
    ptc_fibre%mag%do1mapf = onemap
    ptc_fibre%mag%usef = .true.
    arbre => ptc_fibre%mag%forward
  else
    if (.not. associated(ptc_fibre%mag%backward)) then 
      allocate(ptc_fibre%mag%backward(3))
    else
      call KILL(ptc_fibre%mag%backward)  ! The kill only zeros the Berz-part, it stays associated.
    endif
    call SET_TREE_G_complex(ptc_fibre%mag%backward, ptc_c_damap)
    ptc_fibre%mag%do1mapb = onemap
    ptc_fibre%mag%useb = .true.
    arbre => ptc_fibre%mag%backward
  endif

  !! arbre => ptc_fibre%mag%forward
  call mat_make_unit(arbre(1)%rad)    ! Radiation damping matrix. Unit matrix  => radiation off
  arbre(1)%fix0(1:6) = ref0
  arbre(1)%fixr(1:6) = ref1
  arbre(1)%fix(1:6) = ref1

  if (onemap) then
    arbre(1)%ds = ptc_fibre%mag%p%ld
  else
    arbre(1)%ds = ptc_fibre%mag%p%ld/ptc_fibre%mag%p%nst 
  endif

  arbre(1)%beta0 = ptc_fibre%beta0

  if (ptc_fibre%dir == 1) then
    if (.not.associated(ptc_fibre%magp%forward)) then 
      allocate(ptc_fibre%magp%forward(3))
      allocate(ptc_fibre%magp%usef)
    else
      call KILL(ptc_fibre%magp%forward)
    endif
    !call SET_TREE_G_complex(ptc_fibre%magp%forward,m)
    do i = 1, 3
      call alloc_tree(ptc_fibre%magp%forward(i), ptc_fibre%mag%forward(i)%n, ptc_fibre%mag%forward(i)%np)
      call copy_tree(ptc_fibre%mag%forward(i), ptc_fibre%magp%forward(i))
    enddo
    ptc_fibre%magp%do1mapf = onemap
    ptc_fibre%magp%usef = .true.             ! Use the Taylor map by default 
    arbre => ptc_fibre%magp%forward

  else
    if (.not.associated(ptc_fibre%magp%backward)) then 
      allocate(ptc_fibre%magp%backward(3))
      allocate(ptc_fibre%magp%useb)
    else
      call KILL(ptc_fibre%magp%backward)
    endif
    ! call SET_TREE_G_complex(ptc_fibre%magp%backward,m)
    do i = 1, 3
      call alloc_tree(ptc_fibre%magp%backward(i), ptc_fibre%mag%backward(i)%n, ptc_fibre%mag%backward(i)%np)
      call copy_tree(ptc_fibre%mag%backward(i), ptc_fibre%magp%backward(i))
    enddo
    ptc_fibre%magp%do1mapb = onemap
    ptc_fibre%magp%useb = .true.
    arbre => ptc_fibre%magp%backward
  endif

  ! File names in case a flat file is to be created.
  write (ptc_fibre%mag%filef, '(2a, i0, a)') trim(downcase(ele%name)), '_', ele%ix_ele, '.txf'
  !!write (ptc_fibre%mag%fileb, '(2a, i0, a)') trim(downcase(ele%name)), '_', ele%ix_ele, '.txb'

  !

  call mat_make_unit(arbre(1)%rad)    ! Radiation damping matrix. Unit matrix  => radiation off
  arbre(1)%fix0(1:6) = ref0
  arbre(1)%fixr(1:6) = ref1
  arbre(1)%fix(1:6) = ref1

  if (onemap) then
    arbre(1)%ds = ptc_fibre%mag%p%ld    ! Length of magnet (arc of circle if a bend)
  else
    arbre(1)%ds = ptc_fibre%mag%p%ld/ptc_fibre%mag%p%nst
  endif

  arbre(1)%beta0 = ptc_fibre%beta0
   
  call kill(ptc_c_damap)
endif

! Customization if wanted

call ele_to_fibre_hook (ele, ptc_fibre, param)

end subroutine ele_to_fibre

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine beambeam_fibre_setup(ele, ptc_fibre, param)
!
! Routine to setup a fibre to handle the beambeam interaction.
!
! Input:
!   ele       -- ele_struct: Bmad beambeam element.
!   param     -- lat_param_struct: Lattice parameters.
!
! Output:
!    ptc_fibre  -- Corresponding PTC fibre.
!-

subroutine beambeam_fibre_setup (ele, ptc_fibre, param)

use madx_ptc_module, only: integration_node, alloc, DO_BEAM_BEAM

implicit none

type (ele_struct), target :: ele
type (fibre), target :: ptc_fibre
type (lat_param_struct) param
type (integration_node), pointer :: node

real(rp) sig_x0, sig_y0, beta_a0, beta_b0, alpha_a0, alpha_b0, sig_x, sig_y
real(rp) alpha, beta, s_pos, bbi_const
real(rp), allocatable :: z_slice(:)

integer i, n_slice

character(*), parameter :: r_name = 'beambeam_fibre_setup'

!

DO_BEAM_BEAM = .true.

sig_x0 = ele%value(sig_x$)
sig_y0 = ele%value(sig_y$)
if (sig_x0 == 0 .or. sig_y0 == 0) then
  call out_io (s_error$, r_name, 'STRONG BEAM SIGMAS NOT SET FOR BEAMBEAM ELEMENT: ' // ele%name, &
                                 'SIGMAS WILL BE SET TO 1 METER.')
  sig_x0 = 1
  sig_y0 = 1
endif

beta_a0 = 0;  beta_b0 = 0

if (ele%value(beta_a_strong$) == 0) then
  beta_a0 = ele%a%beta
  alpha_a0 = ele%a%alpha
else
  beta_a0 = ele%value(beta_a_strong$)
  alpha_a0 = ele%value(alpha_a_strong$)
endif

if (ele%value(beta_b_strong$) == 0) then
  beta_b0 = ele%b%beta
  alpha_b0 = ele%b%alpha
else
  beta_b0 = ele%value(beta_b_strong$)
  alpha_b0 = ele%value(alpha_b_strong$)
endif

n_slice = max(1, nint(ele%value(n_slice$)))
allocate (z_slice(n_slice))
call bbi_slice_calc (ele, n_slice, z_slice)

! The fibre has already been initalized to be a drift by ele_to_fibre.
! The beambeam information is put in the third integration node.

node => ptc_fibre%t1%next%next
if (.not. associated(node%bb)) call alloc(node%bb, n_slice, 0.0_dp)

node%bb%n = n_slice

do i = 1, n_slice
  s_pos = z_slice(i) / 2  ! Factor of 2 since strong beam is moving.
  if (beta_a0 == 0 .or. beta_b0 == 0) then
    sig_x = sig_x0
    sig_y = sig_y0
  else
    beta = beta_a0 - 2 * alpha_a0 * s_pos + (1 + alpha_a0**2) * s_pos**2 / beta_a0
    sig_x = sig_x0 * sqrt(beta / beta_a0)
    beta = beta_b0 - 2 * alpha_b0 * s_pos + (1 + alpha_b0**2) * s_pos**2 / beta_b0
    sig_y = sig_y0 * sqrt(beta / beta_b0)
  endif

  bbi_const = -2.0_rp * param%n_part * ele%value(charge$) * classical_radius_factor / ele%value(e_tot$)

  node%bb%bbk(i,:) = 0  ! MAD closed orbit kick. Not used here.
  node%bb%xm(i) = (ele%value(crab_x1$) * s_pos + ele%value(crab_x2$) * s_pos**2 + &
                                                 ele%value(crab_x3$) * s_pos**3) * cos(ele%value(crab_tilt$))
  node%bb%ym(i) = (ele%value(crab_x1$) * s_pos + ele%value(crab_x2$) * s_pos**2 + &
                                                 ele%value(crab_x3$) * s_pos**3) * sin(ele%value(crab_tilt$))
  ! Transverse displacement.
  node%bb%sx(i) = sig_x
  node%bb%sy(i) = sig_y
  node%bb%fk(i) = -bbi_const / n_slice
  node%bb%s(i)  = z_slice(i) / 2   ! Factor of 2 due to strong beam velocity.
enddo

end subroutine beambeam_fibre_setup

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine misalign_ele_to_fibre (ele, use_offsets, ptc_fibre)
!
! Routine to misalign a fibre associated with a Bmad element.
!
! Input:
!   ele -- ele_struct: Bmad element with misalignments.
!   use_offsets -- Logical: Does ptc_fibre include element offsets, pitches and tilt?
!                   This argument is ignored if the element is a patch or floor_shit element.
!
! Output:
!   ptc_fibre -- Fibre: PTC fibre element with misalignments.
!-

subroutine misalign_ele_to_fibre (ele, use_offsets, ptc_fibre)

use madx_ptc_module

implicit none

type (ele_struct) ele
type (floor_position_struct) :: floor0, floor1
type (fibre), pointer :: ptc_fibre
type (fibre) dummy_fibre

real(rp) dr(3), ang(3), exi(3,3)
real(dp) mis_rot(6), beta_start, beta_end
real(dp) omega(3), basis(3,3), angle(3), tiltd
real(rp) x_off, y_off, x_pitch, y_pitch, roll, z_off

logical use_offsets, good_patch

character(*), parameter :: r_name = 'misalign_ele_to_fibre'

! In ptc there is no such thing as a reversed patch. Therefore need to
! use the reverse transformation if the patch is reversed in Bmad.

! Also fibre%dir for a patch must agree with the preceeding element in a layout

if (ele%key == patch$ .or. ele%key == floor_shift$) then
  if (ele%orientation == -1) then
    call ele_geometry (floor0, ele, floor1)
    dr = floor1%r
    ang = [floor1%theta, floor1%phi, floor1%psi]
  else
    dr = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)]
    ang = [ele%value(y_pitch$), ele%value(x_pitch$), ele%value(tilt$)]
  endif

  call bmad_patch_parameters_to_ptc (ang, exi)

  call set_madx_(.true., .false.)
  dummy_fibre = marker('dummy')
  call set_madx_(.false., .false.)

  ptc_fibre%dir = nint(ele%value(upstream_ele_dir$))
  dummy_fibre%dir = nint(ele%value(downstream_ele_dir$))

  call survey (dummy_fibre, exi, dr)
  call find_patch (ptc_fibre, dummy_fibre, next = .false., patching=good_patch)
  if (.not. good_patch) then
    call out_io (s_fatal$, r_name, 'CANNOT COMPUTE PTC PATCH FOR: ' // ele%name)
    return
  endif

  call super_zero_fibre(dummy_fibre, -1)

  ! energy and time patches

  if (ele%key == patch$) then

    beta_start = ele%value(p0c_start$) / ele%value(E_tot_start$)
    beta_end = ele%value(p0c$) / ele%value(E_tot$)

    if (beta_start == beta_end) then
      ptc_fibre%patch%energy = 0
    else
      ptc_fibre%patch%energy = 4 ! Internal entrance patch. Must be done after find_patch call.
      ptc_fibre%patch%p0b = ele%value(p0c_start$) * 1d-9
      ptc_fibre%patch%b0b = beta_start  ! beta velocity
    endif

    ! PTC uses beta_start for the reference time calculation while Bmad uses beta_end so
    ! renormalize the patch length to get PTC to agree with Bmad.

    ptc_fibre%patch%time = 2     ! Subtract off reference time (which affects z in tracking).
    ptc_fibre%patch%b_t = ele%value(l$) / beta_end + ele%value(t_offset$) * c_light 
    ptc_fibre%patch%b_l = ele%value(l$) + ele%value(t_offset$) * c_light * beta_end
  endif

!----------------------------------------------------------------------
! Not patch nor floor_shift element.

elseif (use_offsets) then

  ! Patch elements do not have misalignments

  if (ele%key == fiducial$) return
  if (attribute_index(ele, 'X_OFFSET_TOT') < 1) return

  ! in PTC the reference point for the offsets is the beginning of the element.
  ! In Bmad the reference point is the center of the element..

  x_off = ele%value(x_offset_tot$)
  y_off = ele%value(y_offset_tot$)
  z_off = ele%value(z_offset_tot$)
  x_pitch = ele%value(x_pitch_tot$)
  y_pitch = ele%value(y_pitch_tot$)
  tiltd = ptc_fibre%mag%p%tiltd
  roll = 0
  if (ele%key == sbend$) roll = ele%value(roll_tot$)

  if (x_off /= 0 .or. y_off /= 0 .or. z_off /= 0 .or. x_pitch /= 0 .or. &
                              y_pitch /= 0 .or. roll /= 0 .or. tiltd /= 0) then
    mis_rot = [x_off, y_off, z_off, -y_pitch, -x_pitch, roll]
    angle = 0
    angle(3) = -tiltd
    omega = ptc_fibre%chart%f%o
    basis = ptc_fibre%chart%f%mid
    call geo_rot(basis, angle, 1, basis)                     ! PTC call
    call misalign_fibre (ptc_fibre, mis_rot, omega, basis)   ! PTC call
  endif

  if (ele%value(e_tot_start$) /= ele%value(e_tot$)) then
    ptc_fibre%patch%energy = 4   ! Entrance energy patch
    ptc_fibre%patch%p0b = ele%value(p0c_start$) * 1d-9
    ptc_fibre%patch%b0b = ele%value(p0c_start$) / ele%value(E_tot_start$)  ! beta velocity
  endif

endif

end subroutine misalign_ele_to_fibre

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
! converts patch parameters from bmad to ptc

subroutine bmad_patch_parameters_to_ptc(ang, exi)

use madx_ptc_module

implicit none

real(dp) ent(3,3), exi(3,3)
real(dp) a(3), ang(3)

!

ent=global_frame
exi=ent

a=0.0_dp
a(3)=ang(3)
call geo_rot(ent, a, 1, basis=exi)
!exi=ent

a=0.0_dp
a(1)=-ang(1)
call geo_rot(ent, a, 1, basis=exi)
!exi=ent

a=0.0_dp
a(2)=-ang(2)
call geo_rot(ent, a, 1, basis=exi)
exi=ent

end subroutine bmad_patch_parameters_to_ptc

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+                                
! Subroutine ele_to_ptc_magnetic_an_bn (ele, param, bn, an, n_max)
!
! Routine to compute the a(n) and b(n) magnetic multipole components of a magnet.
! This is used to interface between eles and PTC fibres
!
! Note: On ptc side bn(1) is error field when creating a fibre but 
! is total field when fibre is being modified. This routine returns the error field.
!
! Input:
!   ele                 -- ele_struct: Bmad Element.
!   param               -- lat_param_struct:
!
! Output:
!   bn(1:n_pole_maxx+1) -- real(rp): Normal multipole component.
!   an(1:n_pole_maxx+1) -- real(rp): Skew multipole component.
!   n_max               -- integer, optional: Maximum non-zero multipole component.
!-

subroutine ele_to_ptc_magnetic_an_bn (ele, param, bn, an, n_max)

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param

real(rp) bn(:), an(:)
real(rp) cos_t, sin_t, leng, hk, vk, tilt
real(rp), pointer :: val(:)
real(rp), target, save :: value0(num_ele_attrib$) = 0
real(rp) an0(0:n_pole_maxx), bn0(0:n_pole_maxx)

integer, optional :: n_max
integer n, key, ix_pole_max
logical kick_here, add_kick, add_multipoles

character(*), parameter :: r_name = 'ele_to_ptc_magnetic_an_bn'

!

if (ele%key == taylor$) return
if (ele%key == match$) return

leng = ele%value(l$)

key = ele%key
if (ele%is_on) then
  val => ele%value
else
  val => value0  ! Not is_on -> has zero strength.
endif

bn = 0
an = 0
if (present(n_max)) n_max = 0
add_kick = .true.
add_multipoles = .true.

select case (key)

case (marker$, detector$, fork$, photon_fork$, beginning_ele$, em_field$, patch$, fiducial$, floor_shift$)
  return

case (crab_cavity$)
  if (leng == 0) then
    bn(1) = 1d-9 * e_accel_field(ele, voltage$)
  else
    bn(1) = 1d-9 * e_accel_field(ele, voltage$) / leng
  endif

case (drift$, rfcavity$, lcavity$, ab_multipole$, multipole$, beambeam$, wiggler$, undulator$)
  ! Nothing to be done

case (octupole$)
  bn(4) = val(k3$) / 6

case (quadrupole$) 
  bn(2) = val(k1$)

case (sad_mult$)

case (sbend$)
  if (ele%is_on) then
    bn(1) = ele%value(g_err$)
  else
    bn(1) = -ele%value(g$)
  endif

  ! PTC assumes horizontally_pure
  if (nint(ele%value(exact_multipoles$)) == vertically_pure$ .and. ele%value(g$) /= 0) then
    call multipole_ele_to_ab (ele, .false., ix_pole_max, an, bn, magnetic$, include_kicks$)
    call convert_bend_exact_multipole(ele%value(g$), horizontally_pure$, an, bn)
    if (leng /= 0) then
      an = an / leng
      bn = bn / leng
    endif
    add_kick = .false.   
    add_multipoles = .false.

  else
    bn(2) = val(k1$)
    bn(3) = val(k2$) / 2
  endif

case (sextupole$)
  bn(3) = val(k2$) / 2

case (rcollimator$, ecollimator$, monitor$, instrument$, pipe$)

case (solenoid$)

case (sol_quad$)
  bn(2) = val(k1$)

case (hkicker$, vkicker$)
  if (ele%key == hkicker$) bn(1) = val(kick$) 
  if (ele%key == vkicker$) an(1) = val(kick$) 
  add_kick = .false.

case (kicker$)
  bn(1)  = val(hkick$) 
  an(1) = val(vkick$) 
  add_kick = .false.

case (elseparator$)
  call multipole_ele_to_ab (ele, .false., ix_pole_max, an0, bn0) 
  if (ix_pole_max > -1) then
    call out_io (s_fatal$, r_name, 'MULTIPOLES IN AN ELSEPARATOR NOT SUPPORTED IN A FIBRE.')
    if (global_com%exit_on_error) call err_exit
    an0 = 0; bn0 = 0
  endif
  return

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN ELEMENT CLASS: ' // key_name(ele%key), &
                                 'FOR ELEMENT: ' // trim(ele%name))
  if (global_com%exit_on_error) call err_exit

end select

if (add_kick .and. has_hkick_attributes(ele%key) .and. (val(hkick$) /= 0 .or. val(vkick$) /= 0)) then
  hk = val(hkick$) / leng   ! PTC uses scaled kick for non-kicker elements.
  vk = val(vkick$) / leng
  if (ele%key == sbend$) then
    tilt = ele%value(ref_tilt_tot$)
  else
    tilt = ele%value(tilt_tot$)
  endif
  cos_t = cos(tilt)
  sin_t = sin(tilt)
  bn(1) = bn(1) - hk * cos_t - vk * sin_t
  an(1) = an(1) - hk * sin_t + vk * cos_t
endif

! Bmad an and bn are integrated fields. PTC uses just the field.
! Exception is multipole element.

if (associated(ele%a_pole) .and. add_multipoles) then
  call multipole_ele_to_ab (ele, .false., ix_pole_max, an0, bn0)

  if (leng /= 0) then
    an0 = an0 / leng
    bn0 = bn0 / leng
  endif

  n = n_pole_maxx+1
  an(1:n) = an(1:n) + an0(0:n-1)
  bn(1:n) = bn(1:n) + bn0(0:n-1)
endif

if (present(n_max)) then
  do n = size(bn), 1, -1
    if (an(n) /= 0 .or. bn(n) /= 0) exit
  enddo
  n_max  = n
endif

end subroutine ele_to_ptc_magnetic_an_bn

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine apply_patch_to_ptc_fibre (ele)
!
! Routine to take the patch parameters from a Bmad patch element and
! transfer them to the associated PTC fibre.
!
! Input:
!   ele           -- ele_struct: Patch element.
!
! Output:
!   ele%ptc_fibre -- PTC Fibre which should be a marker.
!-

subroutine apply_patch_to_ptc_fibre (ele)

use s_family

implicit none

type (ele_struct) ele

real(dp) mis_rot(6), dr(3)
real(dp) omega(3), basis(3,3), angle(3), origin(3), frame(3,3)

!

dr = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)]
origin = 0
omega = dr + origin


frame = global_frame
basis = global_frame

angle = [0.0d0, 0.0d0, ele%value(tilt$)]
call geo_rot(basis, angle, 1, basis=frame)     ! PTC call
frame = basis


angle = [0.0d0, ele%value(y_pitch$), 0.0d0]
call geo_rot(basis, angle, 1, basis=frame)     ! PTC call
frame = basis

angle = [ele%value(x_pitch$), 0.0d0, 0.0d0]
call geo_rot(basis, angle, 1, basis=frame)     ! PTC call
frame = basis

basis = global_frame
call find_patch (origin, basis, omega, frame, dr, angle)

ele%ptc_fibre%patch%patch = 2    ! Means entrance patch is not used exit patch is used.
ele%ptc_fibre%patch%b_d = dr
ele%ptc_fibre%patch%b_ang = angle

end subroutine apply_patch_to_ptc_fibre

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine set_ptc_quiet (channel, set, old_val)
!
! Routine to set the lielib_print(:) array or c_verbose logical to suppress informational messages
! that can clutter the output from a program using PTC.
!
! Note: Only suppress printing if bmad_com%ptc_print_info_messages = F.
!
! Input:
!   channel     -- integer: Index in the lielib_print(:) array to set. 0 => c_verbose.
!   set         -- logical: If set$ then set lielib_print(:). If unset$ then undo a previous set$.
!   old_val     -- integer: Old value needed for set = unset$.
!
! Output:
!   old_val     -- integer: Saved value for set = set$.
!-

subroutine set_ptc_quiet (channel, set, old_val)

use c_TPSA, only: c_verbose, lielib_print

implicit none

integer channel, old_val
logical set

!

if (bmad_com%ptc_print_info_messages) return

if (set .eqv. set$) then
  if (channel == 0) then
    old_val = int_logic(c_verbose)
    c_verbose = .false.
  else
    old_val = lielib_print(channel)
    lielib_print(channel) = 0
  endif

!

else
  if (channel == 0) then
    c_verbose = (old_val == 1)
  else
    lielib_print(channel) = old_val
  endif
endif

end subroutine set_ptc_quiet

end module
