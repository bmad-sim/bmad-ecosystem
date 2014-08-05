!+
! Module ptc_interface_mod
!
! Module of basic PTC interface routines.
! Also see: ptc_layout_mod
!-

module ptc_interface_mod

use multipole_mod
use bookkeeper_mod
use bmad_complex_taylor_mod

interface assignment (=)
  module procedure real_8_equal_taylor
  module procedure taylor_equal_real_8
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

use polymorphic_taylor, only: kill, operator(+), real_8

implicit none

type (taylor_struct), intent(in) :: taylor1(:), taylor2(:)
type (taylor_struct) taylor3(size(taylor1))
type (real_8) y1(size(taylor1)), y2(size(taylor1)), y3(size(taylor1))

integer i

! set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

!

call real_8_init(y1)
call real_8_init(y2)
call real_8_init(y3)

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

use polymorphic_taylor, only: kill, operator(-), real_8

implicit none

type (taylor_struct), intent(in) :: taylor1(:), taylor2(:)
type (taylor_struct) taylor3(size(taylor1))
type (real_8) y1(size(taylor1)), y2(size(taylor1)), y3(size(taylor1))

integer i

! set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

!

call real_8_init(y1)
call real_8_init(y2)
call real_8_init(y3)

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
! Modules Needed:                    
!   use ptc_interface_mod
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
! Modules needed:
!   use ptc_interface_mod
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
! Module Needed:
!   use ptc_interface_mod
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
! Module Needed:
!   use ptc_interface_mod
!
! Input:
!   ptc_fibre    -- fibre, pointer: Fibre to type info of.
!   print_coords -- logical, optional: If True then print coordinate and  patch information. 
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
nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%pos    [Index in layout]:   ', int_of(ptc_fibre%pos, 'i0')
nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%loc    [Index in universe]: ', int_of(ptc_fibre%loc, 'i0')
nl=nl+1; write (li(nl), '(2x, a, t40, a)') '%dir    [Direction]:         ', int_of(ptc_fibre%dir, 'i0')
nl=nl+1; write (li(nl), '(2x, a, t39, a)') '%beta0  [Beta velocity]:     ', real_of(ptc_fibre%beta0, 'es13.5')
nl=nl+1; write (li(nl), '(2x, a, t39, a)') '%mass   [Mass]:              ', real_of(ptc_fibre%mass, 'es13.5', 1e9_rp)
nl=nl+1; write (li(nl), '(2x, a, t39, a)') '%charge [Charge]:            ', real_of(ptc_fibre%charge, 'es13.5')

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

  call write_real ('%l      [Len]:    ',  mag%volt,    'es13.5', writeit = .true.)
  call write_real ('%volt             ',  mag%volt,    'es13.5')
  call write_real ('%freq             ',  mag%freq,    'es13.5')
  call write_real ('%phas             ',  mag%phas,    'es13.5')
  call write_real ('%volt             ',  mag%volt,    'es13.5')
  call write_real ('%delta_e          ',  mag%delta_e, 'es13.5')
  call write_real ('%fint             ',  mag%fint,    'es13.5')
  call write_real ('%hgap             ',  mag%hgap,    'es13.5')
  call write_real ('%h1               ',  mag%h1,      'es13.5')
  call write_real ('%h2               ',  mag%h2,      'es13.5')
  call write_real ('%b_sol            ',  mag%b_sol,   'es13.5')
  call write_real ('%va [sad f1]      ',  mag%va,      'es13.5')
  call write_real ('%vs [sad f2]      ',  mag%va,      'es13.5')

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
      call write_reals ('%a_d    [Entrance translation]:', ptch%a_d,   '10f12.9', writeit = .true.)
      call write_reals ('%a_ang  [Entrance angle]:      ', ptch%a_ang, '10f12.9', writeit = .true.)
    endif

    printit = .false.
    if (integer2_value(0, ptch%patch) == 2 .or. integer2_value(0, ptch%patch) == 3) printit = .true.

    if (printit) then
      call write_reals ('%b_d    [Exit translation]:', ptch%b_d,   '10f12.9', writeit = .true.)
      call write_reals ('%b_ang  [Exit angle]:      ', ptch%b_ang, '10f12.9', writeit = .true.)
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
!   kind_str -- Character(40): String representation
!-

function kind_name (this_kind) result (kind_str)

use s_status, only: kind0, kind1, kind2, kind3, kind4, kind5, kind6, kind7, kind8, kind9, &
                    kind10, kind11, kind12, kind13, kind14, kind15, kind16, kind17, kind18, &
                    kind19, kind20, kind21, kind22, kind23, kindwiggler, kindpa

implicit none

integer, pointer :: this_kind
character(40) kind_str

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
! Subroutine set_ptc (e_tot, particle, taylor_order, integ_order, &
!                               n_step, no_cavity, exact_modeling, exact_misalign)
!
! Subroutine to initialize PTC.
! This subroutine uses the FPP/PTC routines:
!     make_states
!     set_mad
!     init
!
! Note: At some point before you use PTC to compute Taylor maps etc.
!   you have to call set_ptc with both e_tot and particle args present. 
!   Always supply both of these args together or not at all. 
! Note: If you just want to use FPP without PTC then call the FPP routine init directly.
! Note: This subroutine cannot be used if you want to have "knobs" (in the PTC sense).
! Note: Use the routine get_ptc_param to get the state of PTC parameters.
!
! Modules needed:
!   use ptc_interface_mod
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
!                     Default = true.
!                     See the PTC guide for more details.
!-

subroutine set_ptc (e_tot, particle, taylor_order, integ_order, &
                                  n_step, no_cavity, exact_modeling, exact_misalign) 

use mad_like, only: make_states, exact_model, always_exactmis, &
              assignment(=), nocavity, default, operator(+), &
              berz, init, set_madx, lp, superkill, TIME0, PHASE0, HIGHEST_FRINGE
use madx_ptc_module, only: ptc_ini_no_append, append_empty_layout, m_u, bmadl

implicit none

integer, optional :: integ_order, particle, n_step, taylor_order
integer this_method, this_steps
integer nd2, t_order

real(rp), optional :: e_tot
real(rp), save :: old_e_tot = 0
real(dp) this_energy

logical, optional :: no_cavity, exact_modeling, exact_misalign
logical, save :: init_needed = .true., init2_needed = .true.
logical params_present

character(16) :: r_name = 'set_ptc'

! ptc cannot be used with photons

if (present(particle)) then
  if (particle == photon$) return
endif

! Some init

if (init2_needed) then
  EXACT_MODEL = .false.
  ALWAYS_EXACTMIS = .true.
  init2_needed = .false.
endif

! More init

HIGHEST_FRINGE = bmad_com%ptc_max_fringe_order

! do not call set_mad

params_present = present(e_tot) .and. present(particle)

if (init_needed .and. params_present) then
  if (particle == positron$ .or. particle == electron$) then
    call make_states(.true._lp)
  else
    call make_states(.false._lp)
  endif
  ! Use PTC time tracking
  DEFAULT = DEFAULT + TIME0
  PHASE0 = 0
endif

if (present (exact_modeling))     EXACT_MODEL = exact_modeling
if (present (exact_misalign)) ALWAYS_EXACTMIS = exact_misalign
if (present(no_cavity))       DEFAULT = DEFAULT+NOCAVITY

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
  if (init_needed .or. old_e_tot /= e_tot .or. present(integ_order) .or. present(n_step)) then
    this_energy = 1e-9 * e_tot
    if (this_energy == 0) then
      call out_io (s_fatal$, r_name, 'E_TOT IS 0.')
      if (global_com%exit_on_error) call err_exit
    endif
    call set_madx (energy = this_energy, method = this_method, step = this_steps)
    old_e_tot  = e_tot
    ! Only do this once
    if (init_needed) then
      call ptc_ini_no_append 
    endif
    init_needed = .false.
  endif
endif

! Do not call init before the call to make_states

if (.not. init_needed) then  ! If make_states has been called
  t_order = 0
  if (present(taylor_order)) t_order = taylor_order
  if (t_order == 0) t_order = bmad_com%taylor_order
  if (t_order == 0) t_order = ptc_com%taylor_order_saved
  if (ptc_com%taylor_order_ptc /= t_order) then
    call init (DEFAULT, t_order, 0, berz, nd2, ptc_com%real_8_map_init)
    ptc_com%taylor_order_ptc = t_order
  endif
endif

! Superkill tells PTC to do a through cleanup when killing a fibre.

SUPERKILL = .false.

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
! Subroutine taylor_equal_real_8 (bmad_taylor, y8)
!
! Subroutine to convert from a real_8 taylor map in Etienne's PTC 
! to a taylor map in Bmad. This does not do any
! conversion between Bmad units (z, dp/p0) and PTC units (dE/p0, c*t).
!
! Subroutine overloads "=" in expressions
!       bmad_taylor = y8
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   y8(:) -- real_8: PTC Taylor map.
!
! Output:
!   bmad_taylor(:) -- Taylor_struct: Input taylor map.
!-

subroutine taylor_equal_real_8 (bmad_taylor, y8)

use polymorphic_taylor, only: assignment (=), universal_taylor, real_8
use definition, only: real_8

implicit none

type (real_8), intent(in) :: y8(:)
type (taylor_struct), intent(inout) :: bmad_taylor(:)
type (universal_taylor) :: u_t(size(y8))

integer i, n_taylor

!

n_taylor = size(y8)
do i = 1, n_taylor
  u_t(i) = 0  ! nullify
  u_t(i) = y8(i)%t
enddo

call universal_to_bmad_taylor (u_t, bmad_taylor)

do i = 1, n_taylor
  u_t(i) = -1  ! deallocate
enddo

end subroutine taylor_equal_real_8

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine real_8_equal_taylor (y8, bmad_taylor)
!
! Subroutine to convert from a taylor map in Bmad to a
! real_8 taylor map in Etienne's PTC. This does not do any
! conversion between Bmad units (z, dp/p0) and PTC units (dE/p0, c*t).
! To convert coordinates, use the taylor_to_real_8 routine.
!
! Subroutine overloads "=" in expressions
!       y8 = bmad_taylor
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   bmad_taylor(:) -- Taylor_struct: Input taylor map.
!
! Output:
!   y8(:) -- real_8: PTC Taylor map.
!-

subroutine real_8_equal_taylor (y8, bmad_taylor)

use polymorphic_taylor, only: kill, assignment(=), real_8, universal_taylor

implicit none

type (real_8), intent(inout) :: y8(:)
type (taylor_struct), intent(in) :: bmad_taylor(:)
type (universal_taylor) :: u_t

integer i, j, n, n_taylor

logical switch

! init

call kill (y8)
call real_8_init (y8, .true.)

!

n_taylor = size(bmad_taylor)
do i = 1, n_taylor

  switch = .true.

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

end subroutine real_8_equal_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_to_real_8 (bmad_taylor, v0, beta0, beta1, y8)
!
! Routine to convert a Bmad Taylor map to PTC real_8 map.
! The conversion includes the conversion between Bmad and PTC time coordinate systems.
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Input taylor map.
!   v0(6)    -- real(rp): Initial vector (in Bmad coords) around which bmad_taylor was made.
!   beta0    -- real(rp): Reference particle velocity at beginning of map
!   beta1    -- real(rp): Reference particle velocity at end of map
!
! Output:
!   y8(6) -- real_8: PTC Taylor map.
!-

subroutine taylor_to_real_8 (bmad_taylor, v0, beta0, beta1, y8)

use s_fibre_bundle

implicit none

type (taylor_struct) :: bmad_taylor(:)
type (real_8) y8(:), rr(6), bet
type (damap) bm, id, si

real(dp) beta0, beta1, fix0(6), v0(:), v_ptc(6)

!

call alloc(bm, id, si)
call alloc(rr)
call alloc (bet)

y8 = bmad_taylor
bm = y8

fix0 = bm     ! Save zeroth order terms
id = 1        ! Identity map

call vec_bmad_to_ptc (v0, beta0, v_ptc)
rr = id + v_ptc

y8 = rr
y8(6) = (2.d0*rr(5)/beta0+rr(5)**2)/(sqrt(1.d0+2.d0*rr(5)/beta0+rr(5)**2)+1.d0)
bet = (1.d0+y8(6))/(1.d0/beta0+rr(5))
y8(5) = -bet*rr(6)

si = y8        ! ptc to bmad map

bm = bm*si
bm = fix0

rr = bm
y8 = rr
y8(5) = (rr(6)**2+2.d0*rr(6))/(1.d0/beta1+sqrt( 1.d0/beta1**2+rr(6)**2+2.d0*rr(6)) )
bet = (1.d0+rr(6))/(1.d0/beta1+y8(5))
y8(6) = -rr(5)/bet

call kill(bm, id, si)
call kill(bet)
call kill(rr)

end subroutine taylor_to_real_8

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine real_8_to_taylor (y8, v0, beta0, beta1, bmad_taylor)
!
! Routine to convert a PTC real_8 map to a Bmad Taylor map.
! The conversion includes the conversion between Bmad and PTC time coordinate systems.
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   y8(6) -- real_8: PTC Taylor map. NOTE: y8 is used as scratch space and therefore trashed.
!   v0(6) -- real(rp): Initial vector (in Bmad coords) around which y8 was made.
!   beta0 -- real(rp): Reference particle velocity at beginning of map
!   beta1 -- real(rp): Reference particle velocity at end of map
!
! Output:
!   bmad_taylor(6) -- Taylor_struct: Bmad Taylor map.
!-

subroutine real_8_to_taylor (y8, v0, beta0, beta1, bmad_taylor)

use s_fibre_bundle

implicit none

type (taylor_struct) :: bmad_taylor(:)
type (real_8) y8(:), rr(6), bet
type (damap) bm, id, si

real(rp) beta0, beta1, fix0(6), v0(:)

!

call alloc (bm, id, si)
call alloc (rr)
call alloc (bet)

bm = y8

fix0 = bm
id = 1

rr = id+v0

y8 = rr 
y8(5) = (rr(6)**2+2.d0*rr(6))/(1.d0/beta0 + sqrt(1.d0/beta0**2+rr(6)**2+2.d0*rr(6)))
bet = (1.d0+rr(6))/(1.d0/beta0+y8(5))
y8(6) = -rr(5)/bet

si=y8  ! bmad to ptc map

bm = bm*si
bm = fix0

rr = bm
y8 = rr
y8(6) = (2.d0*rr(5)/beta1+rr(5)**2)/(sqrt(1.d0+2.d0*rr(5)/beta1+rr(5)**2)+1.d0)
bet = (1.d0+y8(6))/(1.d0/beta1+rr(5))
y8(5) = -bet*rr(6)

bmad_taylor = y8

call kill (rr)
call kill (bet)
call kill (bm, id, si)

end subroutine real_8_to_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_bmad_to_ptc (taylor_bmad, beta0, taylor_ptc, mat6)
!
! Routine to convert a Bmad taylor map to PTC taylor.
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   taylor_bmad(6) -- type(taylor_struct): Bmad Taylor.
!   beta0          -- real(rp): Reference particle velocity
!
! Output:
!   taylor_ptc(6)  -- type(real_8): PTC coordinates.
!-

subroutine taylor_bmad_to_ptc (taylor_bmad, beta0, taylor_ptc)

use s_fibre_bundle

implicit none

type (taylor_struct) taylor_bmad(:)
type (real_8) taylor_ptc(:), t_ptc(6)

real(rp) beta0

! taylor_ptc(5) = (E - E0) / P0c
! taylor_ptc(6) = c (t - t0)
! 1/beta0 + taylor_ptc(5) == E / P0c

t_ptc = taylor_bmad
taylor_ptc = t_ptc
taylor_ptc(5) = (t_ptc(6)**2 + 2*t_ptc(6)) / (1/beta0 + sqrt(1/beta0**2+t_ptc(6)**2+2*t_ptc(6)))
taylor_ptc(6) = -t_ptc(5) * (1/beta0 + taylor_ptc(5)) / (1 + t_ptc(6))

end subroutine taylor_bmad_to_ptc 

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_ptc_to_bmad (taylor_ptc, beta0, taylor_bmad)
!
! Routine to convert a PTC real_8 taylor to a Bmad Taylor.
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   taylor_ptc(6)  -- real_8: PTC taylor.
!   beta0          -- real(rp): Reference particle velocity
!
! Output:
!   taylor_bmad(6) -- taylor_struct: Bmad Taylor.
!-

subroutine taylor_ptc_to_bmad (taylor_ptc, beta0, taylor_bmad)

use s_fibre_bundle

implicit none

type (taylor_struct) taylor_bmad(:)
type (real_8) taylor_ptc(:), t_ptc(6)
real(rp) beta0

!

t_ptc(6) = (2*taylor_ptc(5)/beta0+taylor_ptc(5)**2)/(sqrt(1+2*taylor_ptc(5)/beta0+taylor_ptc(5)**2)+1)
t_ptc(5) = -taylor_ptc(6) * (1 + t_ptc(6)) / (1/beta0 + taylor_ptc(5))

taylor_bmad = t_ptc

end subroutine taylor_ptc_to_bmad 

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine vec_bmad_to_ptc (vec_bmad, beta0, vec_ptc, conversion_mat)
!
! Routine to convert a Bmad vector map to PTC vector,
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   vec_bmad(6) -- real(rp): Bmad coordinates.
!   beta0       -- real(rp): Reference particle velocity
!
! Output:
!   vec_ptc(6)     -- real(rp): PTC coordinates.
!   conversion_mat -- real(rp), optional: Jacobian matrix of Bmad -> PTC conversion map.
!-

subroutine vec_bmad_to_ptc (vec_bmad, beta0, vec_ptc, conversion_mat)

implicit none

real(rp) vec_bmad(:), vec_ptc(:)
real(rp) beta0
real(rp), optional :: conversion_mat(6,6)
real(rp) factor1, factor2

! vec_ptc(5) = (E - E0) / P0c
! vec_ptc(6) = c (t - t0)
! 1/beta0 + vec_ptc(5) == E / P0c

vec_ptc = vec_bmad
vec_ptc(5) = (vec_bmad(6)**2 + 2*vec_bmad(6)) / (1/beta0 + sqrt(1/beta0**2+vec_bmad(6)**2+2*vec_bmad(6)) )
vec_ptc(6) = -vec_bmad(5) * (1/beta0 + vec_ptc(5)) / (1 + vec_bmad(6))

if (present(conversion_mat)) then
  call mat_make_unit(conversion_mat)
  factor1 = sqrt(1/beta0**2+vec_bmad(6)**2+2*vec_bmad(6))
  factor2 = 1+beta0**2*vec_bmad(6)*(2+vec_bmad(6))
  conversion_mat(5,5) = 0
  conversion_mat(5,6) = beta0**2*(1+vec_bmad(6))*factor1/factor2
  conversion_mat(6,5) = -(1/beta0+beta0*vec_bmad(6)*(2+vec_bmad(6))/(1+beta0*factor1))/(1+vec_bmad(6))
  conversion_mat(6,6) = -((beta0**2-1)*vec_bmad(5)*factor1)/((1+vec_bmad(6))**2*factor2)
end if

end subroutine vec_bmad_to_ptc 

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine vec_ptc_to_bmad (vec_ptc, beta0, vec_bmad, conversion_mat)
!
! Routine to convert a PTC orbit vector to a Bmad orbit vector.
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   vec_ptc(6)  -- real(rp): PTC coordinates.
!   beta0       -- real(rp): Reference particle velocity
!
! Output:
!   vec_bmad(6)    -- real(rp): Bmad coordinates.
!   conversion_mat -- Real(rp), optional: Jacobian matrix of PTC -> Bmad conversion map.
!-

subroutine vec_ptc_to_bmad (vec_ptc, beta0, vec_bmad, conversion_mat)

implicit none

real(rp) vec_bmad(:), vec_ptc(:)
real(rp) beta0
real(rp), optional :: conversion_mat(6,6)
real(rp) factor1, factor2

!

vec_bmad = vec_ptc
vec_bmad(6) = (2*vec_ptc(5)/beta0+vec_ptc(5)**2)/(sqrt(1+2*vec_ptc(5)/beta0+vec_ptc(5)**2)+1)
vec_bmad(5) = -vec_ptc(6) * (1 + vec_bmad(6)) / (1/beta0 + vec_ptc(5))

if (present(conversion_mat)) then
  call mat_make_unit(conversion_mat)
  factor1 = sqrt(1+2*vec_ptc(5)/beta0+vec_ptc(5)**2)
  factor2 = beta0+2*vec_ptc(5)+beta0*vec_ptc(5)**2 
  conversion_mat(5,5) = beta0*(beta0**2-1)*factor1*vec_ptc(6)/((1+beta0*vec_ptc(5))**2*factor2)
  conversion_mat(5,6) = -(1+vec_ptc(5)*(2+beta0*vec_ptc(5))/(beta0*(1+factor1)))/(1/beta0+vec_ptc(5))
  conversion_mat(6,5) = (1+beta0*vec_ptc(5))*factor1/factor2
  conversion_mat(6,6) = 0
end if

end subroutine vec_ptc_to_bmad 

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
! Modules needed:
!   use ptc_interface_mod
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
! Modules needed:
!   use ptc_interface_mod
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
! Subroutine real_8_init (y, set_taylor)
!
! Subroutine to allocate a PTC real_8 variable.
! The internal kind parameter will be set to 0.
!
! Note: If this variable has been used before, make sure you have 
! deallocated using:
!   call kill(y)
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   y(:)       -- Real_8: 
!   set_taylor -- Logical, optional :: If present and True then make
!                   y the identity taylor series (kind = 2).
!
! Output:
!   y(:) -- Real_8: Identity map.
!-

subroutine real_8_init (y, set_taylor)

use s_fibre_bundle, only: assignment(=), alloc, real_8

implicit none

type (real_8) :: y(:)
real(dp) :: x(size(y))

logical, optional :: set_taylor

!

call alloc(y)
y = ptc_com%real_8_map_init

if (present(set_taylor)) then
  x = 0
  if (set_taylor) y = x   ! converts y to taylor (kind = 2)
endif

end subroutine real_8_init

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Elemental subroutine universal_to_bmad_taylor (u_taylor, bmad_taylor)
!
! Subroutine to convert from a universal_taylor map in Etienne's PTC 
! to a taylor map in Bmad.
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   u_taylor(:) -- Universal_taylor: Universal_taylor map.
!
! Output:
!   bmad_taylor(:)   -- Taylor_struct:
!-

elemental subroutine universal_to_bmad_taylor (u_taylor, bmad_taylor)

use definition, only: universal_taylor

implicit none

type (universal_taylor), intent(in) :: u_taylor
type (taylor_struct), intent(inout) :: bmad_taylor

integer :: j, k, n

! Remember to suppress any terms that have a zero coef.  

if (associated(bmad_taylor%term)) deallocate(bmad_taylor%term)

n = count(u_taylor%c(:) /= 0)
allocate(bmad_taylor%term(n))

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
! Modules needed:
!   use ptc_interface_mod
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
! Modules needed:
!   use ptc_interface_mod
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
! Modules needed:
!   use ptc_interface_mod
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
     complex_taylor%ref = cmplx(re_taylor%ref, 0.0_rp)
     do n = 1, n1
       complex_taylor%term(n)%coef = cmplx(re_taylor%term(n)%coef, 0.0_rp)
       complex_taylor%term(n)%expn = re_taylor%term(n)%expn
     enddo
 
  else if (n2 > 0) then
  ! purely imaginary
    call init_complex_taylor_series (complex_taylor, n2)
    complex_taylor%ref = cmplx(0.0_rp, im_taylor%ref)
    do n = 1, n2
      complex_taylor%term(n)%coef = cmplx(0.0_rp, im_taylor%term(n)%coef)
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
complex_taylor%ref = cmplx(taylor1%ref, taylor2%ref)

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
  complex_taylor%term(n)%coef = cmplx(re, im)
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
! Subroutine concat_real_8 (y1, y2, y3)
!
! Subroutine to concatinate two real_8 taylor maps.
!       y3 = y2(y1)
! This subroutine assumes that y1, y2, and y3 have been allocated.
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   y1(:) -- real_8: First Input map.
!   y2(:) -- real_8: Second Input map.
!
! Output
!   y3(:) -- real_8: Concatinated map.
!-

subroutine concat_real_8 (y1, y2, y3)

use s_fitting, only: alloc, assignment(=), kill, damap, operator(.o.), real_8

implicit none

type (real_8), intent(in) :: y1(:), y2(:)
type (real_8), intent(inout) :: y3(:)
type (damap) da1, da2, da3

! set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

! Allocate temp vars

call alloc(da1)
call alloc(da2)
call alloc(da3)

! concat

da1 = y1
da2 = y2

da3 = da2 .o. da1  ! concat with constant terms

y3 = da3

! kill temp vars

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
! Subroutine taylor_inverse (taylor_in, taylor_inv, err, ref_pt)
!
! Subroutine to invert a taylor map. Since the inverse map is truncated, it
! is not exact and the reference point about which it is computed matters.
!
! Moudules needed:
!   use ptc_interface_mod
!
! Input:
!   taylor_in(:)  -- Taylor_struct: Input taylor map.
!   ref_pt(:)     -- Real(rp), optional: Reference point about which the 
!                     inverse is taken. Default is zero.
!
! Output:
!   taylor_inv(:) -- Taylor_struct: Inverted taylor map.
!   err           -- Logical, optional: Set True if there is no inverse.
!                     If not present then print an error message.
!-

subroutine taylor_inverse (taylor_in, taylor_inv, err, ref_pt)

use s_fitting, only: assignment(=), alloc, kill, operator(**), damap, real_8

implicit none

type (taylor_struct) :: taylor_in(:)
type (taylor_struct) :: taylor_inv(:)
type (taylor_struct) tlr(size(taylor_in)), tlr2(size(taylor_in))
type (real_8) y(size(taylor_in)), yc(size(taylor_in))
type (damap) da

real(rp), optional :: ref_pt(:)
real(rp) c0(size(taylor_in))

integer i, n_taylor, expn(size(taylor_in))

real(dp) c8(size(taylor_in)), c_ref(size(taylor_in))

logical, optional :: err

character(16) :: r_name = 'taylor_inverse'

! Set the taylor order in PTC if not already done so.

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

call alloc(da)
call alloc(y)

! If ref_pt is present then shift the map to the new reference point.
! Also the inverse operation of PTC ignores constant terms so we have 
! to take them out and then put them back in.

if (present(ref_pt)) then
  y = taylor_in
  call real_8_init(yc)
  c_ref = ref_pt
  yc = c_ref
  call concat_real_8 (yc, y, y)
  tlr2 = y
  call kill (yc)
  call kill (y)
  call remove_constant_taylor (tlr2, tlr, c0, .true.)
  call kill_taylor (tlr2)
else
  call remove_constant_taylor (taylor_in, tlr, c0, .true.)
endif

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

! Put constant terms back in.
! If the Map is written as:
!   R_out = T * R_in + C
! Then inverting:
!   R_in = T_inv * (R_out - C) = T_inv * (I - C) * R_out

if (any(c0 /= 0)) then
  call real_8_init(yc)
  c8 = c0
  yc = -c8                       ! Convert this to taylor map: I - c8
  call concat_real_8 (yc, y, y) 
  call kill (yc)
  taylor_inv%ref = c0
endif

! Transfer inverse to taylor_inv.

taylor_inv = y

! Take out the ref_pt offset if needed

if (present(ref_pt)) then
  do i = 1, n_taylor
    call add_taylor_term (taylor_inv(i), ref_pt(i))
  enddo
endif

! Clean up

call kill (da)
call kill (y)
call kill (yc)
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
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   taylor1(:) -- Taylor_struct: Taylor map.
!   taylor2(:) -- Taylor_struct: Taylor map.
!
! Output
!   taylor3(:) -- Taylor_struct: Concatinated map
!-

subroutine concat_taylor (taylor1, taylor2, taylor3)

use s_fitting, only: assignment(=), kill, real_8

implicit none

type (taylor_struct) :: taylor1(:), taylor2(:)
type (taylor_struct) :: taylor3(:)
type (real_8) y1(size(taylor1)), y2(size(taylor1)), y3(size(taylor1))

! Set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

! Allocate temp vars

call real_8_init (y1)
call real_8_init (y2)
call real_8_init (y3)

! Concat

y1 = taylor1
y2 = taylor2

call concat_real_8 (y1, y2, y3)

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
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   taylor1(6) -- Taylor_struct: Taylor map.
!   ele        -- ele_struct: Element containing the taylor map.
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

! Init

call real_8_init(x_ele)
call real_8_init(x_body)
call real_8_init(x1)
call real_8_init(x3)

beta0 = ele%value(p0c_start$)/ele%value(e_tot_start$)
beta1 = ele%value(p0c$)/ele%value(e_tot$)

! Create a PTC fibre that holds the misalignment info
! and create map corresponding to ele%taylor.

param%particle = positron$  ! Actually this does not matter to the calculation
call ele_to_fibre (ele, fib, param, use_offsets = .true.)

x_dp = 0
x_ele = x_dp  ! x_ele = Identity map 

if (ele%key == sbend$) then
  tilt = ele%value(ref_tilt_tot$)
else
  tilt = ele%value(tilt_tot$)
endif

call dtiltd (1, tilt, 1, x_ele)
call mis_fib (fib, x_ele, DEFAULT, .true., entering = .true.)

call taylor_to_real_8 (ele%taylor, ele%taylor%ref, beta0, beta1, x_body)

call concat_real_8 (x_ele, x_body, x_ele)
call mis_fib (fib, x_ele, DEFAULT, .true., entering = .false.)
call dtiltd (1, tilt, 2, x_ele)

! Concat with taylor1

call taylor_to_real_8 (taylor1, taylor1%ref, beta0, beta0, x1)
call concat_real_8 (x1, x_ele, x3)

! convert x3 to final result taylor3

call real_8_to_taylor(x3, taylor1%ref, beta0, beta1, taylor3)
taylor3(:)%ref = taylor1(:)%ref

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
! Subroutine taylor_propagate1 (bmad_taylor, ele, param)
!
! Subroutine to track (symplectic integration) a taylor map through an element.
! The alternative routine, if ele has a taylor map, is concat_taylor.
!
! This routine will fail if there is no corresponding ptc fibre for this
! element. In general, the transfer_map_calc routine should be used instead.
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Map to be tracked
!   ele            -- Ele_struct: Element to track through
!   param          -- lat_param_struct: 
!
! Output:
!   bmad_taylor(6)  -- Taylor_struct: Map through element
!-

subroutine taylor_propagate1 (bmad_taylor, ele, param)

use s_tracking
use mad_like, only: real_8, fibre, ptc_track => track

implicit none

type (taylor_struct) bmad_taylor(:)
type (real_8), save :: ptc_tlr(6)
type (ele_struct) ele, drift_ele
type (lat_param_struct) param
type (fibre), pointer :: ptc_fibre

real(rp) beta0, beta1, m2_rel

! set the taylor order in PTC if not already done so

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

! Init ptc map with bmad map

beta0 = ele%value(p0c_start$) / ele%value(e_tot_start$)
beta1 = ele%value(p0c$) / ele%value(e_tot$)

call real_8_init (ptc_tlr)
call taylor_to_real_8 (bmad_taylor, bmad_taylor%ref, beta0, beta0, ptc_tlr)

! Track entrance drift if PTC is using a hard edge model

if (tracking_uses_end_drifts(ele)) then
  call create_hard_edge_drift (ele, upstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, ptc_fibre, param, .true.)
  call ptc_track (ptc_fibre, ptc_tlr, default)  ! "track" in PTC
endif

! Init ptc "element" (fibre) and track the map

call ele_to_fibre (ele, ptc_fibre, param, .true.)
call ptc_track (ptc_fibre, ptc_tlr, default)  ! "track" in PTC

! Track exit side drift if PTC is using a hard edge model

if (tracking_uses_end_drifts(ele)) then
  call create_hard_edge_drift (ele, downstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, ptc_fibre, param, .true.)
  call ptc_track (ptc_fibre, ptc_tlr, default)  ! "track" in PTC
endif

! transfer ptc map back to bmad map

call real_8_to_taylor(ptc_tlr, bmad_taylor%ref, beta0, beta1, bmad_taylor)

! cleanup

call kill (ptc_tlr)

end subroutine taylor_propagate1

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ele_to_taylor (ele, param, orb0, taylor_map_includes_offsets)
!
! Subroutine to make a taylor map for an element. 
! The order of the map is set by set_ptc
!
! Modules needed:
!   use ptc_interface_mod
!
! Input:
!   ele   -- Element_struct: 
!     %value(integrator_order$)  -- Order for the symplectic integrator: 2, 4, or 6.
!     %value(ds_step$)          -- Integrater step size.
!     %taylor_map_includes_offsets         -- Make Taylor map with element offsets, pitches, and tilt?
!   orb0  -- Coord_struct, optional: Starting coords around which the Taylor map is evaluated.
!   param -- lat_param_struct: 
!     %e_tot -- Needed for wigglers.
!   taylor_map_includes_offsets -- Logical, optional: If present then overrides 
!                         ele%taylor_map_includes_offsets.
!
! Output:
!   ele -- Element_struct:
!     %taylor(6)  -- Taylor maps.
!-

subroutine ele_to_taylor (ele, param, orb0, taylor_map_includes_offsets)

use s_tracking
use mad_like, only: real_8, fibre, ptc_track => track

implicit none

type (ele_struct) :: ele, drift_ele
type (lat_param_struct) :: param
type (coord_struct), optional, intent(in) :: orb0
type (coord_struct) c0

type (fibre), pointer :: ptc_fibre
type (real_8) y0(6), y2(6), y8(6), bet

real(dp) x(6), beta

logical, optional :: taylor_map_includes_offsets
logical :: warning_given = .false.
logical use_offsets

character(16) :: r_name = 'ele_to_taylor'

! Init

if (ptc_com%taylor_order_ptc == 0) call set_ptc (taylor_order = bmad_com%taylor_order)

call alloc (y8)
call alloc (bet)

call attribute_bookkeeper (ele, param, .true.)

! Match elements are not implemented in PTC so just use the matrix.
! Also Taylor elements already have a taylor map.

if (ele%key == taylor$) return

if (ele%key == match$) then
  c0%vec = 0
  call make_mat6_bmad (ele, param, c0, c0, .true.)
  call mat6_to_taylor (ele%vec0, ele%mat6, ele%taylor)
  if (.not. warning_given) then
    call out_io (s_warn$, r_name, &
      'Note: Taylor maps for Match elements are always 1st order!')
    warning_given = .true.
  endif
  return
endif

! Initial map

use_offsets = logic_option(ele%taylor_map_includes_offsets, taylor_map_includes_offsets)
 

if (present(orb0)) then
  ele%taylor(:)%ref = orb0%vec
  x = orb0%vec  ! y = IdentityMap + const
else
  ele%taylor(:)%ref = 0
  x = 0
endif

call real_8_init(y0)
y0 = x ! y = IdentityMap + const

! Convert to PTC

beta = ele%value(p0c_start$) / ele%value(e_tot_start$)
y8 = y0
y8(5) = (y0(6)**2+2.d0*y0(6))/(1.d0/beta+sqrt( 1.d0/beta**2+y0(6)**2+2.d0*y0(6)) )
bet = (1.d0+y0(6))/(1.d0/beta+y8(5))
y8(6) = -y0(5)/bet

! Track entrance drift if PTC is using a hard edge model

if (tracking_uses_end_drifts(ele)) then
  call create_hard_edge_drift (ele, upstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, ptc_fibre, param, .true.)
  call ptc_track (ptc_fibre, y8, default) ! "track" in PTC
endif

! Track element

call ele_to_fibre (ele, ptc_fibre, param, use_offsets)
call ptc_track (ptc_fibre, y8, default) ! "track" in PTC

! Track exit end drift if PTC is using a hard edge model

if (tracking_uses_end_drifts(ele)) then
  call create_hard_edge_drift (ele, downstream_end$, drift_ele)
  call ele_to_fibre (drift_ele, ptc_fibre, param, .true.)
  call ptc_track (ptc_fibre, y8, default) ! "track" in PTC
endif

! PTC to Bmad

beta = ele%value(p0c$) / ele%value(e_tot$)
y0 = y8
y0(6) = (2.d0*y8(5)/beta+y8(5)**2)/(sqrt(1.d0+2.d0*y8(5)/beta+y8(5)**2)+1.d0)
bet = (1.d0+y0(6))/(1.d0/beta+y8(5))
y0(5) = -bet*y8(6)

! take out the offset

if (any(x /= 0)) then
  call real_8_init(y2)
  y2 = -x  ! y2 = IdentityMap - x
  call concat_real_8 (y2, y0, y0)
  call kill(y2)
endif

! convert to bmad_taylor  

ele%taylor = y0

call kill(y0)
call kill(y8)
call kill(bet)

if (associated (ele%ptc_genfield)) call kill_ptc_genfield (ele%ptc_genfield)

call set_ele_status_stale (ele, mat6_group$)

end subroutine ele_to_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine type_real_8_taylors (y)
!
! Subroutine to type out the taylor map from a real_8 array.
!
! Modules needed:
!   use ptc_interface_mod
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
! Modules needed:
!   use ptc_interface_mod
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
! Modules needed:
!   use ptc_interface_mod
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
! Subroutine ele_to_fibre (ele, ptc_fibre, param, use_offsets, integ_order, steps, for_layout)
!
! Routine to convert a Bmad element to a PTC fibre element.
!
! Note: You need to call set_ptc before using this routine.
!
! Modules Needed:
!   use ptc_interface_mod
!
! Input:
!   ele         -- Ele_struct: Bmad element.
!   param       -- lat_param_struct: 
!   use_offsets -- Logical: Does ptc_fibre include element offsets, pitches and tilt?
!   integ_order -- Integer, optional: Order for the 
!                    sympletic integrator. Possibilities are: 2, 4, or 6
!                    Overrides ele%value(integrator_order$).
!                    default = 2 (if not set with set_ptc).
!   steps       -- Integer, optional: Number of integration steps.
!                    Overrides ele%value(ds_step$).
!   for_layout  -- Logical, optional: If True then fibre will be put in the layout.
!                    Default is False.
!
! Output:
!   ptc_fibre -- Fibre: PTC fibre element.
!+

subroutine ele_to_fibre (ele, ptc_fibre, param, use_offsets, integ_order, steps, for_layout)

use madx_ptc_module

implicit none
 
type (ele_struct), target :: ele
type (lat_param_struct) param
type (ele_struct), pointer :: field_ele, ele2
type (fibre), pointer :: ptc_fibre
type (keywords) ptc_key
type (ele_pointer_struct), allocatable :: field_eles(:)
type (work) energy_work
type (el_list) ptc_el_list

real(rp), allocatable :: dz_offset(:)
real(rp) leng, hk, vk, s_rel, z_patch, phi_tot
real(rp), pointer :: val(:)
real(rp), target, save :: value0(num_ele_attrib$) = 0

integer i, n, key, n_term, exception, n_field, ix
integer, optional :: integ_order, steps

logical use_offsets
logical, optional :: for_layout

character(16) :: r_name = 'ele_to_fibre'

!

call zero_key(ptc_key)  ! init key

select case (ele%ptc_integration_type)
case (drift_kick$);  ptc_key%model = 'DRIFT_KICK'
case (matrix_kick$); ptc_key%model = 'MATRIX_KICK'
case (ripken_kick$); ptc_key%model = 'DELTA_MATRIX_KICK'
end select


leng = ele%value(l$)

ptc_key%list%name = ele%name
ptc_key%list%l    = leng

if (use_offsets) then
  if (ele%key == sbend$) then
    ptc_key%tiltd = ele%value(ref_tilt_tot$)
  else
    ptc_key%tiltd = ele%value(tilt_tot$)
  endif
else
  ptc_key%tiltd = 0
endif

if (present(steps)) then
  ptc_key%nstep = steps
elseif (leng == 0) then
  ptc_key%nstep = 1
else
  if (ele%value(ds_step$) == 0) then
    call out_io (s_fatal$, r_name, 'DS_STEP IS ZERO FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
  endif
  ptc_key%nstep = nint(abs(leng) / ele%value(ds_step$))
  if (ptc_key%nstep == 0) ptc_key%nstep = 1
endif

ptc_key%method = nint(ele%value(integrator_order$))
if (ptc_key%method == 0) ptc_key%method = bmad_com%default_integ_order 
if (present(integ_order)) ptc_key%method = integ_order

!

key = ele%key
if (ele%is_on) then
  val => ele%value
else
  val => value0  ! Not is on then has zero strength.
endif

select case (key)

case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 
  ptc_key%magnet = 'drift'

case (quadrupole$) 
  ptc_key%magnet = 'quadrupole'

case (sad_mult$)
  if (ele%value(l$) == 0) then
    ptc_key%magnet = 'multipole'  ! No Bz field
  else
    ptc_key%magnet = 'solenoid'
    ptc_key%list%bsol = val(ks$)
  endif

case (sbend$) 
  ptc_key%magnet = 'sbend'
  ptc_key%list%b0   = ele%value(g$) * leng ! Yep this is correct. 
  ptc_key%list%t1   = ele%value(e1$)
  ptc_key%list%t2   = ele%value(e2$)
  ptc_key%list%hgap = ele%value(hgap$)
  ptc_key%list%fint = ele%value(fint$)

  if (ele%value(f1$) /= 0) then
    ptc_key%list%fint = ele%value(f1$) / 6
    ptc_key%list%hgap = 0.5
  endif

  if (ele%value(fintx$) /= ele%value(fint$)) then
    call out_io (s_error$, r_name, &
        'FINT AND FINTX ARE NOT THE SAME FOR BEND: ' // ele%name, &
        'PTC CANNOT HANDLE THIS!')
  endif

  if (ele%value(hgapx$) /= ele%value(hgap$)) then
    call out_io (s_error$, r_name, &
        'HGAP AND HGAPX ARE NOT THE SAME FOR BEND: ' // ele%name, &
        'PTC CANNOT HANDLE THIS!')
  endif

  ix = nint(ele%value(ptc_field_geometry$))
  if (ix == straight$) then
    ptc_key%magnet = 'wedgrbend'
    ptc_key%list%t1   = ele%value(e1$) - ele%value(angle$)/2
    ptc_key%list%t2   = ele%value(e2$) - ele%value(angle$)/2
  elseif (ix == true_rbend$) then
    ptc_key%magnet = 'truerbend'
    ptc_key%list%t1   = ele%value(e1$) - ele%value(angle$)/2
    !! ptc_key%list%t2   = ele%value(e2$) - ele%value(angle$)/2 ! Determined by %t1 in this case.
  endif

case (sextupole$)
  ptc_key%magnet = 'sextupole'

case (octupole$)
  ptc_key%magnet = 'octupole'

case (solenoid$)
  ptc_key%magnet = 'solenoid'
  ptc_key%list%bsol = val(ks$)

case (sol_quad$)
  ptc_key%magnet = 'solenoid'
  ptc_key%list%bsol = val(ks$)

case (marker$, detector$, fork$, photon_fork$, beginning_ele$, patch$, floor_shift$, fiducial$)
  ptc_key%magnet = 'marker'
  ptc_key%nstep = 1

case (kicker$, hkicker$, vkicker$)
  ptc_key%magnet = 'kicker'

! 

case (rfcavity$, lcavity$)
  if (ele%value(rf_frequency$) == 0) then
    call out_io (s_fatal$, r_name, 'RF FREQUENCY IS ZERO FOR: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  if (is_true(ele%value(traveling_wave$))) then
    ptc_key%magnet = 'twcavity'
    ptc_key%list%volt = 1e-6 * e_accel_field(ele, voltage$)
  else
    ptc_key%magnet = 'rfcavity'
    ptc_key%list%volt = 2e-6 * e_accel_field(ele, voltage$)
    ptc_key%list%n_bessel = -1   ! Triggers Bmad compatible cavity.
  endif

  ptc_key%list%freq0 = ele%value(rf_frequency$)
  phi_tot = ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_err$) + ele%value(phi0_ref$)
  if (tracking_uses_end_drifts(ele)) ptc_key%list%l = hard_edge_model_length(ele)

  if (ele%key == lcavity$) then
    ptc_key%list%lag = pi / 2 - twopi * phi_tot
  else
    ptc_key%list%lag = twopi * phi_tot
  endif

  ptc_key%list%delta_e = 0     ! For radiation calc.
  ptc_key%list%cavity_totalpath = 1  ! 

case (elseparator$)
  ptc_key%magnet = 'elseparator'
  hk = val(hkick$) / leng
  vk = val(vkick$) / leng
  if (hk == 0 .and. vk == 0) then
    ptc_key%tiltd = 0
  else
    if (param%particle < 0) then
      hk = -hk
      vk = -vk
    endif
    ptc_key%tiltd = -atan2 (hk, vk) + ele%value(tilt_tot$)
  endif
  ptc_key%list%volt = 1e-6 * ele%value(e_tot$) * sqrt(hk**2 + vk**2)

case (ab_multipole$, multipole$)
  ptc_key%magnet = 'multipole'

case (beambeam$)
  ptc_key%magnet = 'beambeam'
  call out_io (s_fatal$, r_name,  'BEAMBEAM ELEMENT NOT YET IMPLEMENTED!')
  if (global_com%exit_on_error) call err_exit

case (wiggler$, undulator$)
  ptc_key%magnet = 'wiggler'

case default
  call out_io (s_fatal$, r_name,  'UNKNOWN ELEMENT CLASS: ' // key_name(ele%key), &
                                  'FOR ELEMENT: ' // trim(ele%name))
  if (global_com%exit_on_error) call err_exit

end select

! Fringe

if (attribute_index(ele, 'FRINGE_TYPE') > 0) then  ! If fringe_type is a valid attribute
  ix = nint(ele%value(fringe_type$))
  ptc_key%list%bend_fringe = (ix == full_bend$ .or. ix == basic_bend$)

  select case (ix)
  case (none$)
    ptc_key%list%permfringe = 0
  case (basic_bend$)
    ptc_key%list%permfringe = 0
  case (full_straight$)
    ptc_key%list%permfringe = 1
  case (full_bend$)
    ptc_key%list%permfringe = 1
  case (sad_nonlin_only$)
    ptc_key%list%permfringe = 1
  case (sad_linear$)
    ptc_key%list%permfringe = 2
  case (sad_full$)
    ptc_key%list%permfringe = 3
  end select

  if (ele%key == sad_mult$ .and. ele%value(l$) == 0) ptc_key%list%permfringe = 0
endif

if (attribute_index(ele, 'FRINGE_AT') > 0) then  ! If fringe_at is a valid attribute
  ix = nint(ele%value(fringe_at$))
  ptc_key%list%kill_ent_fringe = (ix == exit_end$ .or. ix == no_end$)
  ptc_key%list%kill_exi_fringe = (ix == entrance_end$ .or. ix == no_end$)
endif

! Multipole components

call ele_to_an_bn (ele, param, .true., ptc_key%list%k, ptc_key%list%ks, ptc_key%list%nmul)

! Create ptc_fibre
! EXCEPTION is an error_flag. Set to 1 if error. Never reset.

n = lielib_print(12)
lielib_print(12) = 0  ! No printing info messages

if (logic_option(.false., for_layout)) then
  call create_fibre_append (.true., m_u%end, ptc_key, EXCEPTION)   ! ptc routine
  ptc_fibre => m_u%end%end

else
  call set_madx (energy = ele%value(e_tot$), method = ptc_key%method , step = ptc_key%nstep)
  if (associated(bmadl%start)) then
    call kill (bmadl)
    call set_up(bmadl)
  endif
  call create_fibre_append (.false., bmadl, ptc_key, EXCEPTION)   ! ptc routine
  ptc_fibre => bmadl%start
  bmadl%closed=.true.
  call ring_l(bmadl, .true.)
  call survey(bmadl)
  call make_node_layout (bmadl)
endif

ptc_fibre%dir = ele%orientation

lielib_print(12) = n

! sad_mult

if (ele%key == sad_mult$ .or. ele%key == quadrupole$) then
  if (ele%value(l$) /= 0) then
    ptc_fibre%mag%va  = ele%value(f1$)
    ptc_fibre%magp%va = ele%value(f1$)

    ptc_fibre%mag%vs  = ele%value(f2$)
    ptc_fibre%magp%vs = ele%value(f2$)
  endif

endif

! Set reference energy to the exit reference energy.

energy_work = 0
call find_energy (energy_work, p0c =  1d-9 * ele%value(p0c$))
ptc_fibre = energy_work

! wiggler

if (key == wiggler$ .or. key == undulator$) then

  call get_field_ele_list (ele, field_eles, dz_offset, n_field)
  do i = 1, n_field
    ele2 => field_eles(i)%ele
    if (ele2%key == ele%key) then
      s_rel = dz_offset(i)
      exit
    endif
  enddo

  if (hyper_x$ /= hyperbolic_xdollar .or. hyper_y$ /= hyperbolic_ydollar .or. &
                                        hyper_xy$ /= hyperbolic_xydollar) then
    print *, 'ERROR IN ELE_TO_FIBRE: WIGGLER FORM/TYPE MISMATCH!'
    print *, '     ', hyper_y$, hyper_xy$, hyper_x$
    print *, '     ', hyperbolic_ydollar, hyperbolic_xydollar, hyperbolic_xdollar
    if (global_com%exit_on_error) call err_exit
  endif

  n_term = size(ele2%wig%term)
  call init_sagan_pointers (ptc_fibre%mag%wi%w, n_term)   

  if (ele%is_on) then
    ptc_fibre%mag%wi%w%a(1:n_term) = c_light * ele2%value(polarity$) * ele2%wig%term%coef / ele%value(e_tot$)
  else
    ptc_fibre%mag%wi%w%a(1:n_term) = 0
  endif
  ptc_fibre%mag%wi%w%k(1,1:n_term)  = ele2%wig%term%kx
  ptc_fibre%mag%wi%w%k(2,1:n_term)  = ele2%wig%term%ky
  ptc_fibre%mag%wi%w%k(3,1:n_term)  = ele2%wig%term%kz
  ptc_fibre%mag%wi%w%f(1:n_term)    = ele2%wig%term%phi_z + s_rel * ele2%wig%term%kz
  ptc_fibre%mag%wi%w%form(1:n_term) = ele2%wig%term%type

  ! Correct z-position 

  z_patch = ele%value(delta_ref_time$) * c_light * ele%value(p0c$) / ele%value(e_tot$) - ele%value(l$)
  ptc_fibre%mag%wi%internal(6) = z_patch

  call copy (ptc_fibre%mag, ptc_fibre%magp)

endif

! Misalignments and patches...

call misalign_ele_to_fibre (ele, use_offsets, ptc_fibre)

! Set charge

if (associated(ele%branch)) then
  ptc_fibre%charge = ele%branch%param%rel_tracking_charge
endif

end subroutine ele_to_fibre

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
real(dp) omega(3), basis(3,3), angle(3)
real(rp) x_off, y_off, x_pitch, y_pitch, roll

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
  nullify(dummy_fibre%pos) ! To get around PTC initialization bug
  nullify(dummy_fibre%t1) ! To get around PTC initialization bug
  nullify(dummy_fibre%t2) ! To get around PTC initialization bug
  nullify(dummy_fibre%tm) ! To get around PTC initialization bug
  dummy_fibre = marker('dummy')
  nullify(dummy_fibre%loc)
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

    if (ele%value(e_tot_offset$) == 0) then
      ptc_fibre%patch%energy = 0
    else
      ptc_fibre%patch%energy = 4 ! Internal entrance patch. Must be done after find_patch call.
      ptc_fibre%patch%p0b = ele%value(p0c_start$) * 1d-9
      ptc_fibre%patch%b0b = beta_start  ! beta velocity
    endif

    ! PTC uses beta_start for the reference time calculation while Bmad uses beta_end so
    ! renormalize the patch length to get PTC to agree with Bmad.

    ptc_fibre%patch%time = 2     ! Subtract off reference time (which affects z in tracking).
    ptc_fibre%patch%b_t = ele%value(l$) + ele%value(t_offset$) * c_light * beta_end
  endif

!----------------------------------------------------------------------
! Not patch nor floor_shift element.

elseif (use_offsets) then

  ! Patch elements do not have misalignments

  if (ele%key == patch$ .or. ele%key == fiducial$ .or. ele%key == floor_shift$) return
  if (attribute_index(ele, 'X_OFFSET_TOT') < 1) return

  ! in PTC the reference point for the offsets is the beginning of the element.
  ! In Bmad the reference point is the center of the element..

  x_off = ele%value(x_offset_tot$)
  y_off = ele%value(y_offset_tot$)
  x_pitch = ele%value(x_pitch_tot$)
  y_pitch = ele%value(y_pitch_tot$)
  roll = 0
  if (ele%key == sbend$) roll = ele%value(roll_tot$)

  if (x_off /= 0 .or. y_off /= 0 .or. x_pitch /= 0 .or. y_pitch /= 0 .or. roll /= 0) then
    mis_rot = [x_off, y_off, 0.0_rp, -y_pitch, -x_pitch, roll]
    angle = 0
    angle(3) = -ptc_fibre%mag%p%tiltd
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
! Subroutine ele_to_an_bn (ele, param, creating_fibre, k, ks, n_max)
!
! Routine to compute the a(n) and b(n) multipole components of a magnet.
! This is used to interface between eles and PTC fibres
!
! Module needed:
!   ptc_interface_mod
!
! Input:
!   ele                 -- ele_struct: Bmad Element.
!   param               -- lat_param_struct: 
!   creating_fibre      -- integer: Set True if fibre is being created. False if fibre is being modified.
!
! Output:
!   k(1:n_pole_maxx+1)  -- real(rp): Skew multipole component.
!   ks(1:n_pole_maxx+1) -- real(rp): Normal multipole component.
!   n_max               -- integer: If creating fibre: Maximum non-zero multipole component.
!                                   If modifiying fibre: Maximum relavent multipole component.
!-

subroutine ele_to_an_bn (ele, param, creating_fibre, k, ks, n_max)

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param

real(rp) k(:), ks(:)
real(rp) cos_t, sin_t, leng, hk, vk, tilt
real(rp), pointer :: val(:)
real(rp), target, save :: value0(num_ele_attrib$) = 0
real(rp) an0(0:n_pole_maxx), bn0(0:n_pole_maxx)

integer n, n_max, key, n_relavent
logical creating_fibre, kick_here, has_nonzero_pole, add_kick

character(16) :: r_name = 'ele_to_an_bn'

!

leng = ele%value(l$)

key = ele%key
if (ele%is_on) then
  val => ele%value
else
  val => value0  ! Not is_on -> has zero strength.
endif

k = 0
ks = 0
n_max = 0
n_relavent = 0
add_kick = .true.

select case (key)

case (marker$, detector$, fork$, photon_fork$, beginning_ele$, em_field$, patch$, fiducial$, floor_shift$)
  return

case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$, rfcavity$, lcavity$, &
      ab_multipole$, multipole$, beambeam$, wiggler$, undulator$)
  ! Nothing to be done

case (octupole$)
  k(4) = val(k3$) / 6
  n_relavent = 4

case (quadrupole$) 
  k(2) = val(k1$)
  n_relavent = 2

case (sad_mult$)

case (sbend$)
  if (ele%is_on) then
    k(1) = ele%value(g_err$)
  else
    k(1) = -ele%value(g$)
  endif

  ! On ptc side k(1) is error field when creating a fibre but 
  ! is total field when fibre is being modified.

  if (.not. creating_fibre) k(1) = k(1) + ele%value(g$)

  k(2) = val(k1$)
  k(3) = val(k2$) / 2
  n_relavent = 3

case (sextupole$)
  k(3) = val(k2$) / 2
  n_relavent = 3

case (solenoid$)

case (sol_quad$)
  k(2) = val(k1$)
  n_relavent = 2

case (hkicker$, vkicker$)
  if (ele%key == hkicker$) k(1)  = val(kick$) 
  if (ele%key == vkicker$) ks(1) = val(kick$) 
  n_relavent = 1
  add_kick = .false.

case (kicker$)
  k(1)  = val(hkick$) 
  ks(1) = val(vkick$) 
  n_relavent = 1
  add_kick = .false.

case (elseparator$)
  call multipole_ele_to_ab (ele, param, .false., has_nonzero_pole, an0, bn0) 
  if (has_nonzero_pole) then
    call out_io (s_fatal$, r_name, 'MULTIPOLES IN AN ELSEPARATOR NOT SUPPORTED IN A FIBRE.')
    if (global_com%exit_on_error) call err_exit
  endif
  return

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN ELEMENT CLASS: ' // key_name(ele%key), &
                                 'FOR ELEMENT: ' // trim(ele%name))
  if (global_com%exit_on_error) call err_exit

end select

if (add_kick .and. has_hkick_attributes(ele%key) .and. &
                        (val(hkick$) /= 0 .or. val(vkick$) /= 0)) then
  hk = val(hkick$) / leng   ! PTC uses scaled kick for non-kicker elements.
  vk = val(vkick$) / leng
  if (ele%key == sbend$) then
    tilt = ele%value(ref_tilt_tot$)
  else
    tilt = ele%value(tilt_tot$)
  endif
  cos_t = cos(tilt)
  sin_t = sin(tilt)
  k(1)  = k(1)  - hk * cos_t - vk * sin_t
  ks(1) = ks(1) - hk * sin_t + vk * cos_t
  n_relavent = max(1, n_relavent)
endif

! bmad an and bn are integrated fields. PTC uses just the field.
! Exception is multipole element.

if (associated(ele%a_pole)) then
  call multipole_ele_to_ab (ele, param, .false., has_nonzero_pole, an0, bn0)
  if (leng /= 0) then
    an0 = an0 / leng
    bn0 = bn0 / leng
  endif

  n = min(n_pole_maxx+1, size(k))
  if (n-1 < n_pole_maxx) then
    if (any(an0(n:n_pole_maxx) /= 0) .or. any(bn0(n:n_pole_maxx) /= 0)) then
      print *, 'WARNING IN ELE_TO_FIBRE: MULTIPOLE NOT TRANSFERED TO FIBRE'
      print *, '        FOR: ', ele%name
    endif
  endif
   
  ks(1:n) = ks(1:n) + an0(0:n-1)
  k(1:n) = k(1:n) + bn0(0:n-1)
  n_relavent = n_pole_maxx
endif


if (creating_fibre) then
  do n = size(k), 1, -1
    if (ks(n) /= 0 .or. k(n) /= 0) exit
  enddo
  n_max  = n
else
  n_max = min(n_pole_maxx, size(k))
endif

end subroutine ele_to_an_bn

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine apply_patch_to_ptc_fibre (ele)
!
! Routine to take the patch parameters from a Bmad patch element and
! transfer them to the associated PTC fibre.
!
! Module needed:
!   use ptc_interface_mod
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

end module
