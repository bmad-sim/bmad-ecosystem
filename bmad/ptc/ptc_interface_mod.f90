!+
! Module ptc_interface_mod
!
! Module of basic PTC interface routines.
! Also see: ptc_layout_mod
!-

module ptc_interface_mod

use taylor_mod
use attribute_mod

interface assignment (=)
  module procedure damap_equal_bmad_taylor
  module procedure bmad_taylor_equal_damap
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

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine ptc_set_taylor_order_if_needed()
!
! Routine to see if the taylor_order for PTC needs to be set/changed.
! For example, for a change in bmad_com%taylor_order.
!-

subroutine ptc_set_taylor_order_if_needed ()

if ((bmad_com%taylor_order /= 0 .and. bmad_com%taylor_order /= ptc_private%taylor_order_ptc) .or. &
                                                                    ptc_private%taylor_order_ptc == 0) then
  call set_ptc (taylor_order = bmad_com%taylor_order)
endif

end subroutine ptc_set_taylor_order_if_needed

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

! 

call ptc_set_taylor_order_if_needed()

call alloc(y1)
call alloc(y2)
call alloc(y3)

y1 = taylor1
y2 = taylor2

do i = 1, size(taylor1)
  y3(i) = y1(i) + y2(i)
enddo

taylor3 = taylor_struct()   ! Should not need this but gfortran is unhappy without.
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

!

call ptc_set_taylor_order_if_needed()

call alloc(y1)
call alloc(y2)
call alloc(y3)

y1 = taylor1
y2 = taylor2

do i = 1, size(taylor1)
  y3(i) = y1(i) - y2(i)
enddo

taylor3 = taylor_struct()   ! Should not need this but gfortran is unhappy without.
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
!   intern_state -- Internal_state, optional: PTC state. If not present then the 
!                     ptc_private%base_state is used.
!
! Output:
!   lines(:)  -- character(100), optional, allocatable: Character array to hold the output.
!   n_lines   -- integer, optional: Number of lines used in lines(:)
!-

subroutine type_ptc_internal_state (intern_state, lines, n_lines)

implicit none

type (internal_state), optional, target :: intern_state
type (internal_state), pointer :: state_ptr

integer, optional :: n_lines
integer i, nl

character(*), allocatable, optional :: lines(:)
character(100), allocatable :: li(:)

!

state_ptr => ptc_private%base_state
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

    nl=nl+1; write (li(nl), '(2x, a, t40, 2a)') '%energy [Energy Patch at]:     ', patch_name(ptch%energy), ' ! 4 == Entrance, 5 == Exit'
    nl=nl+1; write (li(nl), '(2x, a, t40,  a)') '%time   [Time Patch at]:       ', patch_name(ptch%time)

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
! Subroutine set_ptc_com_pointers ()
!
! Routine to set ptc_com pointers to PTC global variables.
!-

subroutine set_ptc_com_pointers ()

use c_tpsa, only: EPS_EIGENVALUES_OFF_UNIT_CIRCLE
use ptc_spin, only: EXACT_MODEL, ALWAYS_EXACTMIS, VERTICAL_KICK, OLD_INTEGRATOR, HIGHEST_FRINGE

implicit none

logical init_needed

!

init_needed = (.not. associated(ptc_com%exact_model))

ptc_com%vertical_kick    => VERTICAL_KICK
ptc_com%old_integrator   => OLD_INTEGRATOR
ptc_com%exact_model      => EXACT_MODEL
ptc_com%exact_misalign   => ALWAYS_EXACTMIS
ptc_com%max_fringe_order => HIGHEST_FRINGE

if (init_needed) then
  ptc_com%exact_model = .true.
  ptc_com%exact_misalign = .true.  ! Note: Points to ALWAYS_EXACTMIS
  ptc_com%vertical_kick = 1        ! On
  ptc_com%old_integrator = -1 ! Using new integrator. -1 = False, 1 = True.
  EPS_EIGENVALUES_OFF_UNIT_CIRCLE = 1d-3

  ptc_com_default = ptc_com
  allocate (ptc_com_default%vertical_kick, ptc_com_default%old_integrator, ptc_com_default%exact_model, &
            ptc_com_default%exact_misalign, ptc_com_default%max_fringe_order)
  ptc_com_default%vertical_kick    = ptc_com%vertical_kick
  ptc_com_default%old_integrator   = ptc_com%old_integrator
  ptc_com_default%exact_model      = ptc_com%exact_model
  ptc_com_default%exact_misalign   = ptc_com%exact_misalign
  ptc_com_default%max_fringe_order = ptc_com%max_fringe_order
endif

end subroutine set_ptc_com_pointers

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine bmad_taylor_equal_damap (bmad_taylor, da)
!
! Routine to convert from PTC damap to Bmad Taylor map.
! This routine does not do any units conversion
!
! Input
!   da  -- damap: PTC damap.
!
! Output:
!   bmad_taylor(:) -- Taylor_struct: Input taylor map.
!-

subroutine bmad_taylor_equal_damap (bmad_taylor, da)

use s_fitting, only: alloc, kill, assignment(=), damap, real_8

implicit none

type (damap), intent(in) :: da
type (taylor_struct), intent(inout) :: bmad_taylor(:)
type (real_8) :: y8(size(bmad_taylor))

!

call alloc(y8)
y8 = da
bmad_taylor = y8
call kill (y8)

end subroutine bmad_taylor_equal_damap

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine damap_equal_bmad_taylor (da, bmad_taylor)
!
! Routine to convert from Bmad Taylor map to PTC damap.
! This routine does not do any units conversion
!
! Input
!   bmad_taylor(:) -- Taylor_struct: Input taylor map.
!
! Output:
!   da  -- damap: PTC damap.
!-

subroutine damap_equal_bmad_taylor (da, bmad_taylor)

use s_fitting, only: alloc, kill, assignment(=), damap, real_8

implicit none

type (damap), intent(inout) :: da
type (taylor_struct), intent(in) :: bmad_taylor(:)
type (real_8) :: y8(size(bmad_taylor))

!

call kill (da)
call alloc(da)

call alloc(y8)
y8 = bmad_taylor
da = y8
call kill (y8)

end subroutine damap_equal_bmad_taylor

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
! Subroutine form_complex_taylor (re_taylor, im_taylor, complex_taylor)
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
    ix1 = iter1()
    ix2 = iter2()
  else if (ix1 < ix2) then
    ! re term only
    ix1 = iter1()
  else
    ! im term only
    ix2 = iter2()
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
    ix1 = iter1()
    ix2 = iter2()
  
  else if (ix1 < ix2) then
    ! re term only
    expn = taylor1%term(t1)%expn
    re = taylor1%term(t1)%coef
    im = 0.0_rp
    ix1 = iter1()
  
  else
    ! im term only
    expn = taylor2%term(t2)%expn
    re = 0.0_rp
    im = taylor2%term(t2)%coef
    ix2 = iter2()
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

!-----------------------------------------
! Stepping helper routines
contains

function iter1() result (ix1)
implicit none
integer ix1
!
if (t1 < n1) then
  t1 = t1 + 1 
  ix1 = taylor_exponent_index(taylor1%term(t1)%expn)
else
  ix1 = huge(0) ! Set to largest integer
endif
end function iter1

!-----------------------------------------
! contains

function iter2() result (ix2)
implicit none
integer ix2
!
if (t2 < n2) then
  t2 = t2 + 1 
  ix2 = taylor_exponent_index(taylor2%term(t2)%expn)
else
  ix2 = huge(0)
endif
end function iter2

end subroutine form_complex_taylor

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
! In general Bmad treats a map y as being y(r-r_ref) where r 
! are the coordinates and r_ref is the reference coordinates at the beginning of y.
! The reference coordinates at the end of y is thus y(0).
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

! Allocate temp vars

call ptc_set_taylor_order_if_needed()

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

! Remove constant terms from the taylor map first. This is probably
! not needed but we do it to make sure everything is alright.
! Also remove terms that have higher order then bmad_com%taylor_order

call ptc_set_taylor_order_if_needed()

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
      if (sum(taylor_in(i)%term(j)%expn) > ptc_private%taylor_order_ptc) n = n - 1
    endif
  enddo

  if (associated(taylor_out(i)%term)) then
    if (size(taylor_out(i)%term) /= n) deallocate(taylor_out(i)%term)
  endif
  if (.not. associated(taylor_out(i)%term)) allocate (taylor_out(i)%term(n))

  nn = 0
  do j = 1, size(taylor_in(i)%term)
    ss = sum(taylor_in(i)%term(j)%expn)
    if (ss == 0 .or. (remove_higher_order_terms .and. ss > ptc_private%taylor_order_ptc)) cycle
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

character(*), parameter :: r_name = 'taylor_inverse'

!

call ptc_set_taylor_order_if_needed()

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

! Allocate temp vars

call ptc_set_taylor_order_if_needed()

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
! Subroutine concat_ele_taylor (taylor1, ele, taylor3, err_flag)
!
! Routine to concatinate two taylor maps:
!   taylor3[x] = ele_taylor(taylor1[x])
! If ele%taylor_map_includes_offsets = True:  ele_taylor == ele%taylor 
! If ele%taylor_map_includes_offsets = False: ele_taylor == ele%taylor + offset corrections. 
!
! Also see: concat_taylor
!
! Input:
!   taylor1(6)  -- Taylor_struct: Taylor map.
!   ele         -- ele_struct: Element containing a Taylor map.
!
! Output
!   taylor3(6)  -- Taylor_struct: Concatinated map
!   err_flag    -- logical: Set True if there is an error. False otherwise.
!-

Subroutine concat_ele_taylor (taylor1, ele, taylor3, err_flag)

use s_tracking, only: mis_fib, alloc, kill, dtiltd, assignment(=), real_8, fibre

implicit none

type (ele_struct) ele
type (taylor_struct) taylor1(:), taylor3(:)
type (lat_param_struct) param
type (real_8) x_ele(6), x_body(6), x1(6), x3(6)
type (fibre), pointer :: fib

real(rp) beta0, beta1, tilt
real(8) x_dp(6)
logical err_flag
character(*), parameter :: r_name = 'concat_ele_taylor'

! Match elements are not implemented in PTC so just use the matrix.

err_flag = .false.

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

call ptc_set_taylor_order_if_needed()

! Create a PTC fibre that holds the misalignment info
! and create map corresponding to ele%taylor.

param%particle = positron$  ! Actually this does not matter to the calculation
call ele_to_fibre (ele, fib, .true., err_flag)
if (err_flag) then
  call out_io(s_error$, r_name, 'CANNOT USE ELEMENT WITH PTC: ' // ele_full_name(ele))
  return
endif

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
call mis_fib (fib, x_ele, ptc_private%base_state, .true., entering = .true.)

x_body = ele%taylor

call concat_real_8 (x_ele, x_body, x_ele)
call mis_fib (fib, x_ele, ptc_private%base_state, .true., entering = .false.)
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
! The conversion includes the conversion between Bmad and PTC coordinates.
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
!   ref_orb_ptc(6)  -- real(rp), optional: PTC starting reference orbit.
!   exi_orb_ptc(6)  -- real(rp), optional: constant part of the map = orbit at the exit end.
!                        If present, the constant term of ptc_re8 will be removed.
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
call convert_bmad_to_ptc(ref_ptc, beta0, .true.) 

call alloc(diff_orb)
call alloc(start_orb)
call alloc(id)

id = 1
start_orb = id + ref_ptc
if (present(ref_orb_ptc)) ref_orb_ptc = ref_ptc

call convert_ptc_to_bmad(start_orb, beta0, .true.) 

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

call convert_bmad_to_ptc(ptc_re8, beta1, .true.)

! Remove constant.

if (present(exi_orb_ptc)) then
  exi_orb_ptc = ptc_re8
  do i = 1, 6
    ptc_re8(i) = ptc_re8(i) - exi_orb_ptc(i)
  enddo
endif

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
! Subroutine real_8_to_taylor (y8, beta0, beta1, bmad_taylor)
!
! Routine to convert a PTC real_8 map to a Bmad Taylor map.
! The conversion includes the conversion from PTC to Bmad coordinate systems.
! Warning: This routine has never been used and needs to be tested!
!
! Input:
!   y8(6)           -- real_8: PTC Taylor map. NOTE: y8 is used as scratch space and therefore trashed.
!   beta0           -- real(rp): Reference particle velocity at beginning of map
!   beta1           -- real(rp): Reference particle velocity at end of map
!   bmad_taylor(6)  -- Taylor_struct: Only %ref is used at input.
!     %ref            -- Reference orbit
!
! Output:
!   bmad_taylor(6) -- Taylor_struct: Bmad Taylor map.
!-

subroutine real_8_to_taylor (y8, beta0, beta1, bmad_taylor)

use s_fibre_bundle

implicit none

type (taylor_struct) :: bmad_taylor(:)
type (real_8) y8(:), rr(6), bet, ss(6)
type (damap) bm, id, si

real(rp) beta0, beta1, fix0(6)

!

call alloc (bm, id, si)
call alloc (rr)
call alloc (bet)
call alloc (ss)

bm = y8

fix0 = bm
id = 1

rr = id + bmad_taylor%ref

ss = rr 
ss(5) = (rr(6)**2+2.d0*rr(6))/(1.d0/beta0 + sqrt(1.d0/beta0**2+rr(6)**2+2.d0*rr(6)))
bet = (1.d0+rr(6))/(1.d0/beta0+ss(5))
ss(6) = -rr(5)/bet

si=ss  ! bmad to ptc map

bm = bm * si
bm = fix0

rr = bm
ss = rr
ss(6) = (2.d0*rr(5)/beta1+rr(5)**2)/(sqrt(1.d0+2.d0*rr(5)/beta1+rr(5)**2)+1.d0)
bet = (1.d0+ss(6))/(1.d0/beta1+rr(5))
ss(5) = -bet*rr(6)

bmad_taylor = ss

call kill (rr)
call kill (ss)
call kill (bet)
call kill (bm, id, si)

end subroutine real_8_to_taylor

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_propagate1 (bmad_taylor, ele, param, err_flag, ref_in)
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
!   ref_in           -- coord_struct, optional: Particle to be tracked.
!                         Must be present if the particle to be tracked is not the reference particle or
!                         if the direction of propagation is backwards.
!
! Output:
!   bmad_taylor(6)  -- Taylor_struct: Map through element.
!   err_flag        -- logical: Set True if there is an error. False otherwise.
!-

subroutine taylor_propagate1 (bmad_taylor, ele, param, err_flag, ref_in)

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
type (coord_struct), optional :: ref_in

real(rp) beta0, beta1, m2_rel
logical err_flag
character(*), parameter :: r_name = 'taylor_propagate1'

! If the element is a taylor then just concat since this is faster.

if (ele%key == taylor$) then
  call concat_ele_taylor (bmad_taylor, ele, bmad_taylor, err_flag)
  return
endif

! set the taylor order in PTC if not already done so

call ptc_set_taylor_order_if_needed()

! Init ptc map with bmad map

beta0 = ele%value(p0c_start$) / ele%value(e_tot_start$)
beta1 = ele%value(p0c$) / ele%value(e_tot$)

call alloc (ptc_tlr)
ptc_tlr = bmad_taylor

! track the map

call ele_to_fibre (ele, ptc_fibre, .true., err_flag, ref_in = ref_in)
if (err_flag) then
  call out_io(s_error$, r_name, 'CANNOT USE ELEMENT WITH PTC: ' // ele_full_name(ele))
  return
endif

call track_probe_x (ptc_tlr, ptc_private%base_state, fibre1 = bmadl%start)

! transfer ptc map back to bmad map

bmad_taylor = ptc_tlr

! cleanup

call kill (ptc_tlr)

end subroutine taylor_propagate1

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

call indexer (ord, ix)

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
! Subroutine beambeam_fibre_setup(ele, ptc_fibre)
!
! Routine to setup a fibre to handle the beambeam interaction.
!
! Input:
!   ele       -- ele_struct: Bmad beambeam element.
!
! Output:
!    ptc_fibre  -- Corresponding PTC fibre.
!-

subroutine beambeam_fibre_setup (ele, ptc_fibre)

use madx_ptc_module, only: integration_node, alloc, DO_BEAM_BEAM

implicit none

type (ele_struct), target :: ele
type (fibre), target :: ptc_fibre
type (integration_node), pointer :: node

real(rp) sig_x0, sig_y0, beta_a0, beta_b0, alpha_a0, alpha_b0, sig_x, sig_y
real(rp) alpha, beta, s, bbi_const, r
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
  s = z_slice(i) / 2  ! Factor of 2 since strong beam is moving.

  if (beta_a0 == 0 .or. beta_b0 == 0) then
    sig_x = sig_x0
    sig_y = sig_y0
  else
    beta = beta_a0 - 2 * alpha_a0 * s + (1 + alpha_a0**2) * s**2 / beta_a0
    sig_x = sig_x0 * sqrt(beta / beta_a0)
    beta = beta_b0 - 2 * alpha_b0 * s + (1 + alpha_b0**2) * s**2 / beta_b0
    sig_y = sig_y0 * sqrt(beta / beta_b0)
  endif

  bbi_const = -2.0_rp * strong_beam_strength(ele) * classical_radius_factor / ele%value(e_tot$)

  node%bb%bbk(i,:) = 0  ! MAD closed orbit kick. Not used here.

  r = ((((ele%value(crab_x5$) * s + ele%value(crab_x4$)) * s + ele%value(crab_x3$)) * s + & 
                                    ele%value(crab_x2$)) * s + ele%value(crab_x1$)) * s
  node%bb%xm(i) = r * cos(ele%value(crab_tilt$))
  node%bb%ym(i) = r * sin(ele%value(crab_tilt$))

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
! Subroutine misalign_ptc_fibre (ele, use_offsets, ptc_fibre, for_layout)
!
! Routine to misalign a fibre associated with a Bmad element.
!
! Input:
!   ele         -- ele_struct: Bmad element with misalignments.
!   use_offsets -- logical: Does ptc_fibre include element offsets, pitches and tilt?
!                   This argument is ignored if the element is a patch.
!   for_layout  -- logical: If True then fibre is being created as part of a layout as
!                   opposed to a stand-alone fibre
!
! Output:
!   ptc_fibre   -- fibre: PTC fibre element with misalignments.
!-

subroutine misalign_ptc_fibre (ele, use_offsets, ptc_fibre, for_layout)

use madx_ptc_module
use s_frame, only: find_patch_bmad_marker

implicit none

type (ele_struct) ele
type (floor_position_struct) :: floor0, floor1
type (fibre), pointer :: ptc_fibre
type (fibre), target :: dummy_fibre

real(rp) dr(3), ang(3), exi(3,3), beta_start, beta_end
real(rp) x(6), o_chord(3), o_arc(3), basis(3,3), orient(3,3), sagitta
real(rp) r0(3), s_mat(3,3), s0_mat(3,3), r(3), b(3), t, rot(2,2)

logical use_offsets, for_layout, good_patch, addin

character(*), parameter :: r_name = 'misalign_ptc_fibre'

!

if (ele%key == fiducial$) return

if (ele%key == floor_shift$) then
  r0 = 0
  call mat_make_unit(s0_mat)
  dr = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)]
  call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat_inv = s_mat)
  call find_patch(r0, s0_mat, dr, s_mat, r, ang)
  ptc_fibre%patch%b_ang = ang
  ptc_fibre%patch%b_d   = r
  ptc_fibre%patch%patch = 2
  return
endif

! In ptc there is no such thing as a reversed patch. Therefore need to
! use the reverse transformation if the patch is reversed in Bmad.

! Also fibre%dir for a patch must agree with the preceeding element in a layout

if (ele%key == patch$) then
  if (ele%orientation == -1) then
    floor0 = floor_position_struct(vec3_zero$, mat3_unit$, 0.0_rp, 0.0_rp, 0.0_rp)
    call ele_geometry (floor0, ele, floor1)
    dr = floor1%r
    ang = [floor1%phi, floor1%theta, floor1%psi]
  else
    dr = [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)]
    ang = [ele%value(y_pitch$), ele%value(x_pitch$), ele%value(tilt$)]
  endif

  call bmad_patch_parameters_to_ptc (ang, exi)

  ptc_fibre%dir = nint(ele%value(upstream_coord_dir$))

  if (for_layout) then
    call survey_integration_fibre(ptc_fibre%next, exi = transpose(ele%floor%w), b = ele%floor%r)
    call find_patch(ptc_fibre, patching = good_patch)
    call survey(ptc_fibre%next, ptc_fibre%parent_layout%end)
  else
    call find_patch_bmad_marker(ptc_fibre, dr, exi, nint(ele%value(downstream_coord_dir$)), patching = good_patch) 
  endif

  ! energy and time patches

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
  ! For Etienne this can be an annoyance so he sets ptc_com%translate_patch_drift_time = False.

  if (ptc_com%translate_patch_drift_time) then
    ptc_fibre%patch%time = 2     ! Subtract off reference time (which affects z in tracking).
    ptc_fibre%patch%b_t = ele%value(l$) / beta_end + ele%orientation * ele%value(t_offset$) * c_light 
    ptc_fibre%patch%b_l = ptc_fibre%patch%b_t * beta_end
  endif

!----------------------------------------------------------------------
! Not patch nor floor_shift element.

elseif (use_offsets .and. ele%value(dispatch$) /= no_misalignment$) then

  ! Patch elements do not have misalignments

  if (ele%key == fiducial$) return
  if (attribute_index(ele, 'X_OFFSET_TOT') < 1) return

  !

  addin = .false.

  if (ele%key == sbend$) then
    sagitta = ele%value(l_sagitta$)
    orient = ptc_fibre%chart%f%ent
    if (ele%value(ref_tilt$) == 0) then
      call geo_rot(orient, [0.0_rp, ele%value(angle$)/2.0_rp, 0.0_rp], 1, orient)
    else
      call geo_rot(orient, [0.0_rp, 0.0_rp, ele%value(ref_tilt$)], 1, orient)
      call geo_rot(orient, [0.0_rp, ele%value(angle$)/2.0_rp, 0.0_rp], 1, orient)
      call geo_rot(orient, [0.0_rp, 0.0_rp, -ele%value(ref_tilt$)], 1, orient)
    endif

    dr = [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)]

    if (any(dr /= 0)) then
      o_chord = ptc_fibre%mag%p%f%o
      o_arc = o_chord + sagitta * ptc_fibre%mag%p%f%mid(1,1:3)
      basis = orient
      x = 0
      x(1:3) = dr
      call misalign_fibre (ptc_fibre, x, o_arc, basis, add = addin)
      addin = .true.
    endif

    if (ele%value(x_pitch_tot$) /= 0) then
      o_chord = ptc_fibre%mag%p%f%o
      o_arc = o_chord + sagitta * ptc_fibre%mag%p%f%mid(1,1:3)
      basis = orient
      x = 0
      x(5) = -ele%value(x_pitch_tot$)
      call misalign_fibre (ptc_fibre, x, o_arc, basis, add = addin)
      addin = .true.
      call geo_rot(orient, x(4:6), 1, orient)
    endif

    if (ele%value(y_pitch_tot$) /= 0) then
      o_chord = ptc_fibre%mag%p%f%o
      o_arc = o_chord + sagitta * ptc_fibre%mag%p%f%mid(1,1:3)
      basis = orient
      x = 0
      x(4) = -ele%value(y_pitch_tot$)
      call misalign_fibre (ptc_fibre, x, o_arc, basis, add = addin)
      addin = .true.
    endif

    if (ele%value(roll_tot$) /= 0) then
      o_chord = ptc_fibre%mag%p%f%o   ! Chord origin
      basis = orient
      x = 0
      x(6) = ele%value(roll_tot$)
      call misalign_fibre (ptc_fibre, x, o_chord, basis, add = addin)
      addin = .true.
      call geo_rot(orient, x(4:6), 1, orient)
    endif

  else ! Not a bend
    dr = [ele%value(x_offset_tot$), ele%value(y_offset_tot$), ele%value(z_offset_tot$)]
    if (ele%key == sad_mult$ .and. ele%value(l$) == 0) &
                              dr(1:2) = dr(1:2) + [ele%value(x_offset_mult$), ele%value(y_offset_mult$)]

    if (any(dr /= 0)) then
      o_chord = ptc_fibre%mag%p%f%o
      basis = ptc_fibre%mag%p%f%mid
      x = 0
      x(1:3) = dr
      call misalign_fibre (ptc_fibre, x, o_chord, basis, add = addin)
      addin = .true.
    endif

    if (ele%value(y_pitch_tot$) /= 0) then
      o_chord = ptc_fibre%mag%p%f%o
      basis = ptc_fibre%mag%p%f%mid
      x = 0
      x(4) = -ele%value(y_pitch_tot$)
      call misalign_fibre (ptc_fibre, x, o_chord, basis, add = addin)
      addin = .true.
    endif

    if (ele%value(x_pitch_tot$) /= 0) then
      o_chord = ptc_fibre%mag%p%f%o
      basis = ptc_fibre%mag%p%f%mid
      x = 0
      x(5) = -ele%value(x_pitch_tot$)
      call misalign_fibre (ptc_fibre, x, o_chord, basis, add = addin)
      addin = .true.
    endif

    if (ele%value(tilt_tot$) /= 0) then
      o_chord = ptc_fibre%mag%p%f%o   ! Chord origin
      basis = ptc_fibre%mag%p%f%mid
      x = 0
      x(6) = ele%value(tilt_tot$)
      call misalign_fibre (ptc_fibre, x, o_chord, basis, add = addin)
      addin = .true.
    endif

  endif

  !

  if (ele%value(e_tot_start$) /= ele%value(e_tot$)) then
    ptc_fibre%patch%energy = 4   ! Entrance energy patch
    ptc_fibre%patch%p0b = ele%value(p0c_start$) * 1d-9
    ptc_fibre%patch%b0b = ele%value(p0c_start$) / ele%value(E_tot_start$)  ! beta velocity
  endif

  ! Put misalignmnets in patches?

  if (.not. ptc_com%use_orientation_patches) then
    call convert_mis_to_patch(ptc_fibre, .true.)
  endif

endif

end subroutine misalign_ptc_fibre

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
! Subroutine ele_to_ptc_magnetic_bn_an (ele, bn, an, n_max)
!
! Routine to compute the a(n) and b(n) magnetic multipole components of a magnet.
! This is used to interface between eles and PTC fibres
!
! Note: The multipole index uses the PTC convention of starting from 1 instead of zero.
!
! Note: On the PTC side bn(1) is error field when creating a fibre but 
! is the total field when the fibre is being modified. This routine returns the error field.
!
! Input:
!   ele                 -- ele_struct: Bmad Element.
!
! Output:
!   bn(1:n_pole_maxx+1) -- real(rp): Normal multipole component.
!   an(1:n_pole_maxx+1) -- real(rp): Skew multipole component.
!   n_max               -- integer, optional: Maximum non-zero multipole component.
!                           Set to zero if there are no multipoles.
!-

subroutine ele_to_ptc_magnetic_bn_an (ele, bn, an, n_max)

implicit none

type (ele_struct), target :: ele

real(rp) bn(:), an(:)
real(rp) cos_t, sin_t, leng, hk, vk, tilt
real(rp), pointer :: val(:)
real(rp), target, save :: value0(num_ele_attrib$) = 0
real(rp) an0(0:n_pole_maxx), bn0(0:n_pole_maxx)

integer, optional :: n_max
integer n, key, ix_pole_max
logical kick_here, add_kick, add_multipoles

character(*), parameter :: r_name = 'ele_to_ptc_magnetic_bn_an'

!

if (present(n_max)) n_max = 0
bn = 0
an = 0

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

case (marker$, detector$, fork$, photon_fork$, beginning_ele$, em_field$, patch$, fiducial$, floor_shift$, gkicker$)
  return

case (crab_cavity$)
  if (leng == 0) then
    bn(1) = 1d-9 * e_accel_field(ele, voltage$)
  else
    bn(1) = 1d-9 * e_accel_field(ele, voltage$) / leng
  endif

case (rfcavity$, lcavity$, drift$, ab_multipole$, multipole$, beambeam$, wiggler$, undulator$, thick_multipole$)
  ! Nothing to be done

case (octupole$)
  bn(4) = val(k3$) / 6

case (quadrupole$) 
  bn(2) = val(k1$)

case (sad_mult$)

case (sbend$)
  if (ele%is_on) then

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
      bn(1) = bn(1) + ele%value(dg$)
      bn(2) = val(k1$)
      bn(3) = val(k2$) / 2
    endif

  else
    bn(1) = -ele%value(g$)
  endif

case (sextupole$)
  bn(3) = val(k2$) / 2

case (rcollimator$, ecollimator$, monitor$, instrument$, pipe$, ac_kicker$)

case (solenoid$)

case (sol_quad$)
  bn(2) = val(k1$)

case (hkicker$, vkicker$)
  if (ele%key == hkicker$) bn(1) = val(kick$)   ! notice that kickers do not have kick scaled by the length
  if (ele%key == vkicker$) an(1) = val(kick$) 
  add_kick = .false.

case (kicker$)
  bn(1)  = val(hkick$)    ! notice that kicker does not have kick scaled by the length
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

if (add_multipoles) then
  call multipole_ele_to_ab (ele, .false., ix_pole_max, an0, bn0)

  select case (ele%key)
  case (hkicker$, vkicker$, kicker$)
  case (lcavity$, rfcavity$)
    an0 = an0 / ele%value(l_active$)
    bn0 = bn0 / ele%value(l_active$)
  case default
    if (leng /= 0) then
      an0 = an0 / leng
      bn0 = bn0 / leng
    endif
  end select

  n = n_pole_maxx+1
  an(1:n) = an(1:n) + an0(0:n-1)
  bn(1:n) = bn(1:n) + bn0(0:n-1)
endif

if (ele%key == ac_kicker$) then
  an = 1e-9_rp * ele%value(p0c$) * an * ele%ac_kick%frequency(1)%amp
  bn = 1e-9_rp * ele%value(p0c$) * bn * ele%ac_kick%frequency(1)%amp
endif

if (present(n_max)) then
  do n = size(bn), 1, -1
    if (an(n) /= 0 .or. bn(n) /= 0) exit
  enddo
  n_max  = n
endif

end subroutine ele_to_ptc_magnetic_bn_an

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
! Note: Only suppress printing if ptc_com%print_info_messages = F.
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

if (ptc_com%print_info_messages) return

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
