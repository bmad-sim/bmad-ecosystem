#include "CESR_platform.inc"

module ptc_interface_mod


  use bmad_struct
  use bmad_interface
  use multipole_mod

  use mad_like, only: kill, set_up, real_8, layout, fibre, &
        universal_taylor, dp, ptc_track => track, append, ring_l

  interface assignment (=)
    module procedure real_8_equal_taylor
    module procedure taylor_equal_real_8
    module procedure universal_equal_universal
  end interface

  integer, parameter :: bmad_std$ = 1, ptc_std$ = -1

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Function map_coef(y, i, j, k, l, style)
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
!   use bmad
!
! Input:
!   y(:)  -- Real_8: Taylor Map.
!   i     -- Integer: output index.
!   j     -- Integer, optional: 1st input index, needed for 1st order and above.
!   k     -- Integer, optional: 2nd input index, needed for 2nd order and above.
!   l     -- Integer, optional: 3rd input index, needed for 3rd order.
!   style -- Integer, optional: bmad_std$ or ptc_std$.
!             bmad_std$ --> (x, p_x, y, p_y, z, dE/E) for phase space vars.
!             ptc_std$  --> (x, p_x, y, p_y, dE/E, ct ~ -z)
!             default = ptc_std$.
!
! Output:
!   map_coef -- Real*8: Coefficient.
!-

function map_coef (y, i, j, k, l, style) result (m_coef)

  use polymorphic_taylor

  implicit none

  type (real_8) y(:)


  real*8 m_coef

  integer i
  integer, optional :: j, k, l, style
  integer arr(40), n_max, sgn, ii

  character str*40
  character, parameter :: str1(4) = (/ '1', '2', '3', '4' /)

  logical use_bmad

!

  use_bmad = .false.
  if (present(style)) then
    if (style == bmad_std$) use_bmad = .true.
  endif

  arr = 0
  sgn = 1
  str = '0000000000000000000000000000000000000000'
  n_max = 1

  call map_index(j)
  call map_index(k)
  call map_index(l)

  call map_index(i, ii)
  m_coef = sgn * y(ii)%t.sub.str(1:n_max)

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

end subroutine

end function

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine type_layout (lay)
!
! Subroutine to print the global information in a layout
!
! Modules Needed:
!   use bmad
!
! Input:
!   lay - layout: layout to use.
!+

subroutine type_layout (lay)

  implicit none

  type (layout) lay

!

  if (.not. associated(lay%start)) then
    print *, 'Warning from TYPE_LAYOUT: Layout NOT Associated'
    return
  endif

  print *, 'Name:         ', lay%name
  print *, 'N:            ', lay%N,        '  ! Number of Elements'
  print *, 'LatPos:       ', lay%lastpos,  '  ! Last position'

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ring_to_layout (ring, ptc_layout)
!
! Subroutine to create a PTC layout from a BMAD ring.
! Note: If ptc_layout has been already used then you should first do a 
!           call kill(ptc_layout)
! This deallocates the pointers in the layout
!
! Note: Before you call this routine you need to first call:
!    call set_ptc (ring%param, ...)
!
! Modules needed:
!   use bmad
!
! Input:
!   ring -- Ring_struct: 
!
! Output:
!   ptc_layout -- Layout:
!-

subroutine ring_to_layout (ring, ptc_layout)

  implicit none

  type (ring_struct), intent(in) :: ring
  type (layout), intent(inout) :: ptc_layout
  type (fibre), pointer :: fib

  real(8) energy, kinetic, beta0, p0c, brho

  integer i

! setup

  call set_up (ptc_layout)

! transfer energy, etc.

!  ptc_layout%energy = 1e-9 * ring%param%beam_energy
!  call energy_to_kinetic (ring_param%beam_energy, ring%param%particle, &
!                                                 kinetic, beta0, p0c, brho)
!  ptc_layout%kinetic = kinetic
!  ptc_layout%beta0, = beta0
!  ptc_layout%p0c = p0c
!  ptc_layout%brho = brho
!  
!  ptc_layout%circumference = 0


! transfer elements.

  do i = 1, ring%n_ele_ring
    allocate (fib)
    call ele_to_fibre (ring%ele_(i), fib, ring%param)
    call append (ptc_layout, fib)
    call kill (fib)
  enddo

! circular or not?

  if (ring%param%lattice_type == circular_lattice$) then
    ptc_layout%closed = .true.
    call ring_l (ptc_layout, .true.)
  else
    ptc_layout%closed = .false.
  endif

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine type_map1 (y, type0, n_dim, style)
!
! Subroutine to type the transfer map up to first order.
!
! Modules needed:
!   use bmad
!
! Input:
!   y(:)  -- Real_8: 
!   type0 -- Logical: Type zeroth order map
!   n_dim -- Integer: Number of dimensions to type: 4 or 6
!   style -- Integer, optional: bmad_std$ or ptc_std$.
!             bmad_std$ --> (x, p_x, y, p_y, z, dE/E) for phase space vars.
!             ptc_std$  --> (x, p_x, y, p_y, dE/E, ct ~ -z) 
!             default = ptc_std$.
!-

subroutine type_map1 (y, type0, n_dim, style)

  implicit none

  type (real_8), intent(in) :: y(:)

  integer, intent(in) :: n_dim
  integer, optional, intent(in) :: style
  integer :: i, j

  logical, intent(in) :: type0

!

  if (type0) then
    print *, '0th Order Map:'
    print '(6f11.5)', (map_coef(y(:), i, style=style), i = 1, n_dim)
    print *
  endif

  print *, '1st Order Map:'
  do i = 1, n_dim
    print '(6f11.5)', (map_coef(y(:), i, j, style=style), j = 1, n_dim)
  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! function kind_name (this_kind)
!
! function to return the name of a PTC kind.
!
! Input:
!   this_kind -- Integer: 
!
! Output:
!   kind_name -- Character*20: 
!-

function kind_name (this_kind)

  use s_status, only: kind0, kind1, kind2, kind3, kind4, kind5, kind6, kind7, &
        kind8, kind9, kind10, kindfitted, kinduser1, kinduser2

  implicit none

  integer this_kind
  character*20 kind_name

!

  select case (this_kind)
  case (kind0); kind_name  = 'KIND0'
  case (kind1); kind_name  = 'DRIFT1'
  case (kind2); kind_name  = 'DKD2 (Gen Element)' 
  case (kind3); kind_name  = 'KICKT3 (Thin Ele)'
  case (kind4); kind_name  = 'CAV4 (RF Cavity)'
  case (kind5); kind_name  = 'SOL5 (Solenoid)'
  case (kind6); kind_name  = 'KTK (Slow Thick)'
  case (kind7); kind_name  = 'TKTF (Fast Thick)'
  case (kind8); kind_name  = 'NSMI (Normal SMI)'
  case (kind9); kind_name  = 'SSMI (Skew SMI)'
  case (kind10); kind_name = 'TEAPOT (Sector Bend)'
  case (kindfitted); kind_name = 'FITTED KIND'
  case (kinduser1); kind_name = 'USER1 KIND'
  case (kinduser2); kind_name = 'USER2 KIND'
  case default; write (kind_name, '(a, i5)') 'UNKNOWN KIND!', this_kind 
  end select

end function

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine type_fibre (fib)
!
! Subroutine to print the global information in a fibre
!
! Modules Needed:
!   use bmad
!
! Input:
!   fib - fibre: fibre to use.
!-

subroutine type_fibre (fib)

  use s_status

  implicit none

  type (fibre), intent(in) :: fib

  integer i

!
  
  if (.not. associated (fib%mag)) then
    print *, 'Warning from TYPE_FIBRE: Fibre NOT associated with anything.'
    return
  endif

  print *, 'Name:        ', fib%mag%name
  print *, 'Vorname:     ', fib%mag%vorname
  print *, 'Kind:        ', kind_name(fib%mag%kind)
  print *, 'Knob:        ', fib%magp%knob
    print *, 'L:        ', fib%mag%l

  if (fib%mag%kind == kind4) then
    print *, 'Voltage:  ', fib%mag%volt
    print *, 'Frequency:', fib%mag%freq
    print *, 'Voltage:  ', fib%mag%volt
    print *, 'Phase:    ', fib%mag%phas
    print *, 'Delta_e:  ', fib%mag%delta_e
    print *, 'Thin: ', fib%mag%thin
  endif

  if (fib%mag%kind == kind5) then
    print *, 'KS:       ', fib%mag%b_sol
    print *, 'Thin:     ', fib%mag%thin
  endif

  if (fib%mag%kind == kind2 .and. fib%mag%p%b0 /= 0) then
    print *, 'E1:       ', fib%mag%p%edge(1)
    print *, 'E2:       ', fib%mag%p%edge(2)
    print *, 'Rho:      ', fib%mag%p%b0
    print *, 'L_chord:  ', fib%mag%p%lc
  endif

  print *, 'Integration Order: ', fib%mag%p%method
  print *, 'Integration Steps: ', fib%mag%p%nst


  do i = lbound(fib%mag%bn, 1), ubound(fib%mag%bn, 1)
    if (fib%mag%bn(i) /= 0) print '(a, i2, a, 5x, 1pd12.3)', &
                                  ' BN(', i, '):', fib%mag%bn(i)
  enddo  
  do i = lbound(fib%mag%an, 1), ubound(fib%mag%an, 1)
    if (fib%mag%an(i) /= 0) print '(a, i2, a, 5x, 1pd12.3)', &
                                  ' AN(', i, '):', fib%mag%an(i)
  enddo  


end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine set_ptc (param, taylor_order, integ_order, &
!                               num_steps, no_cavity, exact_calc)
!
! Subroutine to initialize ptc.
! Note: This subroutine cannot be used if there are "knobs".
! This subroutine replaces:
!     make_states
!     set_mad
!     init
!
! Modules needed:
!   use bmad
!
! Input:
!   param        -- Param_struct, optional: BMAD parameters:
!     %beam_energy -- Energy.
!     %particle    -- Type of particle.
!   taylor_order -- Integer, optional: Maximum order of the taylor polynomials.
!   integ_order  -- Integer, optional: Default Order for the drift-kick-drift 
!                     sympletic integrator. Possibilities are: 2, 4, or 6
!                     Default = 2
!   num_steps    -- Integer, optional: Default Number of integration steps.
!                     Default = 1
!   no_cavity    -- Logical, optional: RF Cavity exists? 
!                     Default = False.
!   exact_calc   -- logical, optional: Exact Model? 
!                     Default = False.
!-

subroutine set_ptc (param, taylor_order, integ_order, &
                                    num_steps, no_cavity, exact_calc) 

  use mad_like

  implicit none

  type (param_struct), optional :: param

  integer, optional :: integ_order, num_steps, taylor_order
  integer this_method, this_steps
  integer nd2, npara

  real(dp) this_energy

  logical, optional :: no_cavity, exact_calc
  logical, save :: init_needed = .true.
              
! do not call set_mad

  if (init_needed .and. present(param)) then
    if (param%particle == positron$ .or. param%particle == electron$) then
      call make_states(.true.)
    else
      call make_states(.false.)
    endif
    EXACT_MODEL = .false.
    ALWAYS_EXACTMIS = .false.
    init_needed = .false.
  endif

  if (present(exact_calc)) then
    EXACT_MODEL = exact_calc
    ALWAYS_EXACTMIS = exact_calc
  endif
    
  if (present(no_cavity)) default = default+nocavity
  
  if (present (integ_order)) then
    this_method = integ_order
    bmad_com%default_integ_order = integ_order
  else
    this_method = bmad_com%default_integ_order
  endif

  if (present (num_steps)) then
    this_steps = num_steps
    bmad_com%default_num_steps = num_steps
  else
    this_steps = bmad_com%default_num_steps
  endif

  if (present(param)) then
    if (bmad_com%beam_energy /= param%beam_energy .or. &
                        present(integ_order) .or. present(num_steps)) then
      this_energy = 1e-9 * param%beam_energy
      call set_mad (energy = this_energy, method = this_method, &
                                                       step = this_steps)
      bmad_com%beam_energy  = param%beam_energy
    endif
  endif

  if (present(taylor_order)) then  
    if (bmad_com%taylor_order_ptc /= taylor_order) then
      call init (default, taylor_order, 0, berz, nd2, &
                                               bmad_com%real_8_map_init)
      bmad_com%taylor_order_ptc = taylor_order
      bmad_com%taylor_order     = taylor_order
    endif
  endif
  
end subroutine  

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine real_8_equal_taylor (y8, bmad_taylor)
!
! Subroutine to overload "=" in expressions
!       y8 = bmad_taylor
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor(:) -- Taylor_struct: Input taylor series array.
!
! Output:
!   y8(:) -- real_8: PTC Taylor series array.
!-

subroutine real_8_equal_taylor (y8, bmad_taylor)

  implicit none

  type (real_8), intent(out) :: y8(:)
  type (taylor_struct), intent(in) :: bmad_taylor(:)

  call taylor_to_real_8 (bmad_taylor, y8, .true.)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_equal_real_8 (bmad_taylor, y8)
!
! Subroutine to overload "=" in expressions
!       bmad_taylor = y8
!
! Modules needed:
!   use bmad
!
! Input:
!   y8(:) -- real_8: PTC Taylor series array.
!
! Output:
!   bmad_taylor(:) -- Taylor_struct: Input taylor series array.
!-

subroutine taylor_equal_real_8 (bmad_taylor, y8)

  implicit none

  type (real_8), intent(in) :: y8(:)
  type (taylor_struct), intent(out) :: bmad_taylor(:)

  call real_8_to_taylor (y8, bmad_taylor, .true.)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine real_8_to_taylor (y8, bmad_taylor, switch_z)
!
! Subroutine to convert from a real_8 taylor map in Etienne's PTC 
! to a taylor map in BMAD.
! The conversion can also convert from the PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
! to the BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
!
! Modules needed:
!   use bmad
!
! Input:
!   y8(6)       -- Real_8: Taylor map.
!   switch_z    -- Logical, optional: If True then switch coordinate 
!                    conventions. Default is True.
!
! Output:
!   bmad_taylor(6) -- Taylor_struct:
!-

subroutine real_8_to_taylor (y8, bmad_taylor, switch_z)

  use polymorphic_taylor

  implicit none

  type (real_8), intent(in) :: y8(:)
  type (taylor_struct), intent(out) :: bmad_taylor(:)
  type (universal_taylor) :: u_t(6)

  integer i

  logical, optional :: switch_z

!

  do i = 1, 6
    u_t(i) = 0  ! nullify
    u_t(i) = y8(i)%t
  enddo

  call universal_to_bmad_taylor (u_t, bmad_taylor, switch_z)

  do i = 1, 6
    u_t(i) = -1  ! deallocate
  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_to_real_8 (bmad_taylor, y8, switch_z)
!
! Subroutine to convert from a taylor map in BMAD to a
! real_8 taylor map in Etienne's PTC.
! The conversion can also convert from the the BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
! to PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
!
! Modules needed:
!   use bmad
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Taylor map.
!   switch_z       -- Logical, optional: If True then switch coordinate 
!                       conventions. Default is True.
!
! Output:
!   y8(6)       -- Real_8: Taylor map.
!-

subroutine taylor_to_real_8 (bmad_taylor, y8, switch_z)

  use polymorphic_taylor

  implicit none

  type (real_8), intent(out) :: y8(:)
  type (taylor_struct), intent(in) :: bmad_taylor(:)
  type (universal_taylor) :: u_t

  integer i, j, ii, n

  logical, optional :: switch_z
  logical switch

! init

  call kill (y8)
  call real_8_init (y8, .true.)

!

  do i = 1, 6

    switch = .true.
    if (present(switch_z)) switch = switch_z

    ii = i
    if (switch) then
      if (i == 5) ii = 6
      if (i == 6) ii = 5
    endif

    n = size(bmad_taylor(i)%term)
    allocate (u_t%n, u_t%nv, u_t%c(n), u_t%j(n,6))
    u_t%n = n
    u_t%nv = 6

    do j = 1, n
      if (switch) then
        u_t%j(j,:) = bmad_taylor(i)%term(j)%exp((/1,2,3,4,6,5/))
        u_t%c(j) = bmad_taylor(i)%term(j)%coef * &
                                      (-1)**bmad_taylor(i)%term(j)%exp(5)
        if (i == 5) u_t%c(j) = -u_t%c(j)
      else
        u_t%j(j,:) = bmad_taylor(i)%term(j)%exp(:)
        u_t%c(j) = bmad_taylor(i)%term(j)%coef
      endif
    enddo

    y8(ii) = u_t
    u_t = -1   ! deallocate
        
  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine vec_bmad_to_ptc (vec_bmad, vec_ptc)
!
! Subroutine to convert between BMAD and PTC coordinates.
! PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
! BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
!
! Input:
!   vec_bmad(6) -- Real(rp): Input BMAD vector.
!
! Output:
!   vec_ptc(6)  -- Real(dp): Output PTC vector.
!-

subroutine vec_bmad_to_ptc (vec_bmad, vec_ptc)

  implicit none
  
  real(rp), intent(in)  :: vec_bmad(:)
  real(dp), intent(out)   :: vec_ptc(:)
  real(dp) temp_vec(6)

  temp_vec = vec_bmad((/1,2,3,4,6,5/))
  vec_ptc = temp_vec
  vec_ptc(6) = -vec_ptc(6)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine vec_ptc_to_bmad (vec_ptc, vec_bmad)
!
! Subroutine to convert between BMAD and PTC coordinates.
! PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
! BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
!
! Input:
!   vec_ptc(6)  -- Real(rp): Input PTC vector.
!
! Output:
!   vec_bmad(6) -- Real(rp): Output BMAD vector.
!-

subroutine vec_ptc_to_bmad (vec_ptc, vec_bmad)

  implicit none
  
  real(dp), intent(in)    :: vec_ptc(:)
  real(rp), intent(out) :: vec_bmad(:)
  real(rp) temp(6)

  temp = vec_ptc((/1,2,3,4,6,5/))
  vec_bmad = temp
  vec_bmad(5) = -vec_bmad(5)

end subroutine

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
!   use bmad
!
! Input:
!   ut2 -- Universal_taylor:
!
! Output:
!   ut1 -- Universal_taylor:
!-

elemental subroutine universal_equal_universal (ut1, ut2)

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
! Note: If this variable has been used before make sure you have 
! deallocated using:
!   call kill(y)
!
! Modules needed:
!   use bmad
!
! Input:
!   y(:)       -- Real_8: 
!   set_taylor -- Logical, optional :: If present and True then make
!                   y the identity taylor series (kind = 2).
!
! Output:
!   y(:) -- Real_8:
!-

subroutine real_8_init (y, set_taylor)

  use s_fibre_bundle

  implicit none
  
  type (real_8) :: y(:)
  real(dp) :: x(6) = (/ 0, 0, 0, 0, 0, 0 /)

  logical, optional :: set_taylor

!

  call alloc(y)
  y = bmad_com%real_8_map_init

  if (present(set_taylor)) then
    if (set_taylor) y = x   ! converts y to taylor (kind = 2)
  endif

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine universal_to_bmad_taylor (u_taylor, bmad_taylor, switch_z)
!
! Subroutine to convert from a universal_taylor map in Etienne's PTC 
! to a taylor map in BMAD.
! The conversion can also convert from the PTC coordinate convention:
!         (x, P_x, y, P_y, P_z, c*t = -z)
! to the BMAD coordinate convention:
!         (x, P_x, y, P_y, z, P_z)
!
! Modules needed:
!   use bmad
!
! Input:
!   u_taylor(6) -- Universal_taylor: Universal_taylor map.
!   switch_z    -- Logical, optional: If True then switch coordinate 
!                    conventions. Default is True.
!
! Output:
!   bmad_taylor(6)   -- Taylor_struct:
!-

Subroutine universal_to_bmad_taylor (u_taylor, bmad_taylor, switch_z)

  implicit none

  type (universal_taylor), intent(in) :: u_taylor(:)
  type (taylor_struct) :: bmad_taylor(:)

  integer i, j, k, ii, n

  logical, optional :: switch_z
  logical switch

! Remember to suppress any terms that have a zero coef.  

  do i = 1, 6

    switch = .true.
    if (present(switch_z)) switch = switch_z

    ii = i
    if (switch) then
      if (i == 5) ii = 6
      if (i == 6) ii = 5
    endif

    if (associated(bmad_taylor(i)%term)) deallocate(bmad_taylor(i)%term)

    n = count(u_taylor(ii)%c(:) /= 0)
    allocate(bmad_taylor(i)%term(n))

    k = 0
    do j = 1, u_taylor(ii)%n
      if (u_taylor(ii)%c(j) == 0) cycle
      k = k + 1
      if (switch) then
        bmad_taylor(i)%term(k)%exp  = u_taylor(ii)%j(j, (/1,2,3,4,6,5/))
        bmad_taylor(i)%term(k)%coef = u_taylor(ii)%c(j)
        bmad_taylor(i)%term(k)%coef = bmad_taylor(i)%term(k)%coef * &
                                      (-1)**bmad_taylor(i)%term(k)%exp(5)
        if (i == 5) bmad_taylor(i)%term(k)%coef = -bmad_taylor(i)%term(k)%coef
      else
        bmad_taylor(i)%term(k)%exp  = u_taylor(ii)%j(j,:)
        bmad_taylor(i)%term(k)%coef = u_taylor(ii)%c(j)
      endif
    enddo

  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine concat_real_8 (y1, y2, y3)
!
! Subroutine to concatinate two real_8 taylor series.
! This subroutine assumes that y1, y2, and y3 have been allocated.
!
! Modules needed:
!   use bmad
!
! Input:
!   y1(6) -- real_8: Input.
!   y2(6) -- real_8: Input.
!
! Output
!   y3(6) -- real_8: Concatinated output.
!-

subroutine concat_real_8 (y1, y2, y3)

  use s_fitting

  implicit none

  type (real_8), intent(in) :: y1(:), y2(:)
  type (real_8), intent(out) :: y3(:)
  type (damap) da1, da2, da3

! set the taylor order in PTC if not already done so

  if (bmad_com%taylor_order_ptc /= bmad_com%taylor_order) &
                         call set_ptc (taylor_order = bmad_com%taylor_order)

! Allocate temp vars

  call alloc(da1)
  call alloc(da2)
  call alloc(da3)

! concat

  da1 = y1
  da2 = y2

  da3 = da1 .o. da2  ! concat with constant terms
  
  y3 = da3

! kill temp vars

  call kill (da1)
  call kill (da2)
  call kill (da3)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_to_genfield (bmad_taylor, gen_field, c0)
!
! Subroutine to construct a genfield (partially inverted map) from a taylor
! map.
! Note: The constant terms of the taylor map are removed in the process.
! Note: The genfield uses PTC coordinates.
!
! Moudules needed:
!   use bmad
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Input taylor map.
!
! Output:
!   gen_field      -- Genfield: Output partially inverted map.
!   c0(6)          -- Real(rp): The constant part of the bmad_taylor map
!-

subroutine taylor_to_genfield (bmad_taylor, gen_field, c0)

  use s_fitting

  implicit none

  type (taylor_struct), intent(in) :: bmad_taylor(6)
  type (genfield), intent(inout) :: gen_field
  type (taylor_struct) taylor_(6)
  type (damap) da_map
  type (real_8) y(6)

  real(rp), intent(out) :: c0(6)

! set the taylor order in PTC if not already done so

  if (bmad_com%taylor_order_ptc /= bmad_com%taylor_order) &
                         call set_ptc (taylor_order = bmad_com%taylor_order)

! Remove constant terms from the taylor map first. This is probably
! not needed but we do it to make sure everything is alright.
! Also remove terms that have higher order then bmad_com%taylor_order

  call remove_constant_taylor (bmad_taylor, taylor_, c0, .true.)

! allocate pointers

  call alloc (gen_field)
  call alloc (da_map)
  call alloc (y)

! calculate the gen_field

  y = taylor_
  da_map = y
  gen_field = da_map

! cleanup

  call kill (da_map)
  call kill (y)
  call kill_taylor (taylor_)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine remove_constant_taylor (taylor_in, taylor_out, c0, &
!                                                 remove_higher_order_terms)
!
! Subroutine to remove the constant part of a taylor series.
! Optionally terms that are higher order than bmad_com%taylor_order can
! be removed.
!
! Note: It is assumed that taylor_out has been deallocated before the call to
! this routine. Calling this routine with the first two actual arguments the
! same is prohibited.
!
! Moudules needed:
!   use bmad
!
! Input:
!  taylor_in(6) -- Taylor_struct: Input taylor map.
!  remove_higher_order_terms -- Logical: If True then terms that are higher
!                               order than bmad_com%taylor_order are removed.
!
! Output:
!   taylor_out(6)  -- Taylor_struct: Taylor with constant terms removed.
!   c0(6)          -- Real(rp): The constant part of the taylor map
!-

subroutine remove_constant_taylor (taylor_in, taylor_out, c0, &
                                                 remove_higher_order_terms)

  implicit none

  type (taylor_struct), intent(in) :: taylor_in(:)
  type (taylor_struct) taylor_out(:)

  real(rp), intent(out) :: c0(:)

  integer i, j, n, nn, ss

  logical, intent(in) :: remove_higher_order_terms

!

  c0 = 0

  do i = 1, 6

    n = size(taylor_in(i)%term)

    do j = 1, size(taylor_in(i)%term)
      if (all(taylor_in(i)%term(j)%exp == 0)) then
        n = n - 1
        c0(i) = taylor_in(i)%term(j)%coef
      endif
      if (remove_higher_order_terms) then
        if (sum(taylor_in(i)%term(j)%exp) > bmad_com%taylor_order) n = n - 1
      endif
    enddo

    allocate (taylor_out(i)%term(n))

    nn = 0
    do j = 1, size(taylor_in(i)%term)
      ss = sum(taylor_in(i)%term(j)%exp)
      if (ss == 0 .or. (remove_higher_order_terms .and. &
                                            ss > bmad_com%taylor_order)) cycle
      nn = nn + 1
      taylor_out(i)%term(nn) = taylor_in(i)%term(j)
    enddo

  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_inverse (taylor_in, taylor_inv)
!
! Subroutine to invert a taylor map.
!
! Moudules needed:
!   use bmad
!
! Input:
!   taylor_in(6)  -- Taylor_struct: Input taylor map.
!
! Output:
!   taylor_inv(6) -- Taylor_struct: Inverted taylor map.
!-

subroutine taylor_inverse (taylor_in, taylor_inv)

  use s_fitting

  implicit none

  type (taylor_struct) :: taylor_in(:)
  type (taylor_struct) :: taylor_inv(:)
  type (taylor_struct) tlr(6)
  type (real_8) y(6), yc(6)
  type (damap) da

  real(rp) c0(6)
  real(8) c8(6)

! set the taylor order in PTC if not already done so

  if (bmad_com%taylor_order_ptc /= bmad_com%taylor_order) &
                         call set_ptc (taylor_order = bmad_com%taylor_order)

! The inverse operation of PTC ignores constant terms so we have to take
! them out and then put them back in.

  call remove_constant_taylor (taylor_in, tlr, c0, .true.)

  call alloc(da)
  call alloc(y)

! compute inverse

  y = tlr
  da = y
  da = da**(-1)
  y = da

! put constant terms back in.
! If the Map is written as:
!   R_out = T * R_in + C
! Then inverting:
!   R_in = T_inv * (R_out - C)

  if (any(c0 /= 0)) then
    call real_8_init(yc)
    call vec_bmad_to_ptc (c0, c8) ! Convert constant part to PTC coords
    c8 = -c8                      ! negate
    yc = c8                       ! Convert this to PTC taylor map
    call concat_real_8 (y, yc, y) 
    call kill (yc)
  endif

! transfer inverse to taylor_inv

  taylor_inv = y

! clean up

  call kill (da)
  call kill (y)
  call kill_taylor (tlr)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine concat_taylor (taylor1, taylor2, taylor3)
! 
! Subroutine to concatinate two taylor series:
!   taylor3(x) = taylor1(taylor2(x))  
!
! Modules needed:
!   use bmad
!
! Input:
!   taylor1(6) -- Taylor_struct: Taylor series.
!   taylor2(6) -- Taylor_struct: Taylor series.
!
! Output
!   taylor3(6) -- Taylor_struct: Concatinated series
!-

subroutine concat_taylor (taylor1, taylor2, taylor3)

  use s_fitting

  implicit none

  type (taylor_struct) :: taylor1(:), taylor2(:)
  type (taylor_struct) :: taylor3(:)
  type (real_8) y1(6), y2(6), y3(6)
  type (damap) da1, da2, da3

! set the taylor order in PTC if not already done so

  if (bmad_com%taylor_order_ptc /= bmad_com%taylor_order) &
                         call set_ptc (taylor_order = bmad_com%taylor_order)

! Allocate temp vars

  call real_8_init (y1)
  call real_8_init (y2)
  call real_8_init (y3)

! concat

  y1 = taylor1
  y2 = taylor2

  call concat_real_8 (y1, y2, y3)

  taylor3 = y3
  taylor3(:)%ref = taylor1(:)%ref

!

  call kill (y1)
  call kill (y2)
  call kill (y3)

end subroutine  
  
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine taylor_propagate1 (tlr, ele, param)
!
! Subroutine to track a real_8 taylor map through an element.
! The alternative routine if ele has a taylor series is concat_taylor.
!
! Modules needed:
!   use bmad
!
! Input:
!   tlr(6) -- Taylor_struct: Map to be tracked
!   ele    -- Ele_struct: Element to track through
!   param  -- Param_struct: 
!
! Output:
!   tlr(6)  -- Taylor_struct: Map through element
!-

subroutine taylor_propagate1 (tlr, ele, param)

  use s_tracking

  implicit none

  type (taylor_struct) tlr(:)
  type (real_8), save :: y(6)
  type (ele_struct) ele
  type (param_struct) param
  type (fibre), pointer, save :: a_fibre

  logical :: init_needed = .true.

! set the taylor order in PTC if not already done so

  if (bmad_com%taylor_order_ptc /= bmad_com%taylor_order) &
                         call set_ptc (taylor_order = bmad_com%taylor_order)

! init

  if (init_needed) then
    allocate (a_fibre)
    call real_8_init (y)
    init_needed = .false.
  endif

!

  y = tlr

  call alloc (a_fibre)
  call ele_to_fibre (ele, a_fibre, param)
  call ptc_track (a_fibre, y, default, +1)  ! "track" in PTC
  call kill (a_fibre)

  tlr = y

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ele_to_taylor (ele, param, orb0)
!
! Subroutine to make a taylor map for an element. 
! The order of the map is set by set_ptc
!
! Modules needed:
!   use bmad
!
! Input:
!   ele   -- Element_struct: 
!     %integration_ord  -- Order for the symplectic integrator: 2, 4, or 6.
!     %num_steps          -- Number of integrater steps.
!   orb0  -- Coord_struct, optional: Starting coords around which the Taylor series 
!              is evaluated.
!   param -- Param_struct: 
!     %beam_energy -- Needed for wigglers.
!
! Output:
!   ele -- Element_struct:
!     %taylor(6)  -- Taylor maps.
!-

subroutine ele_to_taylor (ele, param, orb0)

  use s_tracking

  implicit none
  
  type (ele_struct), intent(inout) :: ele
  type (param_struct), intent(in) :: param
  type (coord_struct), optional, intent(in) :: orb0

  type (fibre), pointer, save :: a_fibre
  type (real_8) y(6), y2(6)
  type (universal_taylor), save :: u_taylor(6)

  real(dp) x(6)
  
  integer i
  
  logical :: init_needed = .true.

! Init

  if (init_needed) then
    allocate (a_fibre)
    do i = 1, 6
      u_taylor(i) = 0  ! nullify
    enddo
    init_needed = .false.
  endif

  if (bmad_com%taylor_order_ptc /= bmad_com%taylor_order) then
    call set_ptc (taylor_order = bmad_com%taylor_order)
  endif

! Track with offset

  call alloc (a_fibre)
  call ele_to_fibre (ele, a_fibre, param)
 
  if (present(orb0)) then
    ele%taylor(:)%ref = orb0%vec
    call vec_bmad_to_ptc (orb0%vec, x)
  else
    ele%taylor(:)%ref = 0
    x = 0
  endif

  call real_8_init(y)
  y = x
  call ptc_track (a_fibre, y, default, +1) ! "track" in PTC

! take out the offset

  call real_8_init(y2)
  y2 = -x
  call concat_real_8 (y, y2, y)

! convert to bmad_taylor  

  do i = 1, 6
    u_taylor(i) = y(i)%t
  enddo
  
  call universal_to_bmad_taylor (u_taylor, ele%taylor)
  ele%taylor_order = bmad_com%taylor_order_ptc

  call kill(a_fibre)
  call kill(y)
  call kill(y2)

  if (associated (ele%gen_field)) call kill_gen_field (ele%gen_field)

! For wigglers there is a z_patch to take out the non-zero z offset that an
! on axis particle gets.

  if (ele%key == wiggler$) then
    do i = 1, size(ele%taylor(5)%term)
      if (all(ele%taylor(5)%term(i)%exp == 0)) then
        if (all(x == 0)) then    ! an on-axis particle defines the z_patch
          ele%value(z_patch$) = ele%taylor(5)%term(i)%coef
          ele%taylor(5)%term(i)%coef = 0
        else                     ! take out the z_patch
          ele%taylor(5)%term(i)%coef = ele%taylor(5)%term(i)%coef - &
                                                       ele%value(z_patch$)
        endif
        return
      endif
    enddo
  endif

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine type_real_8_taylors (y, switch_z)
!
! Subroutine to type out the taylor series from a real_8 array.
!
! Modules needed:
!   use bmad
!
! Input
!   y(6)     -- Real_8: 6 taylor series: (x, P_x, y, P_y, P_z, -z)
!   switch_z -- Logical, optional: If True then switch from PTC coordinate
!                       convention to BMAD's. Default is True.
!-

subroutine type_real_8_taylors (y, switch_z)

  implicit none

  type (real_8) y(:)
  type (taylor_struct) b_taylor(6)

  logical, optional :: switch_z

!

  call init_taylor (b_taylor)
  call real_8_to_taylor (y, b_taylor, switch_z)
  call type_taylors (b_taylor)
  call kill_taylor (b_taylor)

end subroutine

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
!   use bmad
!
! Input:
!   ut_in     -- Universal_taylor: Unsorted taylor series.
!
! Output:
!   ut_sorted -- Universal_taylor: Sorted taylor series.
!-

subroutine sort_universal_terms (ut_in, ut_sorted)

  use nr

  implicit none

  type (universal_taylor), intent(in)  :: ut_in
  type (universal_taylor) :: ut_sorted

  integer, allocatable :: ix_(:), ord_(:)
  integer i, j, n, nv, expn(6)

! init

  n = ut_in%n
  nv = ut_in%nv

  if (nv /= 6) then
    print *, 'ERROR IN SORT_UNIVERSAL_TERMS: I AM NOT SET UP FOR NV /= 6'
    call err_exit
  endif

  if (associated(ut_sorted%n)) &
              deallocate(ut_sorted%n, ut_sorted%nv, ut_sorted%c, ut_sorted%j)
  allocate(ut_sorted%n, ut_sorted%nv, ut_sorted%c(n), ut_sorted%j(n,nv), &
                                                              ix_(n), ord_(n))

  ut_sorted%n = n
  ut_sorted%nv = nv

!

  do i = 1, n
    expn = ut_in%j(i,:)
    ord_(i) = sum(expn)*10**6 + expn(6)*10**5 + expn(5)*10**4 + &
                expn(4)*10**3 + expn(3)*10**2 + expn(2)*10**1 + expn(1)
  enddo

  call indexx (ord_, ix_)

  do i = 1, n
    ut_sorted%c(i)= ut_in%c(ix_(i))
    ut_sorted%j(i,:)= ut_in%j(ix_(i),:)
  enddo

  deallocate(ord_, ix_)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine type_map (y)
!
! Subroutine to type the transfer maps of a real_8 array.
!
! Modules needed:
!   use bmad
!
! Input:
!   y(:)  -- Real_8: 
!-

subroutine type_map (y)

  use s_fitting

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
      print '(i6, f18.14, 20i3)', sum(ut%j(j,:)), ut%c(j), &
                                          (ut%j(j,k), k = 1, ut%nv)
    enddo
    ut = -1
  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+                                
! Subroutine ele_to_fibre (ele, fiber, param, integ_order, steps)
!
! Subroutine to convert a BMAD element to a PTC fibre element.
! This subroutine allocates fresh storage for the fibre so after calling
! this routine you need to deallocate at some point with:
!       call kill (fiber)
!
! Note: You need to call set_ptc before using this routine.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: BMAD element.
!   param       -- param_struct: 
!     %beam_energy     -- Beam energy (for wigglers).
!   integ_order -- Integer, optional: Order for the 
!                    sympletic integrator. Possibilities are: 2, 4, or 6
!                    Overrides ele%integration_ord.
!                    default = 2 (if not set with set_ptc).
!   steps       -- Integer, optional: Number of integration steps.
!                    Overrides ele%num_steps.
!                    If ele%num_steps = 0 and steps is not present
!                    then the default is used. The default = 10 if not
!                    set with set_ptc. 
!
! Output:
!   fiber -- Fibre: PTC fibre element.
!+

subroutine ele_to_fibre (ele, fiber, param, integ_order, steps)

  use mad_like

  implicit none
 
  type (ele_struct) ele
  type (fibre) fiber
  type (el_list) el
  type (param_struct) :: param

  real(8) mis_rot(6)

  real(rp) an0(0:n_pole_maxx), bn0(0:n_pole_maxx)
  real(rp) cos_t, sin_t, len, hk, vk, x_off, y_off

  integer i, n, metd_temp, nstd_temp, key, n_term
  integer, optional :: integ_order, steps

  character name*16
  logical kick_here

!

  el = 0  ! init: subroutine el_0

  el%name    = ele%name
  el%vorname = ele%type

  el%l    = ele%value(l$)
  el%ld   = ele%value(l$)
  el%lc   = ele%value(l$)

  el%tilt = ele%value(tilt$)

!

  key = ele%key
  len = ele%value(l$)
  if (.not. ele%is_on) key = drift$

  select case (key)

  case (drift$, rcollimator$, ecollimator$, monitor$, instrument$) 
    el%kind = kind1

  case (quadrupole$) 
    el%kind = matrix_kick_matrix  ! kind7
    el%k(2) = ele%value(k1$)

  case (sbend$) 
    el%kind = matrix_kick_matrix  ! kind7
    el%b0 = ele%value(g$)
    el%lc = ele%value(l_chord$)
    el%t1 = ele%value(e1$)
    el%t2 = ele%value(e2$)

  case (sextupole$)
    el%kind = drift_kick_drift  ! kind2
    el%k(3) = ele%value(k2$) / 2

  case (octupole$)
    el%kind = drift_kick_drift ! kind2
    el%k(4) = ele%value(k3$) / 6

  case (solenoid$)
    el%kind = kind17    ! kind5 will be quicker but not exact
    el%bsol = ele%value(ks$)

  case (sol_quad$)
    el%kind = kind17    ! kind5 will be quicker but not exact
    el%bsol = ele%value(ks$)
    el%k(2) = ele%value(k1$)

  case (marker$)
    el%kind = kind0

  case (kicker$, hkicker$, vkicker$)
    el%kind = kind2

  case (rfcavity$)
    el%kind = kind4
    el%volt = ele%value(volt$)
    if (param%total_length == 0) then
      el%freq0 = 1
    else
      el%freq0 = c_light / param%total_length
    endif
    el%lag = ele%value(phi0$)
    el%delta_e = 0

  case (elseparator$)
    el%kind = kind15
    hk = ele%value(hkick$) / len
    vk = ele%value(vkick$) / len
    if (hk == 0 .and. vk == 0) then
      el%tilt = 0
    else
      if (param%particle < 0) then
        hk = -hk
        vk = -vk
      endif
      el%tilt = -atan2 (hk, vk) + ele%value(tilt$)
    endif
    el%volt = 1e-6 * param%beam_energy * sqrt(hk**2 + vk**2)
    call multipole_ele_to_ab (ele, param%particle, an0, bn0, .false.) 
    if (any(an0 /= 0) .or. any(bn0 /= 0)) then
      print *, 'ERROR IN ELE_TO_FIBRE: ', &
                       'MULTIPOLES IN AN ELSEPARATOR NOT SUPPORTED IN A FIBRE.'
      call err_exit
    endif

  case (ab_multipole$, multipole$)
    el%kind = kind3

  case (beambeam$)
    print *, 'ERROR IN ELE_TO_FIBRE: BEAMBEAM ELEMENT NOT YET IMPLEMENTED!'
    call err_exit

  case (wiggler$)
    el%kind = kinduser2    
    if (ele%sub_key == periodic_type$) then
      print *, 'ERROR IN ELE_TO_FIBRE: OLD STYLE WIGGLER: ', ele%name
      print *, '       CANNOT BE USED WITH TAYLOR.'
      call err_exit
    endif

  case default
    print *, 'ERROR IN ELE_TO_FIBRE: UNKNOWN ELEMENT KEY: ', &
                                                 key_name(ele%key)
    print *, '      FOR ELEMENT: ', ele%name
    call err_exit

  end select

! multipole components
! bmad an and bn are integrated fields. PTC uses just the field.

  if (ele%key /= elseparator$) then
    kick_here = .false.
    if (ele%key == hkicker$ .or. ele%key == vkicker$) then
      hk = 0; vk = 0
      if (ele%key == hkicker$) hk = ele%value(kick$) / len
      if (ele%key == vkicker$) vk = ele%value(kick$) / len
      kick_here = .true.
    elseif (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0) then
      hk = ele%value(hkick$) / len
      vk = ele%value(vkick$) / len
      kick_here = .true.
    endif
    if (kick_here) then
      cos_t = cos(ele%value(tilt$))
      sin_t = sin(ele%value(tilt$))
      el%k(1)  = -hk * cos_t - vk * sin_t
      el%ks(1) = -hk * sin_t + vk * cos_t
    endif

    call multipole_ele_to_ab (ele, param%particle, an0, bn0, .false.)
    if (len /= 0) then
      an0 = an0 / len
      bn0 = bn0 / len
    endif

    n = min(n_pole_maxx+1, size(el%k))
    if (n-1 < n_pole_maxx) then
      if (any(an0(n:n_pole_maxx) /= 0) .or. any(bn0(n:n_pole_maxx) /= 0)) then
        print *, 'WARNING IN ELE_TO_FIBRE: MULTIPOLE NOT TRANSFERED TO FIBRE'
        print *, '        FOR: ', ele%name
      endif
    endif
 
    el%ks(1:n) = el%ks(1:n) + an0(0:n-1)
    el%k(1:n) = el%k(1:n) + bn0(0:n-1)

    if (key == sbend$) el%k(1) = el%k(1) + ele%value(g$) + ele%value(delta_g$)

    do n = nmax, 1, -1
      if (el%ks(n) /= 0 .or. el%k(n) /= 0) exit
    enddo
    el%nmul  = n
  endif

  if (el%kind == kind1 .and. el%nmul > 0) el%kind = matrix_kick_matrix  ! kind7

  if (el%kind == matrix_kick_matrix) el%nmul = max(2, el%nmul)

! METD and STEPS are global variables in PTC.

  if (present (integ_order)) then
    el%method = integ_order
  elseif (ele%integration_ord /= 0) then
    el%method = ele%integration_ord
  else
    el%method = METD
  endif

  if (present (steps)) then
    el%nst = steps
  elseif (ele%num_steps /= 0) then
    el%nst = ele%num_steps
  else
    el%nst = NSTD
  endif

  fiber = el

! wiggler

  if (key == wiggler$) then

    if (hyper_x$ /= hyperbolic_x$ .or. hyper_y$ /= hyperbolic_y$ .or. &
                                          hyper_xy$ /= hyperbolic_xy$) then
      print *, 'ERROR IN ELE_TO_FIBRE: WIGGLER FORM/TYPE MISMATCH!'
      print *, '     ', hyper_y$, hyper_xy$, hyper_x$
      print *, '     ', hyperbolic_y$, hyperbolic_xy$, hyperbolic_x$
      call err_exit
    endif

    n_term = size(ele%wig_term)
    call init_wig_pointers (fiber%mag%u2%w, n_term)   

    fiber%mag%u2%w%a(1:n_term) = c_light * &
            ele%value(polarity$) * ele%wig_term%coef / param%beam_energy
    fiber%mag%u2%w%k(1,1:n_term)  = ele%wig_term%kx
    fiber%mag%u2%w%k(2,1:n_term)  = ele%wig_term%ky
    fiber%mag%u2%w%k(3,1:n_term)  = ele%wig_term%kz
    fiber%mag%u2%w%f(1:n_term)    = ele%wig_term%phi_z
    fiber%mag%u2%w%form(1:n_term) = ele%wig_term%type

    call copy (fiber%mag, fiber%magp)
  endif

!  misalignments.
! In PTC the reference point for the offsets is the beginning of the element.
! In BMAD the reference point is the center of the element..

  x_off = ele%value(x_offset$)-ele%value(l$)*ele%value(x_pitch$)
  y_off = ele%value(y_offset$)-ele%value(l$)*ele%value(y_pitch$)

  mis_rot = (/ x_off, y_off, 0.0_rp, &
              -ele%value(y_pitch$), -ele%value(x_pitch$),  0.0_rp /)

  fiber = mis_rot  ! call fibre_mis

end subroutine

end module
