#include "CESR_platform.inc"

module bmad_utils_mod

  use bmad_struct
  use bmad_interface

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine energy_to_kinetic (energy, particle, 
!                                           gamma, kinetic, beta, p0c, brho)
!
! Subroutine to calculate the kinetic energy, etc. from a particle's energy.
!
! Modules needed:
!   use bmad
!
! Input:
!   energy   -- Real(rp): Energy of the particle.
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   gamma   -- Real(rp), optional: Gamma factor.
!   kinetic -- Real(rp), optional: Kinetic energy
!   beta    -- Real(rp), optional: velocity / c_light
!   p0c     -- Real(rp), optional: Particle momentum
!   brho    -- Real(rp), optional: Nominal B_field*rho_bend
!-

subroutine energy_to_kinetic (energy, particle, &
                                            gamma, kinetic, beta, p0c, brho)

  implicit none

  real(rp), intent(in) :: energy
  real(rp), intent(out), optional :: kinetic, beta, p0c, brho, gamma
  real(rp) p0c_, mc2

  integer, intent(in) :: particle

!

  if (particle == positron$ .or. particle == electron$) then
    mc2 = m_electron
  else
    mc2 = m_proton
  endif

  p0c_ = sqrt(energy**2 - mc2**2)
  if (present(p0c))     p0c     = sqrt(energy**2 - mc2**2)
  if (present(beta))    beta    = p0c_ / energy  
  if (present(kinetic)) kinetic = energy - mc2
  if (present(brho))    brho    = p0c_ / c_light
  if (present(gamma))   gamma   = energy / mc2

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine wiggler_vec_potential (ele, energy, here, vec_pot)
!
! Subroutine to calculate the normalized vector potential at 
! a point for a wiggler. The normalized potental a_norm is defined by:
!      p_cononical = p_mv - a_norm
! The Gauge used here is the same one as used in PTC and has A_x = 0.
! 
! Modules needed:
!   use bmad
!
! Input:
!   ele     -- Ele_struct: wiggler element.
!   energy  -- Real(rdef): Particle energy.
!   here    -- Coord_struct: Coordinates for calculating the vector pot.
!
! Output:
!   vec_pot(3) -- Real(rdef): Normalized vector potential
!-

subroutine wiggler_vec_potential (ele, energy, here, vec_pot)

  implicit none

  type (ele_struct), target, intent(in) :: ele
  type (coord_struct), intent(in) :: here
  real(rdef), intent(in) :: energy
  real(rdef), intent(out) :: vec_pot(3)

  type (wig_term_struct), pointer :: t

  real(rdef) c_x, s_x, c_y, s_y, c_z, s_z
  real(rdef) x, y, s, coef

  integer i

!

  if (ele%key /= wiggler$) then
    print *, 'ERROR IN WIGGLER_VEC_POTENTIAL. ELEMENT NOT A WIGGLER: ', &
                                                                 ele%name
    call err_exit
  endif

!

  x = here%vec(1)
  y = here%vec(3)
  s = here%vec(5)

  vec_pot = 0

  do i = 1, size(ele%wig_term)
    t => ele%wig_term(i)

      if (t%type == hyper_y$) then
        c_x = cos(t%kx * x)
        s_x = sin(t%kx * x)
      elseif (t%type == hyper_x$ .or. t%type == hyper_xy$) then
        c_x = cosh(t%kx * x)
        s_x = sinh(t%kx * x)
      else
        print *, 'ERROR IN WIGGLER_VEC_POTENTIAL: UNKNOWN TERM TYPE!'
        call err_exit
      endif

      if (t%type == hyper_y$ .or. t%type == hyper_xy$) then
        c_y = cosh (t%ky * y)
        s_y = sinh (t%ky * y)
      else
        c_y = cos (t%ky * y)
        s_y = sin (t%ky * y)
      endif

      c_z = cos (t%kz * s + t%phi_z)
      s_z = sin (t%kz * s + t%phi_z)

      coef = ele%value(polarity$) * t%coef

      vec_pot(2) = vec_pot(2) - coef  * (t%kz / (t%kx * t%ky)) * s_x * s_y * s_z
      vec_pot(3) = vec_pot(3) - coef  * (1 / t%kx)             * s_x * c_y * c_z
    enddo


end subroutine


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine transfer_ring_parameters (ring_in, ring_out)
!
! Subroutine to transfer the ring parameters (such as ring%name, ring%param, etc.)
! from one ring to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring_in -- Ring_struct: Input ring.
!
! Output:
!   ring_out -- Ring_struct: Output ring with parameters set.
!-

subroutine transfer_ring_parameters (ring_in, ring_out)

  implicit none

  type (ring_struct), intent(in) :: ring_in
  type (ring_struct) :: ring_out

!

  ring_out%name =                 ring_in%name
  ring_out%lattice =              ring_in%lattice
  ring_out%input_file_name =      ring_in%input_file_name
  ring_out%title =                ring_in%title
  ring_out%x =                    ring_in%x
  ring_out%y =                    ring_in%y
  ring_out%z =                    ring_in%z
  ring_out%param =                ring_in%param
  ring_out%version =              ring_in%version
  ring_out%n_ele_ring =           ring_in%n_ele_ring
  ring_out%n_ele_symm =           ring_in%n_ele_symm
  ring_out%n_ele_use =            ring_in%n_ele_use
  ring_out%n_ele_max =            ring_in%n_ele_max
  ring_out%n_control_array =      ring_in%n_control_array
  ring_out%n_ic_array =           ring_in%n_ic_array
  ring_out%input_taylor_order =   ring_in%input_taylor_order

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_taylor (ring_in, ring_out, type_out)
!
! Subroutine to transfer the taylor maps from the elements of one ring to
! the elements of another. The elements are matched between the rings so 
! that the appropriate element in ring_out will get the correct Taylor map
! even if the order of the elements is different in the 2 rings.
!
! Note: The transfered Taylor map will be truncated to bmad_com%taylor_order.
! Note: If the taylor_order of an element in ring_in is less than 
!   bmad_com%taylor_order then it will not be used.  
!
! Modules needed:
!   use bmad
!
! Input:
!   ring_in   -- Ring_struct: Input ring with Taylor maps.
!   type_out  -- Logical: If True then print a message for each Taylor map
!                 transfered.
!
! Output:
!   ring_out  -- Ring_struct: Ring to receive the Taylor maps.
!-

subroutine transfer_taylor (ring_in, ring_out, type_out)

  implicit none

  type (ring_struct), target, intent(in) :: ring_in
  type (ring_struct), target, intent(inout) :: ring_out
  type (ele_struct), pointer :: ele_in, ele_out

  integer i, j, k, it, ix
  integer n_in, ix_in(n_ele_maxx)
 
  logical, intent(in)  :: type_out
  logical vmask(n_attrib_maxx)  

! check global parameters

  if (ring_in%param%beam_energy /= ring_out%param%beam_energy) then
    if (type_out) then
      print *, 'TRANSFER_TAYLOR: THE RING ENERGIES ARE DIFFERENT.'
      print *, '    TAYLOR MAPS NOT TRANSFERED.'
    endif
    return
  endif

! Find the taylor series in the first ring.

  do i = 1, ring_in%n_ele_max
    if (associated(ring_in%ele_(i)%taylor(1)%term)) then
      if (bmad_com%taylor_order > ring_in%ele_(i)%taylor_order) cycle
      n_in = n_in + 1
      ix_in(n_in) = i
    endif
  enddo

! go through ring_out and match elements

  do i = 1, ring_out%n_ele_max

    ele_out => ring_out%ele_(i)

    vmask = .true.
    if (ele_out%key == wiggler$) vmask((/k1$, rho$, b_max$/)) = .false.

    do j = 1, n_in
      ele_in => ring_in%ele_(ix_in(j))
      if (ele_in%key /= ele_out%key) cycle
      if (ele_in%name /= ele_out%name) cycle
      if (any(ele_in%value /= ele_out%value .and. vmask)) cycle
      if (ele_in%num_steps /= ele_out%num_steps) cycle
      if (ele_in%integration_order /= ele_out%integration_order) cycle
      if (associated(ele_in%wig_term) .and. associated(ele_out%wig_term)) then
        if (size(ele_in%wig_term) /= size(ele_out%wig_term)) cycle
        do it = 1, size(ele_in%wig_term)
          if (ele_in%wig_term(it)%coef /= ele_out%wig_term(it)%coef) cycle
          if (ele_in%wig_term(it)%kx /= ele_out%wig_term(it)%kx) cycle
          if (ele_in%wig_term(it)%ky /= ele_out%wig_term(it)%ky) cycle
          if (ele_in%wig_term(it)%kz /= ele_out%wig_term(it)%kz) cycle
          if (ele_in%wig_term(it)%phi_z /= ele_out%wig_term(it)%phi_z) cycle
        enddo
      elseif (associated(ele_in%wig_term) .xor. &
                                          associated(ele_out%wig_term)) then
        cycle
      endif
      exit
    enddo

    if (j == n_in + 1) cycle

! we have a match so transfer the Taylor map.

    if (type_out) print *, &
                    'TRANSFER_TAYLOR: Reusing Taylor for: ', ele_in%name

    do it = 1, 6
      ix = 0
      do k = 1, size(ele_in%taylor(it)%term)
       if (sum(ele_in%taylor(it)%term(k)%exp(:)) <= &
                                      bmad_com%taylor_order) ix = ix + 1
      enddo
      allocate (ele_out%taylor(it)%term(ix))
      ix = 0
      do k = 1, size(ele_in%taylor(it)%term)
       if (sum(ele_in%taylor(it)%term(k)%exp(:)) <= &
                                      bmad_com%taylor_order) then
          ix = ix + 1
          ele_out%taylor(it)%term(ix) = ele_in%taylor(it)%term(k)
        endif      
      enddo
    enddo

    ele_out%taylor_order = bmad_com%taylor_order
    ele_out%taylor(:)%ref = ele_in%taylor(:)%ref

  enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine clear_ring_1turn_mats (ring)
!
! Subroutine to clear the 1-turn matrices in the ring structure:
!   ring%param%t1_mat4
!   ring%param%t1_mat6
! This will force any routine dependent upon these to do a remake.
!
! Modules needed:
!   use bmad
!
! Output:
!   ring -- ring_struct: Ring with 1-turn matrices cleared.
!-

subroutine clear_ring_1turn_mats (ring)

  implicit none

  type (ring_struct) ring

  ring%param%t1_mat4 = 0
  ring%param%t1_mat6 = 0

end subroutine

end module
