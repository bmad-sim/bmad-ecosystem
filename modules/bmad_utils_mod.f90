!+
! Module bmad_utils_mod
!
! Module for subroutines that use bmad_struct structures but do not
! call other routines in bmad_interface.
!
! ALSO: THESE ROUTINES DO NOT HAVE ACCESS TO THE OVERLOADED
! ELE1 = ELE2 AND RING1 = RING2.
!-

#include "CESR_platform.inc"

module bmad_utils_mod

  use bmad_struct

  interface reallocate
    module procedure reallocate_control_
  end interface

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
  elseif (particle == proton$ .or. particle == antiproton$) then
    mc2 = m_proton
  else
    print *, 'ERROR IN ENERGY_TO_KINETIC: UNKNOWN PARTICLE TYPE:', particle
    call err_exit
  endif

  if (energy < mc2) then
    print *, 'ERROR IN ENERGY_TO_KINETIC: ENERGY IS LESS THAN REST MASS:', &
                                                                        energy
    call err_exit
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
! Subroutine wiggler_vec_potential (ele, here, vec_pot)
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
!   here    -- Coord_struct: Coordinates for calculating the vector pot.
!
! Output:
!   vec_pot(3) -- Real(rp): Normalized vector potential
!-

subroutine wiggler_vec_potential (ele, here, vec_pot)

  implicit none

  type (ele_struct), target, intent(in) :: ele
  type (coord_struct), intent(in) :: here
  real(rp), intent(out) :: vec_pot(3)

  type (wig_term_struct), pointer :: t

  real(rp) c_x, s_x, c_y, s_y, c_z, s_z
  real(rp) x, y, s, coef

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
  ring_out%n_ele_use =            ring_in%n_ele_use
  ring_out%n_ele_ring =           ring_in%n_ele_use
  ring_out%n_ele_max =            ring_in%n_ele_max
  ring_out%n_control_max =        ring_in%n_control_max
  ring_out%n_ic_max =             ring_in%n_ic_max
  ring_out%input_taylor_order =   ring_in%input_taylor_order

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)
!
! Subroutine to transfer a Taylor map from one element to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele_in       -- Ele_struct: Element with the Taylor map.
!   taylor_order -- Integer: Order to truncate the Taylor map at.
!
! Output:
!   ele_out      -- Ele_struct: Element receiving the Taylor map truncated to
!                     order taylor_order.
!-

subroutine transfer_ele_taylor (ele_in, ele_out, taylor_order)

  implicit none

  type (ele_struct) ele_in, ele_out
  integer, intent(in) :: taylor_order
  integer it, ix, k 

!

  do it = 1, 6
    ix = 0
    do k = 1, size(ele_in%taylor(it)%term)
     if (sum(ele_in%taylor(it)%term(k)%exp(:)) <= taylor_order) ix = ix + 1
    enddo
    if (.not. associated(ele_out%taylor(it)%term)) &
                            allocate (ele_out%taylor(it)%term(ix))
    if (size(ele_out%taylor(it)%term) /= ix) &
                            allocate (ele_out%taylor(it)%term(ix))
    ix = 0
    do k = 1, size(ele_in%taylor(it)%term)
     if (sum(ele_in%taylor(it)%term(k)%exp(:)) <= taylor_order) then
        ix = ix + 1
        ele_out%taylor(it)%term(ix) = ele_in%taylor(it)%term(k)
      endif      
    enddo
  enddo

  ele_out%taylor_order = taylor_order
  ele_out%taylor(:)%ref = ele_in%taylor(:)%ref

  if (ele_in%key == wiggler$) ele_out%value(z_patch$) = ele_in%value(z_patch$)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_ring (ring, n)
!
! Subroutine to initialize a BMAD ring.
! 
! Modules needed:
!   use bmad
!
! Input:
!   n    -- Integer: Upper bound ring%ele_(0:) array is initialized to.
!
! Output:
!   ring -- Ring_struct: Initialized ring.
!-

subroutine init_ring (ring, n)

  implicit none

  type (ring_struct)  ring
  integer n

!

  call deallocate_ring_pointers (ring)
  call allocate_ring_ele_(ring, n)
  call init_ele (ring%ele_init)

  allocate (ring%control_(1000))
  allocate (ring%ic_(1000))

  ring%title = ' '
  ring%name = ' '
  ring%lattice = ' '
  ring%input_file_name = ' '

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+ 
! Function equivalent_eles (ele1, ele2) result (equiv)
!
! Subroutine to see if to elements are equivalent in terms of attributes so
! that their Taylor Maps would be the same. 
! If the reference orbit about which the Taylor map is made is zero then
! two elements can be equivalent even if the names are different.
!
! Modules needed:
!   use bmad
!
! Input: 
!   ele1 -- Ele_struct: Element with a Taylor map
!   ele2 -- Ele_struct: Element that might receive the Taylor map from ele1.
!
! Output:
!   equiv -- logical: True if elements are equivalent.
!-

function equivalent_eles (ele1, ele2) result (equiv)

  implicit none

  type (ele_struct), intent(in) :: ele1, ele2

  integer it

  logical equiv
  logical vmask(n_attrib_maxx)

!

  equiv = .false.

  if (ele1%key /= ele2%key) return
  if (ele1%sub_key /= ele2%sub_key) return

  if (ele1%name /= ele2%name .and. any(ele1%taylor%ref /= 0)) return

  vmask = .true.
  if (ele1%key == wiggler$ .and. ele1%sub_key == map_type$) &
      vmask((/k1$, rho$, b_max$, z_patch$, p0c$, check_sum$/)) = .false.
  if (any(ele1%value /= ele2%value .and. vmask)) return

  if (ele1%num_steps /= ele2%num_steps) return
  if (ele1%integration_ord /= ele2%integration_ord) return

  if (associated(ele1%wig_term) .neqv. associated(ele2%wig_term)) return
  if (associated(ele1%wig_term)) then
    if (size(ele1%wig_term) /= size(ele2%wig_term)) return
    do it = 1, size(ele1%wig_term)
      if (ele1%wig_term(it)%coef  /= ele2%wig_term(it)%coef)  cycle
      if (ele1%wig_term(it)%kx    /= ele2%wig_term(it)%kx)    cycle
      if (ele1%wig_term(it)%ky    /= ele2%wig_term(it)%ky)    cycle
      if (ele1%wig_term(it)%kz    /= ele2%wig_term(it)%kz)    cycle
      if (ele1%wig_term(it)%phi_z /= ele2%wig_term(it)%phi_z) cycle
    enddo
  endif

  equiv = .true.


end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine clear_ring_1turn_mats (ring)
!
! Subroutine to clear the 1-turn matrices in the ring structure:
!   ring%param%t1_no_RF
!   ring%param%t1_with_RF
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

  ring%param%t1_no_RF = 0
  ring%param%t1_with_RF = 0

end subroutine


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_ele (ele1, ele2)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end ele2's pointers point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1 -- Ele_struct:
!
! Output:
!   ele2 -- Ele_struct:
!-

subroutine transfer_ele (ele1, ele2)

  type (ele_struct), intent(in) :: ele1
  type (ele_struct), intent(out) :: ele2

  ele2 = ele1

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_eles (ele1, ele2)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end ele2's pointers point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1(:) -- Ele_struct:
!
! Output:
!   ele2(:) -- Ele_struct:
!-

subroutine transfer_eles (ele1, ele2)

  type (ele_struct), intent(in) :: ele1(:)
  type (ele_struct), intent(out) :: ele2(:)

  ele2 = ele1

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_ring (ring1, ring2)
!
! Subroutine to set ring2 = ring1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end ring2's pointers point to the same memory as ring1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   ring1 -- ring_struct:
!
! Output:
!   ring2 -- ring_struct:
!-

subroutine transfer_ring (ring1, ring2)

  type (ring_struct), intent(in) :: ring1
  type (ring_struct), intent(out) :: ring2

  ring2 = ring1

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord (coord_, n_coord)
!
! Subroutine to reallocate an allocatable  coord_struct array to at least:
!     coord_(0:n_coord)
! Note: The old coordinates are not saved except for coord_(0).
!  If at input coord_ is not allocated then coord_(0)%vec is set to zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   coord_(:) -- Coord_struct, allocatable: Allocatable array.
!   n_coord   -- Integer: Minimum array upper bound wanted.
!
! Output:
!   coord_(:) -- coord_struct: Allocated array.
!-

subroutine reallocate_coord (coord_, n_coord)

  type (coord_struct), allocatable :: coord_(:)
  type (coord_struct) start

  integer, intent(in) :: n_coord
  integer i

!

  if (allocated (coord_)) then
    if (size(coord_) < n_coord + 1) then
      start = coord_(0)
      deallocate (coord_)
      allocate (coord_(0:n_coord))
      coord_(0) = start
      do i = 1, n_coord
        coord_(i)%vec = 0
      enddo
    endif
  else
    allocate (coord_(0:n_coord))
    coord_(0)%vec = 0
    do i = 1, n_coord
      coord_(i)%vec = 0
    enddo
  endif

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Fuunction reallocate_control_(control_, n) result (new_control_)
!
! Function to reallocate the ring%control_(:) array.
! This data in the array will be saved.
! 
! Modules needed:
!   use bmad
!
! Input:
!   control_(:) -- Control_struct, pointer: Control Array
!   n           -- Integer: Array size for control_(:)
!
! Output:
!   new_control_(:) -- Control_struct, pointer: Allocated array.
!-

function reallocate_control_(control_, n) result (new_control_)

  implicit none

  type (control_struct), pointer :: control_(:), new_control_(:)
  integer, intent(in) :: n
  integer nn

!

  allocate (new_control_(n))
  if (associated(control_)) then
    nn = min(n, size(control_))
    new_control_(1:nn) = control_(1:nn)
  endif

end function

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine deallocate_ele_pointers (ele, nullify_only)
!
! Subroutine to deallocate the pointers in an element.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele -- ele_struct: Element with pointers.
!   nullify_only -- Logical, optional: If present and True then
!               Just nullify. Do not deallocate.
!
! Output:
!   ele -- Ele_struct: Element with deallocated pointers.
!-

subroutine deallocate_ele_pointers (ele, nullify_only)

  implicit none

  type (ele_struct) ele
  logical, optional, intent(in) :: nullify_only

! nullify only

  if (present (nullify_only)) then
    if (nullify_only) then
      nullify (ele%wig_term)
      nullify (ele%const)
      nullify (ele%r)
      nullify (ele%descrip)
      nullify (ele%a, ele%b)
      nullify (ele%wake)
      nullify (ele%taylor(1)%term, ele%taylor(2)%term, ele%taylor(3)%term, &
                ele%taylor(4)%term, ele%taylor(5)%term, ele%taylor(6)%term)
      nullify (ele%gen_field)
      return
    endif
  endif

! Normal deallocate

  if (associated (ele%wig_term)) deallocate (ele%wig_term)
  if (associated (ele%const))    deallocate (ele%const)
  if (associated (ele%r))        deallocate (ele%r)
  if (associated (ele%descrip))  deallocate (ele%descrip)
  if (associated (ele%a))        deallocate (ele%a, ele%b)

  if (associated (ele%wake)) then
    if (associated (ele%wake%sr)) deallocate (ele%wake%sr)
    if (associated (ele%wake%lr)) deallocate (ele%wake%lr)
    deallocate (ele%wake)
  endif

  if (associated (ele%taylor(1)%term)) deallocate &
           (ele%taylor(1)%term, ele%taylor(2)%term, ele%taylor(3)%term, &
           ele%taylor(4)%term, ele%taylor(5)%term, ele%taylor(6)%term)

  call kill_gen_field (ele%gen_field)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine kill_gen_field (gen_field)
!
! Subroutine to kill a gen_field.
!
! Modules needed:
!   use bmad
!
! Input:
!   gen_field -- Genfield, pointer: gen_field to kill.
!
! Output:
!   gen_field -- Genfield, pointer: Killed gen_field.
!-

subroutine kill_gen_field (gen_field)

  use tpsalie_analysis, only: kill 

  implicit none

  type (genfield), pointer :: gen_field

!

  if (associated(gen_field)) then
    call kill (gen_field)
    deallocate (gen_field)
  endif

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine init_ele (ele)
!
! Subroutine to initialize a Bmad element. Element is initialized to be free
! (not a lord or slave) and all %values set to zero.
!
! Modules needed:
!   use bmad
!
! Output:
!   ele -- Ele_struct: Initialized element.
!-

subroutine init_ele (ele)

  implicit none

  type (ele_struct)  ele

!

  ele%type = ' '
  ele%alias = ' '
  ele%name = '<Initialized>'
  ele%attribute_name = ' '
  
  ele%key = 0
  ele%sub_key = 0

  ele%value(:) = 0


  ele%control_type = free$
  ele%ix_value = 0
  ele%ic1_lord = 0
  ele%ic2_lord = -1
  ele%n_lord = 0
  ele%ix1_slave = 0
  ele%ix2_slave = -1
  ele%n_slave = 0
  ele%ix_pointer = 0
  ele%s = 0

  ele%floor%x = 0
  ele%floor%y = 0
  ele%floor%z = 0
  ele%floor%theta = 0
  ele%floor%phi   = 0
  ele%floor%psi   = 0

  ele%mat6_calc_method = bmad_standard$
  ele%tracking_method  = bmad_standard$
  ele%field_calc       = bmad_standard$
  ele%num_steps        = bmad_com%default_num_steps
  ele%integration_ord  = bmad_com%default_integ_order
  ele%ptc_kind = 0

  ele%is_on = .true.
  ele%multipoles_on = .true.
  ele%symplectify = .false.
  ele%exact_rad_int_calc = .false.
  ele%on_an_i_beam = .false.

  ele%field_master = .false.
  ele%aperture_at = exit_end$

  call deallocate_ele_pointers (ele)

! init Twiss

  ele%c_mat = 0
  ele%gamma_c = 1.0

  ele%x%beta  = 0
  ele%x%alpha = 0
  ele%x%gamma = 0
  ele%x%eta   = 0
  ele%x%etap  = 0
  ele%x%phi   = 0
  ele%x%sigma = 0

  ele%y%beta  = 0
  ele%y%alpha = 0
  ele%y%gamma = 0
  ele%y%eta   = 0
  ele%y%etap  = 0
  ele%y%phi   = 0
  ele%y%sigma = 0

! This is needed because of a compiler and/or totalview bug

  allocate (ele%r(1,1))
  ele%r = 0

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_ring_ele_ (ring, des_size)
!
! Subroutine to allocate or re-allocate the ele_ pointer in a ring.
! The upper bound of ring%ele_(0:n) will be des_size if it is present
! or the maximum of: (1000, 1.5*ring%ele_(:)).
!
! Modules needed:
!   use bmad
!
! Input:
!   ring     -- ring_struct: Ring with del_ pointer.
!   des_size -- integer, Optional: Optional desired upper bound for 
!                 ring%ele_(:).
!
! Output:
!   ring     -- ring_struct: Ring with re-allocated %ele_(:) pointer. 
!-

subroutine allocate_ring_ele_ (ring, des_size)

  implicit none

  type (ring_struct) ring
  integer, optional :: des_size

  type (ele_struct), pointer :: temp_ele(:)
  integer curr_n_ele, desired_size, i

! get new size

  desired_size = 1000
  if (associated (ring%ele_)) &
        desired_size = max ((3*size(ring%ele_))/2, desired_size)
  if (present(des_size))  desired_size = des_size

!  save ring%ele_ if present

  if (associated (ring%ele_)) then
    curr_n_ele = size (ring%ele_) - 1
    allocate (temp_ele(0:curr_n_ele))
    call transfer_eles (ring%ele_, temp_ele)
    deallocate (ring%ele_)
    allocate(ring%ele_(0:desired_size))
    call transfer_eles (temp_ele(0:curr_n_ele), ring%ele_(0:curr_n_ele))
    deallocate (temp_ele)
  else
    curr_n_ele = -1
    allocate(ring%ele_(0:desired_size))
  endif

! 

  do i = curr_n_ele+1, desired_size
    call init_ele (ring%ele_(i))
    ring%ele_(i)%ix_ele = i
  end do

  ring%n_ele_maxx = desired_size
  ring%beam_energy => ring%ele_(0)%value(beam_energy$)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine deallocate_ring_pointers (ring)
!
! Subroutine to deallocate the pointers in a ring.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring -- ring_struct: Ring with pointers.
!
! Output:
!   ring -- ring_struct: Ring with deallocated pointers.
!-

subroutine deallocate_ring_pointers (ring)

  implicit none

  type (ring_struct) ring
  integer i

!

  if (associated (ring%ele_)) then

    do i = lbound(ring%ele_, 1), ubound(ring%ele_, 1)
      call deallocate_ele_pointers (ring%ele_(i))
    enddo
    call deallocate_ele_pointers (ring%ele_init)

    deallocate (ring%ele_)
    deallocate (ring%control_)
    deallocate (ring%ic_)

  endif

  ring%n_ele_use  = -1
  ring%n_ele_max  = -1
  ring%n_ele_maxx = -1

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine transfer_mat_from_twiss (ele1, ele2, m)
!
! Subroutine to make a 6 x 6 transfer matrix from the twiss parameters
! at two points.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele1 -- Ele_struct: Element with twiss parameters for the starting point.
!     %x, %y -- a-mode and b-mode Twiss paramters
!       %beta   -- Beta parameter.
!       %alpha  -- Alpha parameter.
!       %phi    -- Phase at initial point.
!       %eta    -- Dispersion at initial point.
!       %etap   -- Dispersion derivative at initial point.
!     %c_mat(2,2) -- Coupling matrix
!   ele2 -- Ele_struct: Element with twiss parameters for the ending point.
!
! Output:
!   m(6,6) -- Real(rp): Transfer matrix between the two points.
!-

subroutine transfer_mat_from_twiss (ele1, ele2, m)

  implicit none

  type (ele_struct) ele1, ele2

  real(rp) m(6,6), v_mat(4,4), v_inv_mat(4,4), det

!

  call mat_make_unit (m)
  call transfer_mat2_from_twiss (ele1%x, ele2%x, m(1:2,1:2))
  call transfer_mat2_from_twiss (ele1%y, ele2%y, m(3:4,3:4))

  call mat_det (ele1%c_mat, det)
  ele1%gamma_c = sqrt(1-det)
  call make_v_mats (ele1, v_mat, v_inv_mat)
  m(1:4,1:4) = matmul (m(1:4,1:4), v_inv_mat)

  call mat_det (ele2%c_mat, det)
  ele2%gamma_c = sqrt(1-det)
  call make_v_mats (ele2, v_mat, v_inv_mat)
  m(1:4,1:4) = matmul (v_mat, m(1:4,1:4))

  m(1:2,6) = (/ ele2%x%eta, ele2%x%etap /) - &
                     matmul (m(1:2,1:2), (/ ele1%x%eta, ele1%x%etap /)) 
  m(3:4,6) = (/ ele2%y%eta, ele2%y%etap /) - &
                     matmul (m(3:4,3:4), (/ ele1%y%eta, ele1%y%etap /)) 

! The m(5,x) terms follow from the symplectic condition.

  m(5,1) = -m(2,6)*m(1,1) + m(1,6)*m(2,1) - m(4,6)*m(3,1) + m(3,6)*m(4,1)
  m(5,2) = -m(2,6)*m(1,2) + m(1,6)*m(2,2) - m(4,6)*m(3,2) + m(3,6)*m(4,2)
  m(5,3) = -m(2,6)*m(1,3) + m(1,6)*m(2,3) - m(4,6)*m(3,3) + m(3,6)*m(4,3)
  m(5,4) = -m(2,6)*m(1,4) + m(1,6)*m(2,4) - m(4,6)*m(3,4) + m(3,6)*m(4,4)


end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine match_ele_to_mat6 (ele, mat6, vec0)
!
! Subroutine to make the 6 x 6 transfer matrix from the twiss parameters
! at the entrance and exit ends of the element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele -- Ele_struct: Match element.
!     %value(beta_x0$) -- Beta_x at the start
!
! Output:
!   mat6(6,6) -- Real(rp): Transfer matrix.
!   vec0(6)   -- Real(rp): Currently just set to zero.
!-

subroutine match_ele_to_mat6 (ele, mat6, vec0)

  implicit none

  type (ele_struct) ele, ele0, ele1

  real(rp) mat6(6,6), vec0(6), v(n_attrib_maxx)

!

  vec0 = 0
  v = ele%value

  ele0%x%beta   = v(beta_x0$)
  ele0%x%alpha  = v(alpha_x0$)
  ele0%x%eta    = v(eta_x0$)
  ele0%x%etap   = v(etap_x0$)
  ele0%x%phi    = 0

  ele0%y%beta   = v(beta_y0$)
  ele0%y%alpha  = v(alpha_y0$)
  ele0%y%eta    = v(eta_y0$)
  ele0%y%etap   = v(etap_y0$)
  ele0%y%phi    = 0

  ele1%x%beta   = v(beta_x1$)
  ele1%x%alpha  = v(alpha_x1$)
  ele1%x%eta    = v(eta_x1$)
  ele1%x%etap   = v(etap_x1$)
  ele1%x%phi    = v(dphi_x$)

  ele1%y%beta   = v(beta_y1$)
  ele1%y%alpha  = v(alpha_y1$)
  ele1%y%eta    = v(eta_y1$)
  ele1%y%etap   = v(etap_y1$)
  ele1%y%phi    = v(dphi_y$)

  ele0%c_mat = 0 
  ele1%c_mat = 0 

  call transfer_mat_from_twiss (ele0, ele1, mat6)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function sr_wake_array_ubound (ele) result (array_ub)
!
! Function to return the size of ele%wake%sr.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele -- Element_struct: Element to test.
!
! Output:
!   array_ub -- Integer: Upper bound. Set to -1 if array not associated.
!-

function sr_wake_array_ubound (ele) result (array_ub)

  type (ele_struct) ele
  integer array_ub

!

  array_ub = -1
  if (.not. associated(ele%wake)) return
  if (associated (ele%wake%sr)) array_ub = ubound(ele%wake%sr, 1)

end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function lr_wake_array_size (ele) result (array_size)
!
! Function to return the size of ele%wake%lr.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele -- Element_struct: Element to test.
!
! Output:
!   array_size -- Integer: Set to 0 if array not associated.
!-

function lr_wake_array_size (ele) result (array_size)

  type (ele_struct) ele
  integer array_size

!

  array_size = 0
  if (.not. associated(ele%wake)) return
  if (associated (ele%wake%lr)) array_size = size(ele%wake%lr)

end function

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_sr_wake (sr_wake, n_ub)
!
! Subroutine to initialize a sr_wake array.
! Note: Lower bound is zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   sr_wake(:) -- Sr_wake_struct: Old structure.
!   n_ub       -- Integer: Upper bound of sr_wake(:) to allocate.
!                   -1 => pointer will be disassociated.
!
! Output:
!   sr_wake(0:n_ub) -- Sr_wake_struct: Initalized structure.
!-

subroutine init_sr_wake (sr_wake, n_ub)

  implicit none

  type (sr_wake_struct), pointer :: sr_wake(:)
  integer n_ub

!

  if (n_ub == -1) then
    if (associated (sr_wake)) deallocate (sr_wake)
    return
  endif

  if (associated (sr_wake)) then
    if (ubound(sr_wake, 1) /= n_ub) then
      deallocate (sr_wake)
      allocate (sr_wake(0:n_ub))
    endif
  else
    allocate (sr_wake(0:n_ub))
  endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_lr_wake (lr_wake, n_term)
!
! Subroutine to initialize a lr_wake array.
! Note: Lower bound is zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   lr_wake(:) -- Lr_wake_struct: Old structure.
!   n_term     -- Integer: Number of terms to allocate. 
!                   0 => pointer will be disassociated.
!
! Output:
!   lr_wake(0:n_term-1) -- Lr_wake_struct: Initalized structure.
!-

subroutine init_lr_wake (lr_wake, n_term)

  implicit none

  type (lr_wake_struct), pointer :: lr_wake(:)
  integer n_term

!

  if (n_term == 0) then
    if (associated (lr_wake)) deallocate (lr_wake)
    return
  endif

  if (associated (lr_wake)) then
    if (size(lr_wake) /= n_term) then
      deallocate (lr_wake)
      allocate (lr_wake(n_term))
    endif
  else
    allocate (lr_wake(n_term))
  endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_wake (wake, n_sr, n_lr)
!
! Subroutine to initialize a wake struct.
!
! Modules needed:
!   use bmad
!
! Input:
!   n_sr    -- Integer: Number of terms: wake%sr(0:n_sr) 
!   n_lr    -- Integer: Number of terms: wake%nr(n_lr)
!
! Output:
!   wake    -- Wake_struct, pointer: Initialized structure
!-

subroutine init_wake (wake, n_sr, n_lr)

  implicit none

  type (wake_struct), pointer :: wake
  integer n_sr, n_lr

!

  if (n_sr == -1 .and. n_lr == 0) then
    if (associated (wake)) then
      if (associated(wake%sr)) deallocate (wake%sr)
      if (associated(wake%lr)) deallocate (wake%lr)
      deallocate (wake)
    endif
  else
    if (.not. associated (wake)) allocate (wake)
    call init_sr_wake (wake%sr, n_sr)
    call init_lr_wake (wake%lr, n_lr)
  endif

end subroutine

end module

