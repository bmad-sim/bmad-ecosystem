!+
! Module bmad_utils_mod
!
! Module for subroutines that use bmad_struct structures but do not
! call other routines in bmad_interface.
!
! ALSO: THESE ROUTINES DO NOT HAVE ACCESS TO THE OVERLOADED
! ELE1 = ELE2 AND LAT1 = LAT2.
!-

#include "CESR_platform.inc"

module bmad_utils_mod

use bmad_struct
use make_mat6_mod

interface reallocate
  module procedure reallocate_control
end interface

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine zero_ele_offsets (ele)
!
! Subroutine to zero the offsets, pitches and tilt of an element.
!
! Modules needed:
!   use bmad
!
! Output:
!   ele -- Ele_struct: Element with no (mis)orientation.
!-

subroutine zero_ele_offsets (ele)

implicit none

type (ele_struct) ele

ele%value(tilt$) = 0
ele%value(x_pitch$) = 0
ele%value(y_pitch$) = 0
ele%value(x_offset$) = 0
ele%value(y_offset$) = 0
ele%value(s_offset$) = 0

ele%value(tilt_tot$) = 0
ele%value(x_pitch_tot$) = 0
ele%value(y_pitch_tot$) = 0
ele%value(x_offset_tot$) = 0
ele%value(y_offset_tot$) = 0
ele%value(s_offset_tot$) = 0

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat6_add_pitch (ele, mat6)
!
! Subroutine to modify a first order transfer matrix to include the affect
! of an element pitch. Note that this routine does not correct the 0th order
! part of the map. It is assumed that on input the transfer map
! does not include the affect of any pitches.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element with pitches
!     %value(x_pitch_tot$) -- Horizontal pitch
!     %value(y_pitch_tot$) -- Vertical pitch
!   mat6(6,6) -- Real(rp): 1st order part of the transfer map (Jacobian).
!
! Output:
!   mat6(6,6) -- Real(rp): 1st order xfer map with pitches.
!-

subroutine mat6_add_pitch (ele, mat6)

implicit none

type (ele_struct) ele
real(rp) mat6(:,:), x_pitch, y_pitch

!

if (ele%value(x_pitch_tot$) == 0 .and. ele%value(y_pitch_tot$) == 0) return

x_pitch = ele%value(x_pitch_tot$)
y_pitch = ele%value(y_pitch_tot$)

mat6(5,6) = mat6(5,6) - mat6(5,2) * x_pitch - mat6(5,4) * y_pitch

mat6(5,1) = mat6(5,1) - x_pitch * (mat6(1,1) - 1) 
mat6(5,2) = mat6(5,2) - x_pitch *  mat6(1,2)
mat6(5,3) = mat6(5,3) - x_pitch *  mat6(1,3)
mat6(5,4) = mat6(5,4) - x_pitch *  mat6(1,4)

mat6(5,1) = mat6(5,1) - y_pitch *  mat6(3,1)
mat6(5,2) = mat6(5,2) - y_pitch *  mat6(3,2)
mat6(5,3) = mat6(5,3) - y_pitch * (mat6(3,3) - 1)
mat6(5,4) = mat6(5,4) - y_pitch *  mat6(3,4)

mat6(1,6) = mat6(5,2) * mat6(1,1) - mat6(5,1) * mat6(1,2) + &
                    mat6(5,4) * mat6(1,3) - mat6(5,3) * mat6(1,4)
mat6(2,6) = mat6(5,2) * mat6(2,1) - mat6(5,1) * mat6(2,2) + &
                    mat6(5,4) * mat6(2,3) - mat6(5,3) * mat6(2,4)
mat6(3,6) = mat6(5,4) * mat6(3,3) - mat6(5,3) * mat6(3,4) + &
                    mat6(5,2) * mat6(3,1) - mat6(5,1) * mat6(3,2)
mat6(4,6) = mat6(5,4) * mat6(4,3) - mat6(5,3) * mat6(4,4) + &
                    mat6(5,2) * mat6(4,1) - mat6(5,1) * mat6(4,2)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine convert_total_energy_to (E_tot, particle, 
!                                         gamma, kinetic, beta, pc, brho)
!
! Routine to calculate the momentum, etc. from a particle's total energy.
!
! Modules needed:
!   use bmad
!
! Input:
!   E_tot    -- Real(rp): Total energy of the particle.
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   gamma   -- Real(rp), optional: Gamma factor.
!   kinetic -- Real(rp), optional: Kinetic energy
!   beta    -- Real(rp), optional: velocity / c_light
!   pc      -- Real(rp), optional: Particle momentum
!   brho    -- Real(rp), optional: Nominal B_field*rho_bend
!-

subroutine convert_total_energy_to (E_tot, particle, &
                                         gamma, kinetic, beta, pc, brho)

  implicit none

  real(rp), intent(in) :: E_tot
  real(rp), intent(out), optional :: kinetic, beta, pc, brho, gamma
  real(rp) pc_new, mc2

  integer, intent(in) :: particle
  character(20) :: r_name = 'convert_total_energy_to'

!

  if (particle == positron$ .or. particle == electron$) then
    mc2 = m_electron
  elseif (particle == proton$ .or. particle == antiproton$) then
    mc2 = m_proton
  else
    call out_io (s_abort$, r_name, &
                    'ERROR: UNKNOWN PARTICLE TYPE:\i4\ ', particle)
    call err_exit
  endif

  if (E_tot < mc2) then
    call out_io (s_abort$, r_name, &
            'ERROR: TOTAL ENERGY IS LESS THAN REST MASS:\f10.0\ ', E_tot)
    call err_exit
  endif

  pc_new = E_tot * sqrt(1.0 - (mc2/E_tot)**2)
  if (present(pc))     pc     = pc_new
  if (present(beta))    beta    = pc_new / E_tot  
  if (present(kinetic)) kinetic = E_tot - mc2
  if (present(brho))    brho    = pc_new / c_light
  if (present(gamma))   gamma   = E_tot / mc2

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho)
!
! Routine to calculate the energy, etc. from a particle's momentum.
!
! Modules needed:
!   use bmad
!
! Input:
!   pc       -- Real(rp): Particle momentum
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   E_tot   -- Real(rp), optional: Total energy of the particle.
!   gamma   -- Real(rp), optional: Gamma factor.
!   kinetic -- Real(rp), optional: Kinetic energy
!   beta    -- Real(rp), optional: velocity / c_light
!   brho    -- Real(rp), optional: Nominal B_field*rho_bend
!-

subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho)

  implicit none

  real(rp), intent(in) :: pc
  real(rp), intent(out), optional :: E_tot, kinetic, beta, brho, gamma
  real(rp) E_tot_new, mc2

  integer, intent(in) :: particle
  character(20) :: r_name = 'convert_pc_to'

!

  if (particle == positron$ .or. particle == electron$) then
    mc2 = m_electron
  elseif (particle == proton$ .or. particle == antiproton$) then
    mc2 = m_proton
  else
    call out_io (s_abort$, r_name, &
                    'ERROR: UNKNOWN PARTICLE TYPE:\i4\ ', particle)
    call err_exit
  endif

  E_tot_new = sqrt(pc**2 + mc2**2)
  if (present(E_tot))   E_tot   = E_tot_new
  if (present(beta))    beta    = pc / E_tot_new
  if (present(kinetic)) kinetic = E_tot_new - mc2
  if (present(brho))    brho    = pc / c_light
  if (present(gamma))   gamma   = E_tot_new / mc2

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
! Subroutine transfer_lat_parameters (lat_in, lat_out)
!
! Subroutine to transfer the lat parameters (such as lat%name, lat%param, etc.)
! from one lat to another. The only stuff that is not transfered are the
! arrays and lat%ele_init.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_out -- lat_struct: Output lat with parameters set.
!-

subroutine transfer_lat_parameters (lat_in, lat_out)

  implicit none

  type (lat_struct), intent(in) :: lat_in
  type (lat_struct) :: lat_out

!

  lat_out%name =                 lat_in%name
  lat_out%lattice =              lat_in%lattice
  lat_out%input_file_name =      lat_in%input_file_name
  lat_out%title =                lat_in%title
  lat_out%a =                    lat_in%a
  lat_out%b =                    lat_in%b
  lat_out%z =                    lat_in%z
  lat_out%param =                lat_in%param
  lat_out%version =              lat_in%version
  lat_out%n_ele_track =          lat_in%n_ele_track
  lat_out%n_ele_max =            lat_in%n_ele_max
  lat_out%n_control_max =        lat_in%n_control_max
  lat_out%n_ic_max =             lat_in%n_ic_max
  lat_out%input_taylor_order =   lat_in%input_taylor_order
  lat_out%beam_start =           lat_in%beam_start

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
! Subroutine init_lat (lat, n)
!
! Subroutine to initialize a BMAD lat.
! 
! Modules needed:
!   use bmad
!
! Input:
!   n    -- Integer: Upper bound lat%ele(0:) array is initialized to.
!
! Output:
!   lat -- lat_struct: Initialized lat.
!-

subroutine init_lat (lat, n)

  implicit none

  type (lat_struct)  lat
  integer n

!

  call deallocate_lat_pointers (lat)
  call allocate_lat_ele(lat, n)
  call init_ele (lat%ele_init)

  allocate (lat%control(1000))
  allocate (lat%ic(1000))

  lat%title = ' '
  lat%name = ' '
  lat%lattice = ' '
  lat%input_file_name = ' '
  lat%param%stable = .true.

  lat%beam_start%vec = 0

  call init_mode_info (lat%a)
  call init_mode_info (lat%b)
  call init_mode_info (lat%z)

  lat%n_ele_track = 0
  lat%n_ele_max = 0
  lat%n_control_max = 0
  lat%n_ic_max = 0

  lat%input_taylor_order = -1
  lat%version = -1


!----------------------------------------
contains

subroutine init_mode_info (t)
  type (mode_info_struct) t
  t%tune = 0
  t%emit = 0
  t%chrom = 0
end subroutine

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

  type (ele_struct) :: ele1, ele2

  integer it

  logical equiv
  logical vmask(n_attrib_maxx)

!

  equiv = .false.

  if (ele1%key /= ele2%key) return
  if (ele1%sub_key /= ele2%sub_key) return
  if (ele1%map_with_offsets .neqv. ele2%map_with_offsets) return
  if (ele1%integrator_order /= ele2%integrator_order) return
  if (ele1%name /= ele2%name .and. any(ele1%taylor%ref /= 0)) return

  vmask = .true.
  if (ele1%key == wiggler$ .and. ele1%sub_key == map_type$) then
    vmask( (/ k1$, rho$, b_max$, z_patch$, p0c$, check_sum$ /) ) = .false.
  endif
  if (.not. ele1%map_with_offsets) then
    vmask( (/ x_offset$, y_offset$, s_offset$, tilt$, x_pitch$, &
              y_pitch$, x_offset_tot$, y_offset_tot$, s_offset_tot$, &
              tilt_tot$, x_pitch_tot$, y_pitch_tot$/) ) = .false.
  endif

  if (any(ele1%value /= ele2%value .and. vmask)) return

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
! Subroutine clear_lat_1turn_mats (lat)
!
! Subroutine to clear the 1-turn matrices in the lat structure:
!   lat%param%t1_no_RF
!   lat%param%t1_with_RF
! This will force any routine dependent upon these to do a remake.
!
! Modules needed:
!   use bmad
!
! Output:
!   lat -- lat_struct: Lat with 1-turn matrices cleared.
!-

subroutine clear_lat_1turn_mats (lat)

  implicit none

  type (lat_struct) lat

  lat%param%t1_no_RF = 0
  lat%param%t1_with_RF = 0

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

  type (ele_struct) :: ele1
  type (ele_struct) :: ele2

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

  type (ele_struct) :: ele1(:)
  type (ele_struct) :: ele2(:)

  ele2 = ele1

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_lat (lat1, lat2)
!
! Subroutine to set lat2 = lat1. 
! This is a plain transfer of information not using the overloaded equal.
! Thus at the end lat2's pointers point to the same memory as lat1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Modules needed:
!   use bmad
!
! Input:
!   lat1 -- lat_struct:
!
! Output:
!   lat2 -- lat_struct:
!-

subroutine transfer_lat (lat1, lat2)

  type (lat_struct), intent(in) :: lat1
  type (lat_struct), intent(out) :: lat2

  lat2 = lat1

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine reallocate_coord (coord, n_coord)
!
! Subroutine to reallocate an allocatable  coord_struct array to at least:
!     coord(0:n_coord)
! Note: The old coordinates are not saved except for coord(0).
! If, at input, coord(:) is not allocated then coord(0)%vec is set to zero.
! In any case, coord(n)%vec for n > 0 is set to zero.
!
! Modules needed:
!   use bmad
!
! Input:
!   coord(:) -- Coord_struct, allocatable: Allocatable array.
!   n_coord   -- Integer: Minimum array upper bound wanted.
!
! Output:
!   coord(:) -- coord_struct: Allocated array.
!-

subroutine reallocate_coord (coord, n_coord)

  type (coord_struct), allocatable :: coord(:)
  type (coord_struct) start

  integer, intent(in) :: n_coord
  integer i

!

  if (allocated (coord)) then
    if (size(coord) < n_coord + 1) then
      start = coord(0)
      deallocate (coord)
      allocate (coord(0:n_coord))
      coord(0) = start
      do i = 1, n_coord
        coord(i)%vec = 0
      enddo
    endif
  else
    allocate (coord(0:n_coord))
    do i = 0, n_coord
      coord(i)%vec = 0
    enddo
  endif

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Fuunction reallocate_control(control, n) result (new_control)
!
! Function to reallocate the lat%control(:) array.
! This data in the array will be saved.
! 
! Modules needed:
!   use bmad
!
! Input:
!   control(:) -- Control_struct, pointer: Control Array
!   n           -- Integer: Array size for control(:)
!
! Output:
!   new_control(:) -- Control_struct, pointer: Allocated array.
!-

function reallocate_control(control, n) result (new_control)

  implicit none

  type (control_struct), pointer :: control(:), new_control(:)
  integer, intent(in) :: n
  integer nn

!

  allocate (new_control(n))
  if (associated(control)) then
    nn = min(n, size(control))
    new_control(1:nn) = control(1:nn)
    deallocate (control)
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
      nullify (ele%a_pole, ele%b_pole)
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
  if (associated (ele%a_pole))        deallocate (ele%a_pole, ele%b_pole)

  if (associated (ele%wake)) then
    if (associated (ele%wake%sr_table))       deallocate (ele%wake%sr_table)
    if (associated (ele%wake%sr_mode_long))  deallocate (ele%wake%sr_mode_long)
    if (associated (ele%wake%sr_mode_trans)) deallocate (ele%wake%sr_mode_trans)
    if (associated (ele%wake%lr))        deallocate (ele%wake%lr)
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
  ele%closed_orb = 0

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
  ele%integrator_order = bmad_com%default_integ_order
  ele%ptc_kind  = 0
  ele%num_steps = 0

  ele%is_on             = .true.
  ele%multipoles_on     = .true.
  ele%symplectify       = .false.
  ele%map_with_offsets  = .true.
  ele%on_an_i_beam      = .false.
  ele%csr_calc_on       = .true.

  ele%field_master = .false.
  ele%aperture_at = exit_end$
  ele%coupler_at  = exit_end$

  call deallocate_ele_pointers (ele)

! init Twiss

  ele%c_mat = 0
  ele%gamma_c = 1.0

  ele%x%eta  = 0
  ele%x%etap = 0

  ele%y%eta  = 0
  ele%y%etap = 0

  ele%a%beta     = 0
  ele%a%alpha    = 0
  ele%a%gamma    = 0
  ele%a%eta      = 0
  ele%a%etap     = 0
  ele%a%phi      = 0
  ele%a%sigma    = 0
  ele%a%emit     = 0

  ele%b%beta     = 0
  ele%b%alpha    = 0
  ele%b%gamma    = 0
  ele%b%eta      = 0
  ele%b%etap     = 0
  ele%b%phi      = 0
  ele%b%sigma    = 0
  ele%b%emit     = 0

  ele%z%beta     = 0
  ele%z%alpha    = 0
  ele%z%gamma    = 0
  ele%z%eta      = 0
  ele%z%etap     = 0
  ele%z%phi      = 0
  ele%z%sigma    = 0
  ele%z%emit     = 0

! This is needed because of a compiler and/or totalview bug

  allocate (ele%r(1,1))
  ele%r = 0.0
end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine allocate_lat_ele (lat, des_size)
!
! Subroutine to allocate or re-allocate the ele pointer in a lat.
! The upper bound of lat%ele(0:n) will be des_size if it is present
! or the maximum of: (1000, 1.3*lat%ele(:)).
!
! Modules needed:
!   use bmad
!
! Input:
!   lat     -- lat_struct: Lat with del pointer.
!   des_size -- integer, Optional: Optional desired upper bound for 
!                 lat%ele(:).
!
! Output:
!   lat     -- lat_struct: Lat with re-allocated %ele(:) pointer. 
!-

subroutine allocate_lat_ele (lat, des_size)

  implicit none

  type (lat_struct) lat
  integer, optional :: des_size

  type (ele_struct), pointer :: temp_ele(:)
  integer curr_n_ele, desired_size, i

! get new size

  desired_size = 1000
  if (associated (lat%ele)) &
        desired_size = max (int(1.3*size(lat%ele)), desired_size)
  if (present(des_size))  desired_size = des_size

!  save lat%ele if present

  if (associated (lat%ele)) then
    curr_n_ele = size (lat%ele) - 1
    allocate (temp_ele(0:curr_n_ele))
    call transfer_eles (lat%ele, temp_ele)
    deallocate (lat%ele)
    allocate(lat%ele(0:desired_size))
    call transfer_eles (temp_ele(0:curr_n_ele), lat%ele(0:curr_n_ele))
    deallocate (temp_ele)
  else
    curr_n_ele = -1
    allocate(lat%ele(0:desired_size))
  endif

! 

  do i = curr_n_ele+1, desired_size
    call init_ele (lat%ele(i))
    lat%ele(i)%ix_ele = i
  end do

  lat%E_TOT => lat%ele(0)%value(E_TOT$)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine deallocate_lat_pointers (lat)
!
! Subroutine to deallocate the pointers in a lat.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- lat_struct: Lat with pointers.
!
! Output:
!   lat -- lat_struct: Lat with deallocated pointers.
!-

subroutine deallocate_lat_pointers (lat)

  implicit none

  type (lat_struct) lat
  integer i

!

  if (associated (lat%ele)) then

    do i = lbound(lat%ele, 1), ubound(lat%ele, 1)
      call deallocate_ele_pointers (lat%ele(i))
    enddo
    call deallocate_ele_pointers (lat%ele_init)

    deallocate (lat%ele)
    deallocate (lat%control)
    deallocate (lat%ic)

  endif

  lat%n_ele_track  = -1
  lat%n_ele_max  = -1

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
!     %a, %b -- a-mode and b-mode Twiss paramters
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
  character(20) :: r_name = 'transfer_mat_from_twiss'

! Error check

  if (ele1%a%beta == 0 .or. ele1%b%beta == 0) then
    call out_io (s_abort$, r_name, 'ZERO BETA IN ELEMENT: ' // ele1%name)
    call err_exit
  endif

  if (ele2%a%beta == 0 .or. ele2%b%beta == 0) then
    call out_io (s_abort$, r_name, 'ZERO BETA IN ELEMENT: ' // ele2%name)
    call err_exit
  endif

! Transfer matrices without coupling or dispersion

  call mat_make_unit (m)
  call transfer_mat2_from_twiss (ele1%a, ele2%a, m(1:2,1:2))
  call transfer_mat2_from_twiss (ele1%b, ele2%b, m(3:4,3:4))

! Add in coupling

  if (any(ele1%c_mat /= 0)) then
    call mat_det (ele1%c_mat, det)
    ele1%gamma_c = sqrt(1-det)
    call make_v_mats (ele1, v_mat, v_inv_mat)
    m(1:4,1:4) = matmul (m(1:4,1:4), v_inv_mat)
  endif

  if (any(ele2%c_mat /= 0)) then
    call mat_det (ele2%c_mat, det)
    ele2%gamma_c = sqrt(1-det)
    call make_v_mats (ele2, v_mat, v_inv_mat)
    m(1:4,1:4) = matmul (v_mat, m(1:4,1:4))
  endif

! Add in dispersion.
! The m(5,x) terms follow from the symplectic condition.

  m(1:2,6) = (/ ele2%a%eta, ele2%a%etap /) - &
                     matmul (m(1:2,1:2), (/ ele1%a%eta, ele1%a%etap /)) 
  m(3:4,6) = (/ ele2%b%eta, ele2%b%etap /) - &
                     matmul (m(3:4,3:4), (/ ele1%b%eta, ele1%b%etap /)) 

  m(5,1) = -m(2,6)*m(1,1) + m(1,6)*m(2,1) - m(4,6)*m(3,1) + m(3,6)*m(4,1)
  m(5,2) = -m(2,6)*m(1,2) + m(1,6)*m(2,2) - m(4,6)*m(3,2) + m(3,6)*m(4,2)
  m(5,3) = -m(2,6)*m(1,3) + m(1,6)*m(2,3) - m(4,6)*m(3,3) + m(3,6)*m(4,3)
  m(5,4) = -m(2,6)*m(1,4) + m(1,6)*m(2,4) - m(4,6)*m(3,4) + m(3,6)*m(4,4)


end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine match_ele_to_mat6 (ele, vec0, mat6)
!
! Subroutine to make the 6 x 6 transfer matrix from the twiss parameters
! at the entrance and exit ends of the element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele -- Ele_struct: Match element.
!     %value(beta_a0$) -- Beta_a at the start
!
! Output:
!   vec0(6)   -- Real(rp): Currently just set to zero.
!   mat6(6,6) -- Real(rp): Transfer matrix.
!-

subroutine match_ele_to_mat6 (ele, vec0, mat6)

  implicit none

  type (ele_struct) ele, ele0, ele1

  real(rp) mat6(6,6), vec0(6), v(n_attrib_maxx)

!

  vec0 = 0
  v = ele%value

  ele0%a%beta   = v(beta_a0$)
  ele0%a%alpha  = v(alpha_a0$)
  ele0%a%eta    = v(eta_a0$)
  ele0%a%etap   = v(etap_a0$)
  ele0%a%phi    = 0

  ele0%b%beta   = v(beta_b0$)
  ele0%b%alpha  = v(alpha_b0$)
  ele0%b%eta    = v(eta_b0$)
  ele0%b%etap   = v(etap_b0$)
  ele0%b%phi    = 0

  ele1%a%beta   = v(beta_a1$)
  ele1%a%alpha  = v(alpha_a1$)
  ele1%a%eta    = v(eta_a1$)
  ele1%a%etap   = v(etap_a1$)
  ele1%a%phi    = v(dphi_a$)

  ele1%b%beta   = v(beta_b1$)
  ele1%b%alpha  = v(alpha_b1$)
  ele1%b%eta    = v(eta_b1$)
  ele1%b%etap   = v(etap_b1$)
  ele1%b%phi    = v(dphi_b$)

  ele0%c_mat = 0 
  ele1%c_mat = 0 

  ele0%name = ele%name
  ele1%name = ele%name

  call transfer_mat_from_twiss (ele0, ele1, mat6)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_wake (wake_in, wake_out)
!
! Subroutine to transfer the wake info from one struct to another.
!
! Modules needed:
!   use bmad
!
! Input:
!   wake_in -- Wake_struct, pointer: Input wake.
!
! Output:
!   wake_out -- Wake_struct, pointer: Output wake.
!-

subroutine transfer_wake (wake_in, wake_out)

  implicit none

  type (wake_struct), pointer :: wake_in, wake_out
  integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr

!

  if (associated (wake_in)) then
    n_sr_table       = size(wake_in%sr_table)
    n_sr_mode_long  = size(wake_in%sr_mode_long)
    n_sr_mode_trans = size(wake_in%sr_mode_trans)
    n_lr        = size(wake_in%lr)
    call init_wake (wake_out, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)
    wake_out%sr_file   = wake_in%sr_file
    wake_out%lr_file   = wake_in%lr_file
    wake_out%z_sr_mode_max  = wake_in%z_sr_mode_max
    wake_out%sr_table       = wake_in%sr_table
    wake_out%sr_mode_long  = wake_in%sr_mode_long
    wake_out%sr_mode_trans = wake_in%sr_mode_trans
    wake_out%lr        = wake_in%lr
  else
    if (associated(wake_out)) call init_wake (wake_out, 0, 0, 0, 0)
  endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_wake (wake, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)
!
! Subroutine to initialize a wake struct.
!
! Modules needed:
!   use bmad
!
! Input:
!   n_sr_table       -- Integer: Number of terms: wake%sr_table(0:n_sr-1).
!   n_sr_mode_long  -- Integer: Number of terms: wake%nr(n_sr_mode_long).
!   n_sr_mode_trans -- Integer: Number of terms: wake%nr(n_sr_mode_trans).
!   n_lr        -- Integer: Number of terms: wake%nr(n_lr)
!
! Output:
!   wake -- Wake_struct, pointer: Initialized structure. 
!               If all inputs are 0 then wake is deallocated.
!-

subroutine init_wake (wake, n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr)

  implicit none

  type (wake_struct), pointer :: wake
  integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, n_lr

! Deallocate wake if all inputs are zero.

  if (n_sr_table == 0 .and. n_sr_mode_long == 0 .and. n_sr_mode_trans == 0 .and. n_lr == 0) then
    if (associated(wake)) then
      deallocate (wake%sr_table)
      deallocate (wake%sr_mode_long)
      deallocate (wake%sr_mode_trans)
      deallocate (wake%lr)
      deallocate (wake)
    endif
    return
  endif

!

  if (associated (wake)) then
    if (size(wake%sr_table) /= n_sr_table) then
      deallocate (wake%sr_table)
      allocate (wake%sr_table(0:n_sr_table-1))
    endif
    if (size(wake%sr_mode_long) /= n_sr_mode_long) then
      deallocate (wake%sr_mode_long)
      allocate (wake%sr_mode_long(n_sr_mode_long))
    endif
    if (size(wake%sr_mode_trans) /= n_sr_mode_trans) then
      deallocate (wake%sr_mode_trans)
      allocate (wake%sr_mode_trans(n_sr_mode_trans))
    endif
    if (size(wake%lr) /= n_lr) then
      deallocate (wake%lr)
      allocate (wake%lr(n_lr))
    endif

  else
    allocate (wake)
    allocate (wake%sr_table(0:n_sr_table-1))
    allocate (wake%sr_mode_long(n_sr_mode_long))
    allocate (wake%sr_mode_trans(n_sr_mode_trans))
    allocate (wake%lr(n_lr))
  endif

end subroutine

end module
