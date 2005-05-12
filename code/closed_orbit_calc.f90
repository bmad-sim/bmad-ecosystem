!+                           
! Subroutine closed_orbit_calc (ring, closed_orb, i_dim, direction)
!
! Subroutine to calculate the closed orbit for a circular machine.
! Closed_orbit_calc uses the 1-turn transfer matrix to converge upon a  
! solution. 
!
! i_dim = 5 simulates the affect of the RF that makes the beam change 
! its energy until the change of path length in the closed orbit over 
! one turn is zero.
!
! Note: This routine uses the 1-turn matrix ring%param%t1_no_RF or 
! ring%param%t1_with_RF in the computations. If you have changed conditions 
! significantly enough you might want to force a remake of the 1-turn matrices
! by calling clear_ring_1turn_mats.
!
! Note: See also closed_orbit_from_tracking as an alternative method
! of finding the closed orbit.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring            -- Ring_struct: Ring to track through.
!   closed_orb(0:)  -- Coord_struct, allocatable: closed_orb(0) is the 
!                       initial guess. closed_orb(0)%vec(6) is used
!                       as the energy around which the closed orbit
!                       is calculated if i_dim = 4 and varialbe_energy = False.
!   i_dim           -- Integer: Dimensions to use
!                     = 4  Transverse closed orbit at constant energy 
!                          (dE/E = closed_orb(0)%vec(6))
!                     = 5 Transverse closed orbit at constant energy with the
!                          energy adjusted so that vec(5) is the same 
!                          at the beginning and at the end.
!                     = 6  Full closed orbit for 6x6 matrix.
!   direction       -- Integer, optional: Direction of tracking. 
!                       +1 --> forwad (default), -1 --> backward.
!                       The closed orbit will be dependent on direction only
!                       in the case that radiation damping is turned on.
!
!   bmad_status -- Bmad status common block
!     %exit_on_error -- If True then subroutine will terminate program
!                         if the orbit does not converge.
!     %type_out      -- If True then the subroutine will type out
!                         a warning message if the orbit does not converge.
!
! Output:
!   closed_orb(0:) -- Coord_struct, allocatable: Closed orbit. closed_orb(i)
!                      is the orbit at the exit end of the ith element.
!   bmad_status    -- Bmad status common block
!     %ok          -- Set False if orbit does not converge, True otherwise.
!-

#include "CESR_platform.inc"

subroutine closed_orbit_calc (ring, closed_orb, i_dim, direction)

  use bmad_struct
  use bmad_interface, except => closed_orbit_calc
  use bookkeeper_mod, only: set_on_off, save_state$, restore_state$, off$

  implicit none

  type (ring_struct), target ::  ring
  type (ele_struct), pointer :: ele
  type (coord_struct)  del_co, del_orb
  type (coord_struct), allocatable, target ::  closed_orb(:)
  type (coord_struct), pointer :: start, end

  real(rp) mat2(6,6), mat(6,6)
  real(rp) :: amp_co, amp_del, t1(6,6), amp_del_old

  integer, optional :: direction
  integer i, n, n_ele, i_dim, i_max, dir, i1_int
  logical fluct_saved, aperture_saved

!----------------------------------------------------------------------
! init
! Random fluctuations must be turned off to find the closed orbit.

  call reallocate_coord (closed_orb, ring%n_ele_max)  ! allocate if needed

  fluct_saved = bmad_com%radiation_fluctuations_on
  bmad_com%radiation_fluctuations_on = .false.  

  aperture_saved = ring%param%aperture_limit_on
  ring%param%aperture_limit_on = .false.

  dir = integer_option(+1, direction)
  n_ele = ring%n_ele_use

  bmad_status%ok = .true.

  if (dir == +1) then
    start => closed_orb(0)
    end   => closed_orb(n_ele)
  else if (dir == -1) then
    start => closed_orb(n_ele)
    end   => closed_orb(0)
  else
    print *, 'ERROR IN CLOSED_ORBIT_CALC: BAD DIRECTION ARGUMENT.'
    call err_exit
  endif

!----------------------------------------------------------------------
! Further init

  select case (i_dim)

! Constant energy case
! Turn off RF voltage if i_dim == 4 (for constant delta_E)

  case (4, 5)
    n = 4
    if (all(ring%param%t1_no_RF == 0)) &
                call transfer_matrix_calc (ring, .false., ring%param%t1_no_RF)
    t1 = ring%param%t1_no_RF
    start%vec(5) = 0
    call set_on_off (rfcavity$, ring, save_state$)
    call set_on_off (rfcavity$, ring, off$)

    call make_mat2

    if (i_dim == 5) then  ! crude I1 integral calculation
      i1_int = 0
      do i = 1, ring%n_ele_use
        ele => ring%ele_(i)
        if (ele%key == sbend$) i1_int = i1_int + ele%value(l$) * &
              ele%value(g$) * (ring%ele_(i-1)%x%eta_lab + ele%x%eta_lab) / 2
      enddo
    endif

! Variable energy case: i_dim = 6

  case (6)
    n = 6
    if (all(ring%param%t1_with_RF == 0)) &
                call transfer_matrix_calc (ring, .true., ring%param%t1_with_RF)
    t1 = ring%param%t1_with_RF

    if (t1(6,5) == 0) then
      print *, 'ERROR IN CLOSED_ORBIT_CALC: CANNOT DO FULL 6-DIMENSIONAL'
      print *, '      CALCULATION WITH NO RF VOLTAGE!'
      bmad_status%ok = .false. 
      call err_exit
    endif

    call make_mat2

! Error

  case default
    print *, 'ERROR IN CLOSED_ORBIT_CALC: BAD I_DIM ARGUMENT:', i_dim
    bmad_status%ok = .false. 
    call err_exit
  end select
          
! Orbit correction = (T-1)^-1 * (orbit_end - orbit_start)
!                  = mat2     * (orbit_end - orbit_start)


!--------------------------------------------------------------------------
! Because of nonlinearities we may need to iterate to find the solution

  amp_del_old = 1e20  ! something large
  i_max = 100  

  do i = 1, i_max

    if (dir == +1) then
      call track_all (ring, closed_orb)
    else
      call track_many (ring, closed_orb, n_ele, 0, -1)
    endif

    if (i == i_max .or. ring%param%lost) then
      if (bmad_status%type_out) then
        if (ring%param%lost) then
          print *, 'ERROR IN CLOSED_ORBIT_CALC: ORBIT DIVERGING TO INFINITY!'
        else
          print *, 'ERROR IN CLOSED_ORBIT_CALC: NONLINEAR ORBIT NOT CONVERGING!'
        endif
      endif
      if (bmad_status%exit_on_error) call err_exit
      bmad_status%ok = .false.
      exit
    endif

    del_orb%vec = end%vec - start%vec
    del_co%vec(1:n) = matmul(mat2(1:n,1:n), del_orb%vec(1:n)) 
    if (i_dim == 5) del_co%vec(5) = del_orb%vec(5) / i1_int

    amp_co = sum(abs(start%vec(1:i_dim)))
    amp_del = sum(abs(del_co%vec(1:i_dim)))                                  

    if (amp_del < amp_co * bmad_com%rel_tollerance + bmad_com%abs_tollerance) exit

    if (amp_del < amp_del_old) then
      start%vec(1:i_dim) = start%vec(1:i_dim) + del_co%vec(1:i_dim)
    else  ! not converging so remake mat2 matrix
      call ring_make_mat6 (ring, -1, closed_orb)
      call transfer_matrix_calc (ring, .true., t1)
      call make_mat2 
      amp_del_old = 1e20  ! something large
    endif

    amp_del_old = amp_del

  enddo

! return rf cavities to original state

  if (n == 4) call set_on_off (rfcavity$, ring, restore_state$)

  bmad_com%radiation_fluctuations_on = fluct_saved  ! restore state
  ring%param%aperture_limit_on = aperture_saved

!------------------------------------------------------------------------------
contains

subroutine make_mat2 

  if (dir == -1) call mat_inverse (t1(1:n,1:n), t1(1:n,1:n))
  call mat_make_unit (mat(1:n,1:n))
  mat(1:n,1:n) = mat(1:n,1:n) - t1(1:n,1:n)
  call mat_inverse(mat(1:n,1:n), mat2(1:n,1:n))

end subroutine

end subroutine
