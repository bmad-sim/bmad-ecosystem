!+                           
! Subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, exit_on_error)
!
! Subroutine to calculate the closed orbit for a circular machine.
! Closed_orbit_calc uses the 1-turn transfer matrix to converge upon a  
! solution. 
!
! i_dim = 5 simulates the affect of the RF that makes the beam change 
! its energy until the change of path length in the closed orbit over 
! one turn is zero.
!
! Note: This routine uses the 1-turn matrix lat%param%t1_no_RF or 
! lat%param%t1_with_RF in the computations. If you have changed conditions 
! significantly enough you might want to force a remake of the 1-turn matrices
! by calling clear_lat_1turn_mats.
!
! Note: See also closed_orbit_from_tracking as an alternative method
! of finding the closed orbit.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat            -- lat_struct: Lat to track through.
!   closed_orb(0:) -- Coord_struct, allocatable: closed_orb(n0) 
!                      is the initial guess where n0 = 0 for direction = 1 and 
!                      n0 = lat%n_ele_track for direction = -1. Additionally, 
!                      if i_dim = 4, then closed_orb(n0)%vec(6) is used as the energy 
!                      around which the closed orbit is calculated.
!   i_dim          -- Integer: Dimensions to use:
!                     = 4  Transverse closed orbit at constant energy 
!                          (dE/E = closed_orb(n0)%vec(6))
!                     = 5 Transverse closed orbit at constant energy with the
!                          energy adjusted so that vec(5) is the same 
!                          at the beginning and at the end.
!                     = 6  Full closed orbit for 6x6 matrix.
!   direction      -- Integer, optional: Direction of tracking. 
!                       +1 --> forwad (default), -1 --> backward.
!                       The closed orbit will be dependent on direction only
!                       in the case that radiation damping is turned on.
!   exit_on_error  -- Logical, optional: If True then subroutine will terminate 
!                         programif the orbit does not converge. Default is
!                         determined by bmad_status%exit_on_error
!
!   bmad_status    -- Bmad status common block
!     %exit_on_error -- Default for exit_on_error argument.
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

subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, exit_on_error)

  use bmad_struct
  use bmad_interface, except => closed_orbit_calc
  use bookkeeper_mod, only: set_on_off, save_state$, restore_state$, off$

  implicit none

  type (lat_struct), target ::  lat
  type (ele_struct), pointer :: ele
  type (coord_struct)  del_co, del_orb
  type (coord_struct), allocatable, target ::  closed_orb(:)
  type (coord_struct), pointer :: start, end

  real(rp) mat2(6,6), t1(6,6)
  real(rp) :: amp_co, amp_del, amp_del_old, i1_int

  integer, optional :: direction
  integer i, n, n_ele, i_dim, i_max, dir, nc

  logical, optional :: exit_on_error
  logical fluct_saved, aperture_saved, damp_saved

  character(20) :: r_name = 'closed_orbit_calc'

!----------------------------------------------------------------------
! init
! Random fluctuations must be turned off to find the closed orbit.

  call reallocate_coord (closed_orb, lat%n_ele_max)  ! allocate if needed

  fluct_saved = bmad_com%radiation_fluctuations_on
  bmad_com%radiation_fluctuations_on = .false.  

  aperture_saved = lat%param%aperture_limit_on
  lat%param%aperture_limit_on = .false.

  dir = integer_option(+1, direction)
  n_ele = lat%n_ele_track

  bmad_status%ok = .true.

  if (dir == +1) then
    start => closed_orb(0)
    end   => closed_orb(n_ele)
  else if (dir == -1) then
    start => closed_orb(n_ele)
    end   => closed_orb(0)
  else
    call out_io (s_error$, r_name, 'BAD DIRECTION ARGUMENT.')
    call err_exit
  endif

!----------------------------------------------------------------------
! Further init

  n = i_dim  ! dimension of transfer matrix
  nc = i_dim ! number of dimensions to compare.

  select case (i_dim)

! Constant energy case
! Turn off RF voltage if i_dim == 4 (for constant delta_E)

  case (4, 5)

    damp_saved  = bmad_com%radiation_damping_on
    bmad_com%radiation_damping_on = .false.  ! Want constant energy

    if (all(lat%param%t1_no_RF == 0)) &
                call transfer_matrix_calc (lat, .false., lat%param%t1_no_RF)
    t1 = lat%param%t1_no_RF
    start%vec(5) = 0
    call set_on_off (rfcavity$, lat, save_state$)
    call set_on_off (rfcavity$, lat, off$)

    call make_mat2 

    if (i_dim == 5) then  ! crude I1 integral calculation
      n = 4   ! Still only compute the transfer matrix for the transverse
      nc = 6  ! compare all 6 coords.
      i1_int = 0
      do i = 1, lat%n_ele_track
        ele => lat%ele(i)
        if (ele%key == sbend$) then
          i1_int = i1_int + ele%value(l$) * &
              ele%value(g$) * (lat%ele(i-1)%x%eta + ele%x%eta) / 2
        endif
      enddo
    endif

! Variable energy case: i_dim = 6

  case (6)
    if (all(lat%param%t1_with_RF == 0)) &
                call transfer_matrix_calc (lat, .true., lat%param%t1_with_RF)
    t1 = lat%param%t1_with_RF

    if (t1(6,5) == 0) then
      call out_io (s_error$, r_name, 'CANNOT DO FULL 6-DIMENSIONAL', &
                                     'CALCULATION WITH NO RF VOLTAGE!')
      bmad_status%ok = .false. 
      call err_exit
    endif

    call make_mat2 

! Error

  case default
    call out_io (s_error$, r_name, 'BAD I_DIM ARGUMENT: \i4\ ', i_dim)
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
      call track_all (lat, closed_orb)
    else
      call track_many (lat, closed_orb, n_ele, 0, -1)
    endif

    if (i == i_max .or. lat%param%lost) then
      if (bmad_status%type_out) then
        if (lat%param%lost) then
          call out_io (s_error$, r_name, 'ORBIT DIVERGING TO INFINITY!')
        else
          call out_io (s_error$, r_name, 'NONLINEAR ORBIT NOT CONVERGING!')
        endif
      endif
      if (logic_option(bmad_status%exit_on_error, exit_on_error)) call err_exit
      bmad_status%ok = .false.
      exit
    endif

    del_orb%vec = end%vec - start%vec
    del_co%vec(1:n) = matmul(mat2(1:n,1:n), del_orb%vec(1:n)) 
    if (i_dim == 5) then
      del_co%vec(5) = 0
      del_co%vec(6) = del_orb%vec(5) / i1_int      
    endif

    amp_co = sum(abs(start%vec(1:nc)))
    amp_del = sum(abs(del_co%vec(1:nc)))                                  

    if (amp_del < amp_co * bmad_com%rel_tolerance + bmad_com%abs_tolerance) exit

    if (amp_del < amp_del_old) then
      start%vec(1:nc) = start%vec(1:nc) + del_co%vec(1:nc)
    else  ! not converging so remake mat2 matrix
      call lat_make_mat6 (lat, -1, closed_orb)
      call transfer_matrix_calc (lat, .true., t1)
      call make_mat2 
      amp_del_old = 1e20  ! something large
    endif

    amp_del_old = amp_del

  enddo

! return rf cavities to original state

  if (n == 4 .or. n == 5) then
    call set_on_off (rfcavity$, lat, restore_state$)
    bmad_com%radiation_damping_on = damp_saved   ! restore state
  endif

  bmad_com%radiation_fluctuations_on = fluct_saved  ! restore state
  lat%param%aperture_limit_on = aperture_saved

!------------------------------------------------------------------------------
contains

subroutine make_mat2 

  real(rp) mat(6,6)
  logical ok1, ok2

!

  ok1 = .true.
  if (dir == -1)  call mat_inverse (t1(1:n,1:n), t1(1:n,1:n), ok1)
  call mat_make_unit (mat(1:n,1:n))
  mat(1:n,1:n) = mat(1:n,1:n) - t1(1:n,1:n)
  call mat_inverse(mat(1:n,1:n), mat2(1:n,1:n), ok2)

  if (.not. ok1 .or. .not. ok2) then 
    if (bmad_status%type_out) &
                 call out_io (s_error$, r_name, 'MATRIX INVERSION FAILED!')
    if (logic_option(bmad_status%exit_on_error, exit_on_error)) call err_exit
    bmad_status%ok = .false.
    return
  endif


end subroutine

end subroutine
