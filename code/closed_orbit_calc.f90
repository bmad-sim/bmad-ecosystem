!+                           
! Subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag)
!
! Subroutine to calculate the closed orbit for a circular machine.
! Closed_orbit_calc uses the 1-turn transfer matrix to converge upon a  
! solution. 
!
! For i_dim = 4 this routine tracks through the lattice with the RF turned off.
! and the particle energy will be determined by closed_orb(0)%vec(6) for direction = 1
! and closed_orb(n0)%vec(6) with n0 = lat%n_ele_track for direction = -1.
! The z component (closed_orb(i)%vec(5)) will be set to zero at the beginning
! of tracking and will generally not be zero at the end.
!
! i_dim = 5 simulates the affect of the RF that makes the beam change 
! its energy until the change of path length in the closed orbit over 
! one turn is zero. Tracking is done with RF off. This can be useful in
! determining the affect of kicks on the beam energy (this can, of course,
! also be done with i_dim = 6).
!
! i_dim = 6 finds the closed orbit with energy variation. The RF needs to be
! turned on in this case. Additionally, to simulate cases where the RF frequency
! is not a multiple of the revolution harmonic (EG in a dispersion measurement), 
! lat%absolute_time_tracking needs to be set to True and the phi0_ref attributes 
! of the RF cavities should be adjusted using auto_scale_field_phase_and_amp.
!
! Note: This routine uses the 1-turn matrix lat%param%t1_no_RF or 
! lat%param%t1_with_RF in the computations. If you have changed conditions 
! significantly enough you might want to force a remake of the 1-turn matrices
! by calling clear_lat_1turn_mats.
!
! The closed orbit calculation stops when the following condition is satisfied:
!   amp_del < amp_co * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking
! Where:
!   amp_co = sum(abs(closed_orb[at beginning of lattice]))
!   amp_del = sum(abs(closed_orb[at beginning] - closed_orb[at end]))
!   closed_orb = vector: (x, px, y, py, z, pz)
! closed_orb[at end] is calculated by tracking through the lattice starting 
! from closed_orb[at start].
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
!                     = 4  Transverse closed orbit at constant energy (RF off).
!                          (dE/E = closed_orb(0)%vec(6))
!                     = 5 Transverse closed orbit at constant energy (RF off) with 
!                          the energy adjusted so that vec(5) is the same 
!                          at the beginning and at the end.
!                     = 6 True closed orbit.
!   direction      -- Integer, optional: Direction of tracking. 
!                       +1 --> forwad (default), -1 --> backward.
!                       The closed orbit will be dependent on direction only
!                       in the case that radiation damping is turned on.
!   ix_branch      -- Integer, optional: Lattice branch to find the closed orbit of. 
!                       Default is 0 (main branch).
!
!   global_com    -- Global_common_struct: Bmad status common block
!     %type_out      -- If True then the subroutine will type out
!                         a warning message if the orbit does not converge.
!   bmad_com       -- Bmad_common_struct: Bmad common block.
!     %rel_tol_tracking -- Relative error. See above. Default = 1e-8
!     %abs_tol_tracking -- Absolute error. See above. Default = 1e-10
!
! Output:
!   closed_orb(0:) -- Coord_struct, allocatable: Closed orbit. closed_orb(i)
!                      is the orbit at the exit end of the ith element.
!   err_flag       -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine closed_orbit_calc (lat, closed_orb, i_dim, direction, ix_branch, err_flag)

use bmad_interface, except_dummy => closed_orbit_calc
use bookkeeper_mod, only: set_on_off, save_state$, restore_state$, off$

implicit none

type (lat_struct), target ::  lat
type (ele_struct), pointer :: ele, ele_start
type (branch_struct), pointer :: branch
type (coord_struct)  del_co, del_orb
type (coord_struct), allocatable, target ::  closed_orb(:)
type (coord_struct), pointer :: start, end

real(rp) mat2(6,6), t1(6,6)
real(rp) :: amp_co, amp_del, amp_del_old, i1_int, rf_freq, dt

integer, optional :: direction, ix_branch
integer i, n, n_ele, i_dim, i_max, dir, nc, track_state

logical, optional, intent(out) :: err_flag
logical fluct_saved, aperture_saved, damp_saved, err

character(20) :: r_name = 'closed_orbit_calc'

!----------------------------------------------------------------------
! init
! Random fluctuations must be turned off to find the closed orbit.

if (present(err_flag)) err_flag = .true.
branch => lat%branch(integer_option(0, ix_branch))

call reallocate_coord (closed_orb, branch%n_ele_max)  ! allocate if needed

fluct_saved = bmad_com%radiation_fluctuations_on
bmad_com%radiation_fluctuations_on = .false.  

aperture_saved = branch%param%aperture_limit_on
branch%param%aperture_limit_on = .false.

dir = integer_option(+1, direction)
n_ele = branch%n_ele_track

if (dir == +1) then
  start => closed_orb(0)
  end   => closed_orb(n_ele)
  ele_start => branch%ele(0)
else if (dir == -1) then
  start => closed_orb(n_ele)
  end   => closed_orb(0)
  ele_start => branch%ele(n_ele)
else
  call out_io (s_error$, r_name, 'BAD DIRECTION ARGUMENT.')
  return
endif

!----------------------------------------------------------------------
! Further init

n = i_dim  ! dimension of transfer matrix
nc = i_dim ! number of dimensions to compare.

select case (i_dim)

! Constant energy case
! Turn off RF voltage if i_dim == 4 (for constant delta_E)

case (4, 5)

  ! Check if rf is on and if so issue a warning message

  if (rf_is_on(branch)) then
    call out_io (s_warn$, r_name, 'Inconsistant calculation: RF ON with i_dim = \i4\ ', i_dim)
  endif

  !

  damp_saved  = bmad_com%radiation_damping_on
  bmad_com%radiation_damping_on = .false.  ! Want constant energy

  if (all(branch%param%t1_no_RF == 0)) &
              call transfer_matrix_calc (lat, .false., lat%param%t1_no_RF, ix_branch = branch%ix_branch)
  t1 = branch%param%t1_no_RF
  start%vec(5) = 0

  ! Save %old_is_on state in %bmad_logic to preserve it in case a calling routine is using it.

  branch%ele%bmad_logic = branch%ele%old_is_on
  call set_on_off (rfcavity$, lat, save_state$, ix_branch = branch%ix_branch)
  call set_on_off (rfcavity$, lat, off$, ix_branch = branch%ix_branch)

  call make_mat2 (err)
  if (err) then
    call end_cleanup
    return
  endif

  if (i_dim == 5) then  ! crude I1 integral calculation
    n = 4   ! Still only compute the transfer matrix for the transverse
    nc = 6  ! compare all 6 coords.
    i1_int = 0
    do i = 1, branch%n_ele_track
      ele => branch%ele(i)
      if (ele%key == sbend$) then
        i1_int = i1_int + ele%value(l$) * &
            ele%value(g$) * (branch%ele(i-1)%x%eta + ele%x%eta) / 2
      endif
    enddo
  endif

! Variable energy case: i_dim = 6

case (6)
  if (all(branch%param%t1_with_RF == 0)) &
              call transfer_matrix_calc (lat, .true., lat%param%t1_with_RF, ix_branch = branch%ix_branch)
  t1 = branch%param%t1_with_RF

  if (t1(6,5) == 0) then
    call out_io (s_error$, r_name, 'CANNOT DO FULL 6-DIMENSIONAL', &
                                   'CALCULATION WITH NO RF VOLTAGE!')
    return
  endif

  call make_mat2 (err)
  if (err) then
    call end_cleanup
    return
  endif

  ! Assume that frequencies are comensurate otherwise a closed orbit does not exist.

  if (branch%lat%absolute_time_tracking) then
    rf_freq = 0
    do i = 1, branch%n_ele_track
      if (branch%ele(i)%key /= rfcavity$) cycle
      if (.not. branch%ele(i)%is_on) cycle
      rf_freq = max(rf_freq, branch%ele(i)%value(rf_frequency$))
    enddo
  endif

! Error

case default
  call out_io (s_error$, r_name, 'BAD I_DIM ARGUMENT: \i4\ ', i_dim)
  return
end select
        
! Orbit correction = (T-1)^-1 * (orbit_end - orbit_start)
!                  = mat2     * (orbit_end - orbit_start)


!--------------------------------------------------------------------------
! Because of nonlinearities we may need to iterate to find the solution

amp_del_old = 1e20  ! something large
i_max = 100  

do i = 1, i_max

  if (dir == +1) then
    call track_all (lat, closed_orb, branch%ix_branch, track_state)
  else
    call track_many (lat, closed_orb, n_ele, 0, -1, branch%ix_branch, track_state)
  endif

  if (i == i_max .or. track_state /= moving_forward$) then
    if (global_com%type_out) then
      if (track_state /= moving_forward$) then
        call out_io (s_error$, r_name, 'ORBIT DIVERGING TO INFINITY!')
      else
        call out_io (s_error$, r_name, 'NONLINEAR ORBIT NOT CONVERGING!')
      endif
    endif
    call end_cleanup
    return
  endif

  del_orb%vec = end%vec - start%vec

  if (i_dim == 6 .and. branch%lat%absolute_time_tracking) then
    dt = (end%t - start%t) - nint((end%t - start%t) * rf_freq) / rf_freq
    del_orb%vec(5) = -end%beta * c_light * dt
  endif

  del_co%vec(1:n) = matmul(mat2(1:n,1:n), del_orb%vec(1:n)) 

  if (i_dim == 5) then
    del_co%vec(5) = 0
    del_co%vec(6) = del_orb%vec(5) / i1_int      
  endif

  amp_co = sum(abs(start%vec(1:nc)))
  amp_del = sum(abs(del_co%vec(1:nc)))

  ! We want to do at least one iteration to prevent problems when 
  ! only a small change is made to the machine and the closed orbit recalculated.

  if (i > 1 .and. amp_del < amp_co * bmad_com%rel_tol_tracking + bmad_com%abs_tol_tracking) exit

  if (amp_del < amp_del_old .and. modulo(i, 10) /= 0) then
    start%vec(1:nc) = start%vec(1:nc) + del_co%vec(1:nc)
    call init_coord (start, start, ele_start, upstream_end$)
    amp_del_old = amp_del
  else  ! not converging so remake mat2 matrix
    call lat_make_mat6 (lat, -1, closed_orb, branch%ix_branch)
    call transfer_matrix_calc (lat, .true., t1, ix_branch = branch%ix_branch)
    call make_mat2 (err)
    if (err) then
      call end_cleanup
      return
    endif

    amp_del_old = 1e20  ! something large
  endif

enddo

! Cleanup

call end_cleanup
if (present(err_flag)) err_flag = .false.

!------------------------------------------------------------------------------
! return rf cavities to original state

contains

subroutine end_cleanup

if (n == 4 .or. n == 5) then
  call set_on_off (rfcavity$, lat, restore_state$, ix_branch = branch%ix_branch)
  branch%ele%old_is_on = branch%ele%bmad_logic
  bmad_com%radiation_damping_on = damp_saved   ! restore state
endif

bmad_com%radiation_fluctuations_on = fluct_saved  ! restore state
lat%param%aperture_limit_on = aperture_saved

end subroutine

!------------------------------------------------------------------------------
! contains

subroutine make_mat2 (err)

real(rp) mat(6,6)
logical ok1, ok2, err

!

err = .true.

ok1 = .true.
if (dir == -1)  call mat_inverse (t1(1:n,1:n), t1(1:n,1:n), ok1)
call mat_make_unit (mat(1:n,1:n))
mat(1:n,1:n) = mat(1:n,1:n) - t1(1:n,1:n)
call mat_inverse(mat(1:n,1:n), mat2(1:n,1:n), ok2)

if (.not. ok1 .or. .not. ok2) then 
  if (global_com%type_out) call out_io (s_error$, r_name, 'MATRIX INVERSION FAILED!')
  return
endif
  
err = .false.

end subroutine

end subroutine
