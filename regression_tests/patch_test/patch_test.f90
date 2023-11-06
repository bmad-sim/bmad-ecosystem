!+
! Program patch_test
!
! This program is part of the Bmad regression testing suite.
!-

program patch_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, slave, slave2
type (branch_struct), pointer :: branch
type (coord_struct) :: start_orb, start2_orb, end_orb, end2_orb, end_orb_bs

integer ip, ie
character(40) fmt

!

em_field_custom_ptr => em_field_custom_test

open (1, file = 'output.now')
fmt = '(3a, 3f20.15, 5x, 3f20.15)'

call bmad_parser ('patch_test.bmad', lat)

!

branch => lat%branch(1)
call init_coord (start_orb, lat%particle_start, branch%ele(0), upstream_end$)

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)

  call track1 (start_orb, ele, branch%param, end_orb)
  call reverse_orbit(end_orb, start2_orb)
  call track1 (start2_orb, ele, branch%param, end2_orb)
  call reverse_orbit(end2_orb, end2_orb)
  write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-BS"  ABS 1E-12     ', end_orb%vec
  write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-BS-dif"  ABS 1E-12 ', end2_orb%vec - start_orb%vec

  if (valid_tracking_method(ele, positron$, runge_kutta$)) then
    end_orb_bs = end_orb
    ele%tracking_method = runge_kutta$

    call track1 (start_orb, ele, branch%param, end_orb)
    call reverse_orbit(end_orb, start2_orb)
    call track1 (start2_orb, ele, branch%param, end2_orb)
    call reverse_orbit(end2_orb, end2_orb)
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RK"  ABS 1E-12     ', end_orb%vec
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RK-dif"  ABS 1E-12 ', end2_orb%vec - start_orb%vec

    write (1, '(3a, f12.6)') '"', trim(ele%name), '-L"  ABS 1E-12 ', ele%value(l$)
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RK-BS"  ABS 1E-12  ', end_orb%vec - end_orb_bs%vec
    write (1, *)

    ele%field_calc = custom$
    call track1 (start_orb, ele, branch%param, end_orb)
    call reverse_orbit(end_orb, start2_orb)
    call track1 (start2_orb, ele, branch%param, end2_orb)
    call reverse_orbit(end2_orb, end2_orb)
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RKC"  ABS 1E-12     ', end_orb%vec
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RKC-dif"  ABS 1E-12 ', end2_orb%vec - start_orb%vec
  endif

  write (1, *)
enddo

!

do ip = 1, 3
  ele => lat%ele(ip)
  if (ele%key == marker$) cycle
  call init_coord (start_orb, lat%particle_start, ele, upstream_end$)
  call track1 (start_orb, ele, lat%param, end_orb)
  write (1, '(3a, 6es14.6)') '"', trim(ele%name), '" ABS 0', end_orb%vec
  if (ele%key == patch$) then
    write (1, '(a, f20.14)') '"L" REL 1E-12 ', ele%value(l$)
  endif
enddo

ele => lat%ele(4)
write (1, '(a, 6es14.6)') '"Flexible" ABS 1E-18 ', ele%floor%r, ele%floor%theta, ele%floor%phi, ele%floor%psi

!

branch => lat%branch(2)
do ip = 1, 3
  ele => branch%ele(ip)
  write (1, '(a, i0, a, f20.15)') '"L', ip, '-ref" ABS 1E-15 ', ele%value(l$)
enddo

!

close (1)

!------------------------------------------------------
contains

subroutine reverse_orbit(orb_in, orb_out)

type (coord_struct) orb_in, orb_out

orb_out = orb_in
orb_out%vec(2) = -orb_out%vec(2)
orb_out%vec(4) = -orb_out%vec(4)
orb_out%vec(5) = -orb_out%vec(5)
orb_out%direction = -1
orb_out%species = antiparticle(orb_out%species)

end subroutine reverse_orbit

!------------------------------------------------------
! contains

!+
! Subroutine em_field_custom (ele, param, s_rel, orbit, local_ref_frame, field, calc_dfield, err_flag)
!
! Routine for handling custom (user supplied) EM fields.
! This routine is called when ele%field_calc = custom$ or when ele is a custom element (ele%key = custom$)
! In order to be used, this stub file must be modified appropriately. See the Bmad manual for more details.
!
! Note: Unlike all other elements, "s_rel" and "here" areguments for a patch element are with respect to 
! the exit reference frame of the element. See the Bmad manual for more details.
!
! Note: Fields should not have any unphysical discontinuities. 
! Discontinuities may cause Runge-Kutta integration to fail resulting in particles getting marked as "lost".
! The mode of failure here is that RK will try smaller and smaller steps to integrate through the 
! discontinuity until the step size gets lower than bmad_com%min_ds_adaptive_tracking. At this
! point the particle gets marked as lost.
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below.
!
! Input:
!   ele         -- Ele_struct: Custom element.
!   param       -- lat_param_struct: Lattice parameters.
!   s_rel       -- Real(rp): Longitudinal position relative to the start of the element.
!   orbit       -- Coord_struct: Coords with respect to the reference particle.
!   local_ref_frame 
!               -- Logical, If True then take the 
!                     input coordinates and output fields as being with 
!                     respect to the frame of referene of the element. 
!   calc_dfield -- Logical, optional: If present and True then the field 
!                     derivative matrix is wanted by the calling program.
!
! Output:
!   field    -- Em_field_struct: Structure hoding the field values.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

recursive subroutine em_field_custom_test (ele, param, s_rel, orbit, local_ref_frame, field, calc_dfield, err_flag, &
                                       calc_potential, use_overlap, grid_allow_s_out_of_bounds, rf_time, used_eles)

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct), intent(in) :: orbit
type (em_field_struct) :: field
type (ele_pointer_struct), allocatable, optional :: used_eles(:)

real(rp), intent(in) :: s_rel
real(rp) f
logical local_ref_frame
real(rp), optional :: rf_time
logical, optional :: err_flag, grid_allow_s_out_of_bounds
logical, optional :: calc_dfield, calc_potential, use_overlap
character(*), parameter :: r_name = 'em_field_custom'

!

if (s_rel < -1 .or. s_rel > 10) then
  if (present(err_flag)) err_flag = .true.
  call out_io (s_fatal$, r_name, 'OUT OF RANGE!')
  call err_exit
endif

f = -0.001 * (1 + s_rel)
field%e = 0
field%b = [f*orbit%vec(3), -f*orbit%vec(1), 100*f*orbit%vec(1)]

if (present(err_flag)) err_flag = .false.

end subroutine em_field_custom_test

end program
