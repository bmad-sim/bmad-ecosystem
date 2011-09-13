module correct_ring_mod

  use super_recipes_mod
  use bmad
  use nr
  !use sim_meas_mod
  use sim_bpm_mod

  implicit none

  private lmfunc, ring_to_y

  type, private :: parameter_ele_struct
     character*40 name, param_name
     integer loc, param
     real(rp) orig_val, wt
  end type parameter_ele_struct

  type meas_struct
     character*40 name, param_name
     integer loc, param
     real(rp) wt !, rotation, x_offset, y_offset
  end type meas_struct

  integer, parameter, private :: n_det_max = 3000

  type (parameter_ele_struct), private, save :: p_ele(n_det_max)
  type (meas_struct) :: det(n_det_max)
  real(rp), private, pointer, save :: rms_param_ptr
  integer, private, save :: n_det, n_p_ele


!==========================================================================
! Define some logicals and structures to be used here and elsewhere

  ! Measurement parameters
  integer, parameter :: orbit_x$ = 1, orbit_y$ = 2, eta_xx$ = 3, eta_yy$ = 4
  integer, parameter :: phi_aa$ = 5, phi_bb$ = 6, cbar12$ = 7
  character(40), parameter :: allowed_meas_names(7) = (/'orbit_x', 'orbit_y', 'eta_xx ', 'eta_yy ', 'phi_aa ', 'phi_bb ', 'cbar12 ' /)

  type detcor_grp
     character*40 :: mask = "" ! note: 'mask' for detectors is superceded by correct(:)%bpm_mask.
     character(40) :: param_name = ""
     integer param
     real(rp) wt
  end type detcor_grp

  integer, parameter :: n_detcor_groups = 4
  type correct_struct
     character(40) :: bpm_mask = "^DET\_[0-9]{2}[ewEW]$" ! default to using CESR naming convention
     type(detcor_grp) :: det(n_detcor_groups)
     type(detcor_grp) :: cor(n_detcor_groups)
  end type correct_struct

  type correct_ring_params_struct
     real(rp) :: det_abs_res, det_diff_res, det_rot_res, n_lm_iterations, &
          eta_delta_e_e, sigma_cutoff
     logical :: write_elements, skip_dup_ele
  end type correct_ring_params_struct

  type(correct_ring_params_struct), save :: correct_ring_params

  type(lat_struct), save :: cr_model_ring
  type(coord_struct), allocatable, save :: cr_model_co(:)

contains

!==========================================================================
! Routine to correct a misaligned ring

  subroutine correct_ring(ring, bpm, correct, rms_param)
    implicit none
    
    type(lat_struct), intent(inout) :: ring
    type(correct_struct) correct
    real(rp), target :: rms_param
    type(det_struct) :: bpm(:)
    real(rp) :: current = 750000

    type measurement_struct
       real(rp) :: orbit_x, orbit_y, eta_x, eta_y, phi_a, phi_b, cbar12
    end type measurement_struct
    type (measurement_struct) meas(n_det_max)

    integer loc, ix, jx, n_bpms, bpm_ix

    type(lat_struct) :: hybrid_ring
    type(coord_struct), allocatable :: co(:), eta1(:), eta2(:), eta(:)
    type(ele_struct), pointer :: ele, ma_ele
    integer i_ele, i_det_grp, i_cor_grp, i_det, i, i_ele2
    integer size_y, iy, iter, ip
    integer status
    real(rp) harvest(7), cbar(2,2), chisq, alamda, old_chisq, d_chisq
    real(rp), allocatable, dimension(:) :: x, y, sig, a
    real(rp), allocatable, dimension(:,:) :: alpha, covar
    logical, allocatable, dimension(:) :: maska
    character(40) attrib_name

    logical :: verbose = .true.


    if (.not. allocated(eta1)) allocate(eta1(0:ring%n_ele_max), eta2(0:ring%n_ele_max), eta(0:ring%n_ele_max))

    n_bpms = size(bpm)

    ! Locate detectors

    rms_param_ptr => rms_param
    n_det  = 0
    do i_det_grp = 1, n_detcor_groups

       if (len(trim(correct%det(i_det_grp)%param_name)) == 0) then ! no name provided
          !write(*,'(a,i0,a)') "No parameter name provided for detector group ", i_det_grp,". Cycling."
          cycle
       endif

       ele_loop: do i_det = 1, size(bpm)

          ele => ring%ele(bpm(i_det)%ix_lat)

          n_det = n_det + 1

          ! look up index of desired parameter
          call match_word(trim(correct%det(i_det_grp)%param_name), allowed_meas_names, ix)
          correct%det(i_det_grp)%param = ix
          if (correct%det(i_det_grp)%param .eq. 0) then
             write(*,*) "Parameter ", trim(correct%det(i_det_grp)%param_name), &
                  " not found for detectors. Cycling." 
             cycle ele_loop
          endif

          ! can't say det(i_det)%blah because n_det > size(bpm)
          det(n_det)%name = ele%name
          det(n_det)%loc = bpm(i_det)%ix_lat
          det(n_det)%param = correct%det(i_det_grp)%param
          det(n_det)%param_name = correct%det(i_det_grp)%param_name
          det(n_det)%wt  = correct%det(i_det_grp)%wt

       end do ele_loop

    end do

    ! Locate elements to be varied for a particular correction
    n_p_ele = 0
    do i_cor_grp = 1, n_detcor_groups
       
       if (len(trim(correct%cor(i_cor_grp)%param_name)) == 0) cycle
       
       ele_loop2: do i_ele = 1, ring%n_ele_max ! steerings are overlays, not tracked elements
          
          ele => ring%ele(i_ele)
          
          if (match_reg(trim(ele%name), trim(correct%cor(i_cor_grp)%mask))) then
             
             if (len(trim(correct%cor(i_cor_grp)%param_name)) == 0) then
                !write(*,'(a,i0,a)') "No parameter name provided for corrector group ", i_cor_grp,". Cycling."
                cycle
             endif
             
             correct%cor(i_cor_grp)%param = attribute_index(ele, trim(correct%cor(i_cor_grp)%param_name))
             
             if (correct%cor(i_cor_grp)%param .eq. 0) then
                !write(*,*) "Parameter ", trim(correct%cor(i_cor_grp)%param_name), &
                !     " not found for element ", trim(ele%name), ". Cycling." 
                cycle ele_loop2
             endif
             
             ! By checking if attribute is free, we get around having to do a second loop
             ! to make sure the element is unique. -JSh
             
             if (.not. attribute_free(i_ele, correct%cor(i_cor_grp)%param_name, ring, .false.)) then
                !write(*,*) "WARNING: Attribute ", trim(attrib_name), " not free for element ", trim(ele%name), &
                !     " - skipping"
                cycle
             endif             
             
             n_p_ele = n_p_ele + 1
             p_ele(n_p_ele)%name  = ele%name
             p_ele(n_p_ele)%loc   = i_ele
             p_ele(n_p_ele)%param = correct%cor(i_cor_grp)%param
             p_ele(n_p_ele)%param_name = correct%cor(i_cor_grp)%param_name
             p_ele(n_p_ele)%wt    = correct%cor(i_cor_grp)%wt
             
             if (correct%cor(i_cor_grp)%param == -1) then
                call multipole_init(ring%ele(i_ele))
                call multipole_init(cr_model_ring%ele(i_ele))
             end if
             if (correct%cor(i_cor_grp)%param > 0) then
                p_ele(n_p_ele)%orig_val = cr_model_ring%ele(i_ele)%value(p_ele(n_p_ele)%param)
             else
                p_ele(n_p_ele)%orig_val = cr_model_ring%ele(i_ele)%a_pole(1)
             end if
          endif
       end do ele_loop2
       
       ! Check for too many parameter elements
       if (n_p_ele == size(p_ele)) then
          write(*,*) "Too many parameter elements.  Increase size(p_ele)."
          call err_exit
       end if
    end do
    
!==========================================================================
! Write out list of detectors and parameter elements

    if (correct_ring_params%write_elements) then
       write(*,*) " "
       write(*,'(a,i0,a)') "Found ", n_det, " detectors."
       if (verbose .eqv. .true.) then
          do i_ele = 1, n_det
             write(*,'(a10,a,a10,g10.2)', advance='no') &
                  trim(det(i_ele)%name), '  ', &
                  det(i_ele)%param_name, det(i_ele)%wt
             if (mod(i_ele,3) == 0) then
                write(*,*)
             else
                write(*,'(a3)', advance='no') " | "
             end if
          end do
          write(*,*)
       endif ! if verbose
          
       write(*,'(a,i0,a)') "Found ", n_p_ele, " parameter elements."
       if (verbose .eqv. .true.) then
          do i_ele = 1, n_p_ele
             write(*,'(a10,a,a10,g10.2)', advance='no') trim(p_ele(i_ele)%name), '  ',&
                  p_ele(i_ele)%param_name, p_ele(i_ele)%wt
             if (mod(i_ele,3) == 0) then
                write(*,*)
             else
                write(*,'(a3)', advance='no') " | "
             end if
          end do
          write(*,*)
       else ! if NOT verbose
          write(*,*) "*************Verbose mode disabled!*****************"
          write(*,*) "To see full list of correctors, set verbose = .true."
       endif
       write(*,*) " " 
    end if
       
!==========================================================================
! Store simulated measurements with BPM errors

    call twiss_and_track(ring, co)

    ! dispersion (orbit 1):
    bpm(:)%vec(1) = co(bpm(:)%ix_lat)%vec(1) + (ring%ele(bpm(:)%ix_lat)%x%eta * correct_ring_params%eta_delta_e_e)/2.
    bpm(:)%vec(3) = co(bpm(:)%ix_lat)%vec(3) + (ring%ele(bpm(:)%ix_lat)%y%eta * correct_ring_params%eta_delta_e_e)/2.
    do jx = 1, n_bpms
       call rotate_meas(bpm(jx)%vec(1),bpm(jx)%vec(3),bpm(jx)%tilt)
       call apply_bpm_errors(bpm(jx), current)  ! note: also transfers vec(6)-->amp(4) button signals
       
       eta1(bpm(jx)%ix_lat)%vec = bpm(jx)%vec

    enddo


    ! dispersion (orbit 2):
    bpm(:)%vec(1) = co(bpm(:)%ix_lat)%vec(1) - (ring%ele(bpm(:)%ix_lat)%x%eta * correct_ring_params%eta_delta_e_e)/2.
    bpm(:)%vec(3) = co(bpm(:)%ix_lat)%vec(3) - (ring%ele(bpm(:)%ix_lat)%y%eta * correct_ring_params%eta_delta_e_e)/2.
    do jx = 1, n_bpms
       call rotate_meas(bpm(jx)%vec(1),bpm(jx)%vec(3),bpm(jx)%tilt)
       call apply_bpm_errors(bpm(jx), current)  ! note: also transfers vec(6)-->amp(4) button signals

       eta2(bpm(jx)%ix_lat)%vec = bpm(jx)%vec
       eta(jx)%vec = (eta1(jx)%vec - eta2(jx)%vec) / correct_ring_params%eta_delta_e_e

    enddo


    ! orbit:
    bpm(:)%vec(1) = co(bpm(:)%ix_lat)%vec(1)
    bpm(:)%vec(3) = co(bpm(:)%ix_lat)%vec(3)
    do jx = 1, n_bpms
       call rotate_meas(bpm(jx)%vec(1),bpm(jx)%vec(3),bpm(jx)%tilt)
       call apply_bpm_errors(bpm(jx), current)  ! note: also transfers vec(6)-->amp(4) button signals

       ! transfer back to co for now:
       co(bpm(jx)%ix_lat)%vec(:) = bpm(jx)%vec
    enddo


    ! load relevent values into meas structure:
    do i_det = 1, n_det

       loc = det(i_det)%loc
       ele => ring%ele(loc)

       ! there might be an easier way to align bpm(:) and det(:)...
       do ix = 1, size(bpm)
          if (bpm(ix)%ix_lat .eq. loc) then
             bpm_ix = ix
             exit
          endif
       enddo

       meas(i_det)%orbit_x = co(loc)%vec(1)
       meas(i_det)%orbit_y = co(loc)%vec(3)
       meas(i_det)%eta_x   = eta(loc)%vec(1)
       meas(i_det)%eta_y   = eta(loc)%vec(3)
       meas(i_det)%phi_a = ele%a%phi
       meas(i_det)%phi_b = ele%b%phi
       call c_to_cbar(ele, cbar)
       call rotate_cbar(cbar, bpm(bpm_ix)%tilt)
       meas(i_det)%cbar12 = cbar(1,2)
    end do

!==========================================================================
! Assemble the "target" values for the fitter. These are just the parameter
! element values followed by the simulated measurements.

    size_y = n_det + n_p_ele
    allocate(x(size_y), y(size_y), sig(size_y))
    allocate(a(n_p_ele), maska(n_p_ele), alpha(n_p_ele,n_p_ele), covar(n_p_ele,n_p_ele))

    iy = 1
    y = 0.
    a = 0.

  ! Load element values
    do i_ele = 1, n_p_ele
       x(iy) = iy
       sig(iy) = 1./sqrt(p_ele(i_ele)%wt)
       maska(iy) = .true.
       iy = iy + 1
    end do

  ! Load measurements
    do i_ele = 1, n_det
       x(iy) = i_ele
       sig(iy) = 1./sqrt(det(i_ele)%wt)
       select case (det(i_ele)%param)
       case(orbit_x$)
          y(iy) = meas(i_ele)%orbit_x
       case(orbit_y$)
          y(iy) = meas(i_ele)%orbit_y
       case(eta_xx$)
          y(iy) = meas(i_ele)%eta_x
       case(eta_yy$)
          y(iy) = meas(i_ele)%eta_y
       case(phi_aa$)
          y(iy) = meas(i_ele)%phi_a
       case(phi_bb$)
          y(iy) = meas(i_ele)%phi_b
       case(cbar12$)
          y(iy) = meas(i_ele)%cbar12
       case default
          write(*,*) "Invalid detector key"
          call err_exit
       end select
       iy = iy + 1
    end do

!==========================================================================
! Make specified number of calls to MRQMIN, then cleanup. Run at least
! n_lm_iterations times, then stop as soon as chisq doesn't change.

  alamda = -1.
!  call super_mrqmin(x,y,sig,a,maska,covar,alpha,chisq,lmfunc,alamda)
  call super_mrqmin(y,sig,a,covar,alpha,chisq,lmfunc,alamda,status,maska)
  iter = 0
  write (*,'(A8,2A12)') "Iter", "Chisq", "D_Chisq"
  write(*,'(i8,es12.4)') iter, chisq
  old_chisq = 999.
  alamda = .1
  do while (iter < correct_ring_params%n_lm_iterations)
     iter = iter + 1
     old_chisq = chisq
!     call super_mrqmin(x,y,sig,a,maska,covar,alpha,chisq,lmfunc,alamda)
     call super_mrqmin(y,sig,a,covar,alpha,chisq,lmfunc,alamda, status, maska)
     if (status /= 0) then
        print*, "*** Singular matrix encountered--nothing more to do with this seed."
        exit
     end if
     d_chisq = (chisq-old_chisq)
     write(*,'(i8,2es12.4)') iter, chisq, d_chisq
  end do
  alamda = 0.
!  call super_mrqmin(x,y,sig,a,maska,covar,alpha,chisq,lmfunc,alamda)
  call super_mrqmin(y,sig,a,covar,alpha,chisq,lmfunc,alamda,status,maska)

!==========================================================================
! Apply correction to misaligned ring, i.e., put in the negative of the
! parameter values calculated by the minimization.

  write(*,*) "Loading correction..."
  do i_ele = 1, n_p_ele
     select case(p_ele(i_ele)%param)
     case(hkick$,vkick$,k1$)
        ring%ele(p_ele(i_ele)%loc)%value(p_ele(i_ele)%param) = &
        ring%ele(p_ele(i_ele)%loc)%value(p_ele(i_ele)%param) - a(i_ele)
     case (-1)
        ring%ele(p_ele(i_ele)%loc)%a_pole(1) = &
        ring%ele(p_ele(i_ele)%loc)%a_pole(1) - a(i_ele)
     case default
        write(*,*) "Invalid parameter key"
        call err_exit
     end select
  end do
  call twiss_and_track(ring, co)
  write(*,*) "Correction loaded."

  ! Clean up
  deallocate(x, y, sig, a, maska, alpha, covar, co)

contains
!==========================================================================
! Two routines for calculating rotated BPM measurements. The first is
! just a simple rotation matrix. The second (for cbar) is more
! complicated, but a complete derivation can be found in my forthcoming
! thesis.

  subroutine rotate_meas(x, y, angle)
    implicit none
    real(rp), intent(inout) :: x, y
    real(rp), intent(in) :: angle
    real(rp) :: x0, y0

    x0 = x
    y0 = y
    x = cos(angle)*x0 - sin(angle)*y0
    y = sin(angle)*x0 + cos(angle)*y0
  end subroutine rotate_meas

  subroutine rotate_cbar(cbar, theta)
    implicit none
    real(rp), intent(inout) :: cbar(2,2)
    real(rp), intent(in) :: theta
    real(rp) oldcbar12, oldcbar22

    ! This is the small angle formula
    oldcbar12 = cbar(1,2)
    oldcbar22 = cbar(2,2)
    cbar(1,2) = oldcbar12 * (1 - 2 * oldcbar22*theta)
    cbar(2,2) = oldcbar22 + (oldcbar12**2 - oldcbar22**2 - 1) * theta
  end subroutine rotate_cbar

end subroutine correct_ring

!==========================================================================
! Routine to provide the values for the chi^2 that MRQMIN will attempt to
! minimize. This is essentially just a wrapper for passing values to
! RING_TO_Y.

subroutine lmfunc(a,y,dyda, status)
  implicit none
  real(rp), dimension(:), intent(in) :: a
  real(rp), dimension(:), intent(out) :: y
  real(rp), dimension(:,:), intent(out) :: dyda

  real(rp), dimension(size(a)) :: a2
  real(rp), dimension(size(y)) :: y1, y2
  real(rp), parameter :: delta = 1.e-7
  integer status
  integer ia

  call ring_to_y(a, y)
  if (any(y==999)) then
     dyda = 999
     return
  end if

  a2 = a
  do ia = 1, size(a)
     a2(ia) = a(ia) + delta
     call ring_to_y(a2, y2)
     ! Central difference
     ! a2(ia) = a(ia) - delta
     ! call ring_to_y(a2, y1)
     ! dyda(:,ia) = (y2-y1)/(2.*delta)

     ! Forward difference
     dyda(:,ia) = (y2-y)/delta

     a2(ia) = a(ia)
  end do
end subroutine lmfunc

!==========================================================================
! Routine used in minimization for calculating how the relevant output
! parameters (Y) change as a function of the input parameters (A).

subroutine ring_to_y(a, y)
  implicit none

  real(rp), dimension(:), intent(in) :: a
  real(rp), dimension(:), intent(out) :: y

  integer ia, iy, ip, i_det
  real(rp) cbar(2,2)
  logical everything_ok

  ! Put the values A into the right places in the ring
  do ia = 1, size(a)
     select case(p_ele(ia)%param)
     case(hkick$, vkick$, k1$)
        cr_model_ring%ele(p_ele(ia)%loc)%value(p_ele(ia)%param) = &
             p_ele(ia)%orig_val + a(ia)
     case(-1)
        cr_model_ring%ele(p_ele(ia)%loc)%a_pole(1) = &
             p_ele(ia)%orig_val + a(ia)
     case default
        write(*,*) "Invalid parameter key"
        call err_exit
     end select
  end do

  ! Recalculate the ring parameters
  if (allocated(cr_model_co)) cr_model_co(0)%vec = 0.
  call twiss_and_track(cr_model_ring, cr_model_co, everything_ok)

  if (.not. everything_ok) then
     y = 999
     return
  end if

  ! Stick everything back into Y
  iy = 1
  ! Load element values
  do ip = 1, n_p_ele
     y(iy) = a(ip)
     iy = iy + 1
  end do
  
  ! Load measurement information
  do i_det = 1, n_det
     select case (det(i_det)%param)
     case(orbit_x$)
        y(iy) = cr_model_co(det(i_det)%loc)%vec(1)
     case(orbit_y$)
        y(iy) = cr_model_co(det(i_det)%loc)%vec(3)
     case(eta_xx$)
        y(iy) = cr_model_ring%ele(det(i_det)%loc)%x%eta
     case(eta_yy$)
        y(iy) = cr_model_ring%ele(det(i_det)%loc)%y%eta
     case(phi_aa$)
        y(iy) = cr_model_ring%ele(det(i_det)%loc)%a%phi
     case(phi_bb$)
        y(iy) = cr_model_ring%ele(det(i_det)%loc)%b%phi
     case(cbar12$)
        call c_to_cbar(cr_model_ring%ele(det(i_det)%loc), cbar)
        y(iy) = cbar(1,2)
     case default
        write(*,*) "Invalid detector key"
        call err_exit
     end select
     iy = iy + 1
  end do

  ! Send this back if it's wanted
  rms_param_ptr = sqrt(dot_product(a, a)/size(a))

end subroutine ring_to_y

end module correct_ring_mod
