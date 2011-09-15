!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Module sim_bpm_mod
!
! This module contains subroutines to simulate: position resolution in terms 
! of additive button signal; physical misalignment and gain errors on BPMs; 
! and applying random resolution, gain errors, tilt, shear, etc. to a 
! measurement.
!--------------------------------------------------------------------------
module sim_bpm_mod

  use bmad
  use sim_utils
  use random_mod
  use cesr_basic_mod
  use nonlin_bpm_mod

  implicit none
 
  type det_struct
     real(rp) :: vec(6) ! 6-vector of positions and momenta
     real(rp) :: amp(4) = 0. ! button signals
     real(rp) :: gain(4) = 1.
     real(rp) :: timing(4) = 0.
     integer :: ix_db = 0
     integer :: ix_lat = 0
     real(rp) :: shear_x = 0
     real(rp) :: butn_res_sigma = 0
     real(rp) :: button_offset(4) = 0
     real(rp) :: tilt = 0.
     real(rp) :: x_offset = 0.
     real(rp) :: y_offset = 0.
     real(rp)    wt, param ! compatability with ring_ma
  end type det_struct
  
  type det_error_struct
     real(rp) :: bpm_offset = 0.
     real(rp) :: bpm_rotation = 0.
     real(rp) :: shear_x = 0.
     real(rp) :: gain_sigma = 0.
     real(rp) :: timing_sigma = 0.
  end type det_error_struct

!---------------------------------------------

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine find_bpms_old(cesr, bpm)
! OBSOLETE - use "find_bpms" below.
!
! From cesr_struct, locate BPMs of interest. 
!
! Input:
!   cesr -- cesr_struct
!
! Output:
!   bpm(1:n_bpms) -- det_struct(:); Holds information about BPMs of interest.
!--------------------------------------------------------------------------
  subroutine find_bpms_old(cesr, bpm)
    
    implicit none
    type(cesr_struct) :: cesr
    type(det_struct), allocatable :: bpm(:)
    integer :: i, bpm_index, n_bpms = 0
    
    do i=0, ubound(cesr%det,1)
       if (match_reg(cesr%det(i)%name, 'DUMMY') .or. match_reg(cesr%det(i)%name,'[aA]')) cycle
       if(cesr%det(i)%ix_db .gt. 100) cycle
       n_bpms = n_bpms + 1 ! index starts at zero, but we want the element before the first dummmy
    enddo
    
    allocate(bpm(1:n_bpms))
    bpm_index = 1
    
    do i=0, ubound(cesr%det,1)
       if (match_reg(cesr%det(i)%name, 'DUMMY')  .or. match_reg(cesr%det(i)%name,'[aA]')) cycle
       if (cesr%det(i)%ix_db .gt. 100) cycle
       bpm(bpm_index)%ix_lat = cesr%det(i)%ix_lat ! det() index starts at zero
       bpm(bpm_index)%ix_db = cesr%det(i)%ix_db
       bpm_index = bpm_index + 1
    enddo
    
  end subroutine find_bpms_old
  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine find_bpms(ring, bpm_mask, bpm)
!
! From lat_struct, locate BPMs of interest. 
!
! Input:
!      lat - lat_struct
! bpm_mask - char(40); regex mask for locating BPMs.
!
! Output:
!   bpm(1:n_bpms) -- det_struct(:); Holds information about BPMs of interest.
!--------------------------------------------------------------------------
 subroutine find_bpms(ring, bpm_mask, bpm)

  use cesr_basic_mod
   
   implicit none
   type(cesr_struct) cesr
   type(lat_struct) :: ring
   type(det_struct), allocatable :: bpm(:)
   integer :: i, ix, jx, bpm_index, n_bpms = 0
   character(40) bpm_mask
   
   do i=1, ring%n_ele_max
      if (match_reg(ring%ele(i)%name, trim(bpm_mask))) n_bpms = n_bpms + 1 
   enddo
   
   allocate(bpm(1:n_bpms))

   bpm_index = 1
   
   ! note: bpm(:)%ix_db not populated in this routine anymore.
   bpm(:)%ix_db = -1

   do i=1, ring%n_ele_max
      if (match_reg(ring%ele(i)%name, trim(bpm_mask))) then
         bpm(bpm_index)%ix_lat = ring%ele(i)%ix_ele
         bpm_index = bpm_index + 1
      endif
   enddo
  

  ! if CESR-type lattice, assign bpm(:)%ix_db indices
  if (match_reg(ring%name, "CESR")) then
     call bmad_to_cesr(ring, cesr)
     jx = 1
     do i = 1, ubound(cesr%det,1)
        if (match_reg(cesr%det(i)%name, 'DUMMY')  .or. match_reg(cesr%det(i)%name,'[aA]')) cycle
        if (cesr%det(i)%ix_db .gt. 100) cycle
        bpm(jx)%ix_db = cesr%det(i)%ix_db
        jx = jx + 1
     enddo
  else ! non-CESR-type lattice
     do ix = 1, n_bpms
        bpm(ix)%ix_db = ix
     enddo
  endif
 
 end subroutine find_bpms
 


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine resolution_to_button_error(bpm, current, res)
!
! Determine the desired resolution in terms of an additive error on button 
! signal amplitudes.
!
! Input:
!   bpm -- det_struct, holds information about a single BPM
!   current -- real(rp), in "nonlin_bpm" units
!   res -- real(rp). resolution, in meters
!
! Output:
!   bpm
!      %butn_res_sigma -- resolution in terms of additive button signal error
!--------------------------------------------------------------------------
  subroutine resolution_to_button_error(bpm, current, res)

    implicit none

    type(det_struct) :: bpm
    real(rp), intent(in) :: res, current
    real(rp) x, y, b1(4), b2(4), d(4,0:2,0:2)
    real(rp) :: butn_res_sigma = 0
    integer i


    call nonlin_bpm_set_pointers(25) !use BPM #25 as an example
    x = 1.e-3 !1mm orbit in x and y
    y = 1.e-3
    call nonlin_bpm_interpolate((/x,y,current/),d)
    b1(1:4) = d(1:4,0,0)
    
    ! Now add 1\sigma to the x and y positions, see what the effect is on the button signal amplitudes
    
    x = 1.e-3 + res
    y = 1.e-3 + res
    call nonlin_bpm_interpolate((/x,y,current/),d) 
    b2(1:4) = d(1:4,0,0)
    
    ! butn_res_sigma initialized to zero already
    
    do i=1,4
       butn_res_sigma = butn_res_sigma + (b2(i) - b1(i))**2/4
    enddo
    butn_res_sigma = sqrt(butn_res_sigma)
    bpm%butn_res_sigma = butn_res_sigma
    
  end subroutine resolution_to_button_error
  
  !------------------------------------------


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine bpm_errors(this_bpm, bpm_error_sigmas)
!
! Subroutine to assign random Gaussian errors to a BPM. All supplied sigmas are real(rp).
!
! Input:
!  this_bpm         -- det_struct: holds information about misalignments specific to that BPM.
!  bpm_error_sigmas -- structure holding all BPM misalignment sigmas
!     %bpm_offset   -- Sigma for random horizontal/vertical offsets (meters)
!     %bpm_rotation -- Sigma for BPM tilts (radians)
!     %shear_x      -- Sigma for a "shear" effect where top and bottom BPM button blocks are 
!                         transversely misaligned (in x). (meters)
!     %gain_sigma   -- Sigma for single-button gain errors. (unitless: 0.01 = 1%)
!     %timing_sigma -- Sigma for timing errors (seconds: 10.0e-12)
!
! Output:
!   this_bpm         -- Now with a unique, Gaussian set of misalignements/errors:
!      %x_offset     -- 
!      %y_offset     -- 
!      %tilt         -- 
!      %shear_x      --
!--------------------------------------------------------------------------
  subroutine bpm_errors(this_bpm, bpm_error_sigmas)
        
    implicit none

    type(det_struct) :: this_bpm
    type(det_error_struct) :: bpm_error_sigmas
    real(rp) :: harvest
    integer jx
    
    ! 2011.05.23 - should we incorporate a sigma_cutoff like ring_ma? -JSh

    ! apply tilt first
    call ran_gauss(harvest)
    !this_bpm%tilt = this_bpm%tilt + harvest *  bpm_error_sigmas%bpm_rotation
    this_bpm%tilt = harvest *  bpm_error_sigmas%bpm_rotation


    call ran_gauss(harvest)
    !this_bpm%x_offset = this_bpm%x_offset + harvest * bpm_error_sigmas%bpm_offset
    this_bpm%x_offset = harvest * bpm_error_sigmas%bpm_offset
    call ran_gauss(harvest)
    !this_bpm%y_offset = this_bpm%y_offset + harvest *  bpm_error_sigmas%bpm_offset
    this_bpm%y_offset = harvest *  bpm_error_sigmas%bpm_offset

    
    ! shear is in the rotated coordinate system, so this is ok
    call ran_gauss(harvest)
    !this_bpm%shear_x = this_bpm%shear_x + harvest *  bpm_error_sigmas%shear_x  
    this_bpm%shear_x = harvest *  bpm_error_sigmas%shear_x  
 
    do jx = 1,4 ! loop over buttons at one BPM
       call ran_gauss(harvest)
       !this_bpm%gain(jx) = this_bpm%gain(jx) + harvest * bpm_error_sigmas%gain_sigma
       this_bpm%gain(jx) = 1. + harvest * bpm_error_sigmas%gain_sigma
       call ran_gauss(harvest)  
       !this_bpm%timing(jx) = this_bpm%timing(jx) + harvest * bpm_error_sigmas%timing_sigma
       this_bpm%timing(jx) = harvest * bpm_error_sigmas%timing_sigma
    enddo
    
    
  end subroutine bpm_errors

  !------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine rotate_meas(x,y,angle)
!
! Subroutine to rotate a measurement.
!
! Input:
!   x,y   -- Horizontal and vertical position
!   angle -- angle of rotation
!
! Output:
!   x,y   -- Rotated measurements
!--------------------------------------------------------------------------
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

  
!---------------------------------------------


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine apply_bpm_errors(this_bpm, current)
!
! Subroutine to assign previously-determined random Gaussian errors to 
! BPM data.
!
! Input:
!   this_bpm -- det_struct; holds misalignment information
!   current  -- real(rp); in "nonlin_bpm" units
!
! Output:
!   this_bpm  -- det_struct, now with applied errors.
!     %vec(1)  -- horizontal orbit, after BPM errors are applied
!     %vec(3)  -- vertical orbit, after BPM errors have been applied
!--------------------------------------------------------------------------
  subroutine apply_bpm_errors(this_bpm, current)
    implicit none
    real(rp) :: harvest
    real(rp), intent(in) :: current
    real(rp) :: x_offset, y_offset
    real(rp) :: g(4), b(4), d(4,0:2,0:2), x2(3),x3(3), &
         x, y, butn_res_sigma, shear_x, temp
    integer bpm_ix,kx
    type(det_struct) :: this_bpm
    
    real(rp) timing_factor(4)

    butn_res_sigma = this_bpm%butn_res_sigma
    shear_x = this_bpm%shear_x
    g(1:4) = this_bpm%gain(1:4)
    bpm_ix = this_bpm%ix_db
    b(1:4) = this_bpm%amp(1:4)


    ! 2011.07.27 - 
    ! When errors are being applied, we apply tilts before offsets. Therefore, 
    ! we must apply the offsets in the "un-rotated" coordinate system.
    x_offset = this_bpm%x_offset
    y_offset = this_bpm%y_offset
    call rotate_meas(x_offset,y_offset,-1.*this_bpm%tilt)

    x = this_bpm%vec(1) + this_bpm%x_offset
    y = this_bpm%vec(3) + this_bpm%y_offset

    call bpm_time_to_gain(this_bpm, timing_factor)

! add shear in 2 steps-- buttons 1 and 2, simulate shifting to the right (therefore, the BEAM
! would shift to the LEFT)
!NOTE: this is not entirely equivalent to a button shift! this is actually shifting the BEAM,
!therefore the image charges will not be the same!
!buttons 3/4, simulate shifting to the left: (ie, beam shifts to the right)


    call nonlin_bpm_set_pointers(bpm_ix)

    x2(1)=x-shear_x
    x2(2)=y
    x2(3)=current

    call nonlin_bpm_interpolate(x2,d)
    b(1) = g(1)*d(1,0,0)*timing_factor(1)
    b(2) = g(2)*d(2,0,0)*timing_factor(2)

    x2(1)=x+shear_x
    x2(2)=y
    x2(3)=current

    call nonlin_bpm_interpolate(x2,d)
    b(3) = g(3)*d(3,0,0)*timing_factor(3)
    b(4) = g(4)*d(4,0,0)*timing_factor(4)
    

    do kx=1,4
       call ran_gauss(harvest)
       b(kx) = b(kx) + harvest*butn_res_sigma
    enddo

    call nonlin_orbit(-1,b,x2) ! disable reading BPM offsets
    call nonlin_orbit(bpm_ix,b,x2)

    this_bpm%vec(1) = x2(1)
    this_bpm%vec(3) = x2(2)
    this_bpm%amp(1:4) = b(1:4)
  
  end subroutine apply_bpm_errors

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine bpm_time_to_gain(this_bpm, timing_factor)
!
! Subroutine to convert the timing errors to a gain factor
!
! Input:
!   this_bpm -- det_struct; holds misalignment information
!
! Output:
!   timing_factor -- real(rp); holds the gain factors for 4 buttons
!--------------------------------------------------------------------------
subroutine bpm_time_to_gain(this_bpm, timing_factor)

    type(det_struct) this_bpm
    real(rp) timing_factor(4)
    integer j

    real(rp) z_offset, t_offset, timing_err(4)
    real(rp) no_err_signal, err_signal

!   quadratic fit parameters from an actual fit of the button pulse
!   ( f(t)=a0*t**2+a1*t+a2 )
    real(rp) :: a0= -14.1676, a1= 7586.69, a2= -993141.0

    z_offset=this_bpm%vec(5)
    t_offset=z_offset / c_light

    timing_err(:)=this_bpm%timing(:)+t_offset
    no_err_signal= a2-a1**2/a0/4

    do j=1,4
       err_signal=a0*(timing_err(j)/10.0e-12)**2+a2-a1**2/a0/4
       timing_factor(j)=err_signal/no_err_signal
    enddo

end subroutine bpm_time_to_gain


end module sim_bpm_mod
