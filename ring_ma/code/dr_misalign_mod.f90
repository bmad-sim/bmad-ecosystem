!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Module dr_misalign_mod
!
! This module contains a routine for assigning random misalignments to
! ellements in a lattice, and a structure for defining the misalignment
! parameters
!--------------------------------------------------------------------------
module dr_misalign_mod
  use bmad_struct
  implicit none

  type ma_struct
     integer key
     character*40 mask
     integer param
     real amp
!    logical detector
  end type ma_struct

  type dr_misalign_params_struct
     real(rp) sigma_cutoff, alignment_multiplier
     logical accumulate_errors, tie_dup_ele
  end type dr_misalign_params_struct
  type(dr_misalign_params_struct) dr_misalign_params

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Subroutine dr_misalign_mod(ring, ma_params)
!
! Subroutine to assign random Gaussian errors to elements in a lattice.
! This works by modifying the x_offset, y_offset, etc. parameters in each
! applicable ele_struct. Certain "global" parameters are set in
! dr_misalign_params.
!
! Input:
!   ma_params -- dr_misalign_params_struct(:). List of parameters for
!                  misalignment
!
! Output:
!   ring      -- lat_struct. Ring whose elements will be misaligned
!     %ele(:) -- Modified elements in ring
!--------------------------------------------------------------------------

  subroutine dr_misalign(ring, ma_params)
    use bmad_struct
    use cesr_utils
    use nr

    implicit none
    type(lat_struct), intent(inout), target :: ring
    type(ma_struct), intent(in), target :: ma_params(:)

    integer, parameter :: allowed_params(8) = (/ x_offset$, y_offset$, tilt$, k1$, roll$, &
                                               s_offset$, x_pitch$, y_pitch$/)
    integer, parameter :: allowed_keys(5) = (/ sbend$, quadrupole$, sextupole$, wiggler$, &
                                             marker$ /)

    type(ele_struct), pointer :: ele
    type(ma_struct), pointer :: ma
    integer i_param, i_ele
    real(rp) harvest, multiplier

    ! Check that the parameters make sense
    if (dr_misalign_params%alignment_multiplier < 0.) then
       write(*,*) "dr_misalign: alignment_multiplier must be non-negative"
       call err_exit
    end if
    if (dr_misalign_params%sigma_cutoff <= 0.) then
       write(*,*) "dr_misalign: sigma_cutoff must be positive"
       call err_exit
    end if

    do i_param = 1, size(ma_params)
       if (ma_params(i_param)%amp == 0) cycle

       ma => ma_params(i_param)

       ! Check a few things
       if (.not. any(allowed_keys==ma%key)) then
          write(*,*) "dr_misalign: Key ", ma%key, " not allowed."
          call err_exit
       else if (ma%key == marker$ .and. ma%param > 0) then
          ! Technically, someone might want to do this. Maybe it should just be a warning.
          write(*,*) "dr_misalign: You are misaligning a marker with a positive parameter."
          call err_exit
       else if (.not. any(allowed_params==ma%param)) then
          write(*,*) "dr_misalign: Parameter ", ma%param, " not allowed."
          call err_exit
       else if (ma%key==sbend$ .and. ma%param==tilt$) then
          write(*,*) "dr_misalign: sbends should use roll, not tilt."
          call err_exit
       end if

       ! Loop over elements
       do i_ele = 1, ring%n_ele_track
          ele => ring%ele(i_ele)
          if (ele%key /= ma%key) cycle
          if (.not. (ma%mask == "" .or. match_reg(ele%name, ma%mask))) cycle

          ! Get a random number within the cutoff. If we're tying elements, and this isn't
          ! the first element, and it's the same as last element, then reuse the old numbers.
          if ((.not. dr_misalign_params%tie_dup_ele) .or. &
              (i_ele == 1) .or. (ring%ele(i_ele)%name .ne. ring%ele(i_ele-1)%name)) then
             harvest = 1000
             do while (abs(harvest) >= dr_misalign_params%sigma_cutoff)
                call gasdev(harvest)
             end do
          end if

          ! Don't apply the multiplier to quad strengths or markers
          if (ma%param == k1$ .or. ma%param < 0) then
             multiplier = 1.
          else
             multiplier = dr_misalign_params%alignment_multiplier
          end if

          ! Quads have to "accumulate" because there initial values are nonzero.
          ! This could be fixed by storing the original quad strength.
          if (ma%param == k1$ .or. dr_misalign_params%accumulate_errors) then
             ele%value(ma%param) = ma%amp * harvest * multiplier + ele%value(ma%param)

          else
             ele%value(ma%param) = ma%amp * harvest * multiplier
          end if

       end do
    end do
  end subroutine dr_misalign
end module dr_misalign_mod
