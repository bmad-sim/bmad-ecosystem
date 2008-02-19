!+
! Subroutine twiss_propagate_all (lat, set_match)
!
! Subroutine to propagate the twiss parameters from the start to the end.
!
! Note: lat%ele(:)%a Twiss parameters are associated with the "A" mode and
! the lat%ele(:)%b Twiss parameters are associated with the "B" mode.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat%ele(0) -- lat_struct: Twiss parameters at the start
!   set_match    -- Logical, optional: If True then when the routine gets 
!                     to a Match element the Twiss values at the end of 
!                     the previous element are transfered to the starting 
!                     Twiss attributes of the Match element. This ensures
!                     that the Twiss values at the end of the Match element
!                     will be the same as the ending Twiss attributes of the
!                     Match element. Default is False.
!   bmad_status  -- Common block status structure:
!       %type_out      -- If True then will type a message if the modes are flipped.
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   lat         -- lat_struct: Lat
!   bmad_status  -- Common block status structure:
!       %ok         -- Set False if an input beta is zero. True otherwise
!-

#include "CESR_platform.inc"

subroutine twiss_propagate_all (lat, set_match)

  use bmad_struct
  use bmad_interface, except_dummy => twiss_propagate_all

  implicit none

  type (lat_struct)  lat

  integer n, n_use
  logical, optional :: set_match
  logical do_set

! Propagate twiss

  n_use = lat%n_ele_track

  bmad_status%ok = .true.
  do_set = logic_option(.false., set_match)

  do n = 1, n_use

    if (do_set .and. lat%ele(n)%key == match$) then
      lat%ele(n)%value(beta_a0$)  = lat%ele(n-1)%a%beta
      lat%ele(n)%value(beta_b0$)  = lat%ele(n-1)%b%beta
      lat%ele(n)%value(alpha_a0$) = lat%ele(n-1)%a%alpha
      lat%ele(n)%value(alpha_b0$) = lat%ele(n-1)%b%alpha
      lat%ele(n)%value(eta_a0$)   = lat%ele(n-1)%a%eta
      lat%ele(n)%value(eta_b0$)   = lat%ele(n-1)%b%eta
      lat%ele(n)%value(etap_a0$)  = lat%ele(n-1)%a%etap
      lat%ele(n)%value(etap_b0$)  = lat%ele(n-1)%b%etap
    endif

    call twiss_propagate1 (lat%ele(n-1), lat%ele(n))
    if (.not. bmad_status%ok) return

  enddo

! Make sure final mode is same as initial mode

  if (lat%param%lattice_type == circular_lattice$) then
    if (lat%ele(0)%mode_flip .neqv. lat%ele(n_use)%mode_flip) then
      call do_mode_flip (lat%ele(n_use))
      if (bmad_status%type_out .and. .not. bmad_status%ok) then
        print *, 'ERROR IN TWISS_PROPAGATE_ALL: CANNOT MAKE FINAL FLIP STATE'
        print *, '      EQUAL TO THE INITIAL'
      endif
    endif
  endif

end subroutine
