!+
! Subroutine tao_taper_cmd (include_solenoids, uni_names)
!
! Routine to adjust magnet strenghts to counteract the radiation damping sawtooth effect.
!
! Input:
!   include_solenoids   -- logical: Include solenoid ks fields in the scaling?
!   uni_names           -- character(*): Universes to taper.
!-

subroutine tao_taper_cmd(include_solenoids, uni_names)

use tao_interface, dummy => tao_taper_cmd

implicit none

type (tao_universe_pointer_struct), allocatable, target :: unis(:)
type (tao_universe_struct), pointer :: u

integer iu
logical include_solenoids, err
character(*) uni_names
character(*), parameter :: r_name = 'tao_taper_cmd'

!

if (.not. bmad_com%radiation_damping_on) then
  call out_io (s_warn$, r_name, 'bmad_com%radiation_damping_on is set False.', &
                                'The tapering adjustment is fairly useless in this case.')
endif

call tao_pointer_to_universes(uni_names, unis, err); if (err) return

do iu = 1, size(unis)
  u => unis(iu)%u
  call taper_mag_strengths(u%model%lat, u%base%lat, include_solenoids)
  u%calc%lattice = .true.
enddo

end subroutine
