!+
! Subroutine calc_emit_from_beam_init (beam_init, ele, species)
!
! Routine to calculate the emittances from the beam_init structure.
!
! Input:
!   beam_init   -- beam_init_struct: 
!   ele         -- ele_struct:
!   param       -- lat_param_struct:
!
! Ouput:
!   ele%a%emit  -- Real(rp): a emittance
!   ele%b%emit  -- Real(rp): b emittance
!-

subroutine calc_emit_from_beam_init (beam_init, ele, species)

use bmad_routine_interface, dummy => calc_emit_from_beam_init

implicit none

type (beam_init_struct) beam_init
type (ele_struct) ele

real(rp) ran_g(2)
integer species

character(16) :: r_name = 'calc_emit_from_beam_init'

! Check

if ((beam_init%a_norm_emit /= 0 .and. beam_init%a_emit /= 0) .or. &
    (beam_init%b_norm_emit /= 0 .and. beam_init%b_emit /= 0)) then
  call out_io (s_fatal$, r_name, 'SETTING BOTH NORM_EMIT AND EMIT IN BEAM_INIT STRUCTURE IS NOT ALLOWED.')
  if (global_com%exit_on_error) call err_exit
endif

!

if (beam_init%a_norm_emit /= 0  .and. ele%value(p0c$) /= 0) then
  ele%a%emit = beam_init%a_norm_emit * mass_of(species) / ele%value(p0c$)
else
  ele%a%emit = beam_init%a_emit
endif

if (beam_init%b_norm_emit /= 0 .and. ele%value(p0c$) /= 0) then
  ele%b%emit = beam_init%b_norm_emit * mass_of(species) / ele%value(p0c$)
else
  ele%b%emit = beam_init%b_emit 
endif

! Add jitter if needed

if (any(beam_init%emit_jitter /= 0)) then
  call ran_gauss(ran_g) ! ran(3:4) for z and e jitter used below
  ele%a%emit = ele%a%emit * (1 + beam_init%emit_jitter(1) * ran_g(1))
  ele%b%emit = ele%b%emit * (1 + beam_init%emit_jitter(2) * ran_g(2))
endif

end subroutine calc_emit_from_beam_init 

