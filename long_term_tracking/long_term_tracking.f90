program long_term_tracking

use lt_tracking_mod

implicit none

type (ltt_params_struct) lttp
type (beam_init_struct) beam_init
type (ltt_com_struct), target :: ltt_com
type (bunch_struct), pointer :: bunch
real(rp) del_time
integer i, iu
character(4) prefix

!

call ltt_init_params(lttp, ltt_com, beam_init)
call ltt_init_tracking (lttp, ltt_com)
call ltt_print_inital_info (lttp, ltt_com)

call run_timer ('START')

select case (lttp%simulation_mode)
case ('CHECK');  call ltt_run_check_mode(lttp, ltt_com, beam_init)  ! A single turn tracking check
case ('SINGLE'); call ltt_run_single_mode(lttp, ltt_com, beam_init) ! Single particle tracking
case ('BEAM');   call ltt_run_beam_mode(lttp, ltt_com, beam_init)   ! Beam tracking
case ('STAT');   call ltt_run_stat_mode(lttp, ltt_com)              ! Lattice statistics (radiation integrals, etc.).
case default
  print *, 'BAD SIMULATION_MODE: ' // lttp%simulation_mode
end select

call run_timer ('READ', del_time)
call ltt_write_master('Tracking time (min): ' // real_str(del_time/60, 4, 2), lttp, append = .true.)

! Regression output is used in conjunction with the Bmad regression test suite.

if (lttp%regression_test) then
  prefix = ltt_com%master_input_file(1:4)
  open (1, file = prefix // '.dat')
  bunch => ltt_com%bunch
  do i = 1, size(bunch%particle)
    write (1, '(3a, i0, a, 6es16.8)') '"', prefix, '-vec-', i, '" ABS 1e-10 ', bunch%particle(i)%vec
    write (1, '(3a, i0, a, 3f16.10)') '"', prefix, '-spin-', i, '" ABS 1e-10 ', bunch%particle(i)%spin
  enddo
  close (1)
endif  

end program

