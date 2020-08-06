program long_term_tracking

use lt_tracking_mod

implicit none

type (ltt_params_struct) lttp
type (beam_init_struct) beam_init
type (ltt_com_struct) ltt_com

real(rp) del_time

!

call ltt_init_params(lttp, ltt_com, beam_init)
call ltt_init_tracking (lttp, ltt_com)
call ltt_print_inital_info (lttp, ltt_com)

call run_timer ('START')

select case (lttp%simulation_mode)
case ('CHECK');  call ltt_run_check_mode(lttp, ltt_com, beam_init)  ! A single turn tracking check
case ('SINGLE'); call ltt_run_single_mode(lttp, ltt_com, beam_init) ! Single particle tracking
case ('BUNCH');  call ltt_run_bunch_mode(lttp, ltt_com, beam_init)  ! Beam tracking
case ('STAT');   call ltt_run_stat_mode(lttp, ltt_com)              ! Lattice statistics (radiation integrals, etc.).
case default
  print *, 'BAD SIMULATION_MODE: ' // lttp%simulation_mode
end select

call run_timer ('READ', del_time)
print '(a, f8.2)', 'Tracking time (min)', del_time/60

end program

