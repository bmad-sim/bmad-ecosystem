program tune_scan

use ts_mod

implicit none

type (ts_params_struct) ts
type (ts_com_struct) ts_com
type (ts_data_struct), allocatable, target :: ts_dat(:,:,:)

real(rp) del_time, time0
integer ja, jb, jz

!---------------------------------------------
! Init

call ts_init_params (ts, ts_com)
allocate (ts_dat(0:ts_com%n_a, 0:ts_com%n_b, 0:ts_com%n_z))
call run_timer ('START')
time0 = 0

!---------------------------------------------
! Main loop

do ja = 0, ts_com%n_a
do jb = 0, ts_com%n_b
do jz = 0, ts_com%n_z
  call ts_track_particle (ts, ts_com, [ja, jb, jz], ts_dat(ja,jb,jz))
  call run_timer ('READ', del_time)

  if (del_time - time0 > ts%timer_print_dtime) then
    print '(a, f10.2, a, 3i5)', 'Time (min): ', del_time/60, '  At point ja, jb, jz = ', ja, jb, jz
    time0 = del_time
  endif
enddo
enddo
enddo

!---------------------------------------------
! Write results

call ts_write_results (ts, ts_com, ts_dat)

end program
