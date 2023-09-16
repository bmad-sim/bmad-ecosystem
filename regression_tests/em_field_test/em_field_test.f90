!+
! Program em_field_test
!
! This program is part of the Bmad regression testing suite.
!-

program em_field_test

use bmad
use runge_kutta_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct) :: orb, dorb, orb2
type (em_field_struct) field0, fp, fm, ff, ftp, ftm

real(rp) del, rf_time, ds, dt, dr_ds_kick(11), dr_ds_track(11), dE_dt(3), dB_dt(3)
integer i, ib, ie, j, nargs
logical err, print_extra, rf_on

character(100) lat_file

namelist / params / del, rf_time

!

print_extra = .false.
lat_file = 'em_field_test.bmad'

nargs = command_argument_count()
if (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit

elseif (nargs > 0)then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.
endif

call bmad_parser (lat_file, lat)

rf_time = 1.0
open (1, file = lat_file)
read (1, nml = params)
close (1)

!

open (1, file = 'output.now')

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key == marker$) cycle

    if (ele%key == sad_mult$) then
      print *, 'SKIPPING SAD_MULT!'
      cycle
    endif

    call init_coord (orb, lat%particle_start, ele, inside$)
    orb%vec(2) = 0
    orb%vec(4) = 0
    orb%vec(6) = 0

    call em_field_calc (ele, branch%param, orb%vec(5), orb, .false., field0, .true., err, .true., rf_time = rf_time)

    ff = em_field_struct()

    do i = 1, 3
      if (i == 3) then
        call em_field_calc (ele, lat%param, orb%vec(5)+del, orb, .false., fp, .false., err, .true., rf_time = rf_time)
        call em_field_calc (ele, lat%param, orb%vec(5)-del, orb, .false., fm, .false., err, .true., rf_time = rf_time)
      else
        j = 2*i - 1
        dorb = orb
        dorb%vec(j) = orb%vec(j) + del
        call em_field_calc (ele, lat%param, dorb%vec(5), dorb, .false., fp, .false., err, .true., rf_time = rf_time)
        dorb%vec(j) = orb%vec(j) - del
        call em_field_calc (ele, lat%param, dorb%vec(5), dorb, .false., fm, .false., err, .true., rf_time = rf_time)
      endif

      ff%dE(:,i) = (fp%e - fm%e) / (2 * del)
      ff%dB(:,i) = (fp%b - fm%b) / (2 * del)
      select case (i)
      case (1)
        ff%B(2) = ff%B(2) - (fp%A(3) - fm%A(3)) / (2 * del)
        ff%B(3) = ff%B(3) + (fp%A(2) - fm%A(2)) / (2 * del)
      case (2)
        ff%B(3) = ff%B(3) - (fp%A(1) - fm%A(1)) / (2 * del)
        ff%B(1) = ff%B(1) + (fp%A(3) - fm%A(3)) / (2 * del)
      case (3)
        ff%B(1) = ff%B(1) - (fp%A(2) - fm%A(2)) / (2 * del)
        ff%B(2) = ff%B(2) + (fp%A(1) - fm%A(1)) / (2 * del)
      end select
    enddo

    rf_on = .false.
    if (ele%key == rfcavity$ .or. ele%key == lcavity$) then
      dt = 1e-4_rp / ele%value(rf_frequency$)
      call em_field_calc (ele, lat%param, dorb%vec(5), orb, .false., ftp, .false., err, .true., rf_time = rf_time+dt)
      call em_field_calc (ele, lat%param, dorb%vec(5), orb, .false., ftm, .false., err, .true., rf_time = rf_time-dt)
      dE_dt = (ftp%E - ftm%E) / (2 * dt) / c_light**2
      dB_dt = (ftp%B - ftm%B) / (2 * dt)
      rf_on = .true.
    endif

    !

    ds = 1d-6
    call twiss_and_track_intra_ele (ele, branch%param, orb%vec(5), orb%vec(5)+ds, .false., .false., orb, orb2)
    dr_ds_track(2) = (orb2%vec(2) - orb%vec(2)) / ds
    dr_ds_track(4) = (orb2%vec(4) - orb%vec(4)) / ds

    orb2 = orb
    call offset_particle (ele, set$, orb2, set_hvkicks = .false., drift_to_edge = no$, s_pos = orb%vec(5))
    call kick_vector_calc (ele, branch%param, orb%vec(5), orb2, dr_ds_kick, err)
    orb2%vec(2:4:2) = dr_ds_kick(2:4:2)
    call tilt_coords (-ele%value(tilt_tot$), orb2%vec)
    dr_ds_kick(2:4:2) = orb2%vec(2:4:2) 

    !

    if (print_extra) then
      print '(a, 3f10.3)', 'At:', orb%vec(1:5:2)

      if (any(field0%E /= 0)) then
        print *
        print '(2x, a, t49, a, t95, a)', 'dE(em_field)', 'dE(diff)', 'dE(em_field) - dE(diff)'
        do i = 1, 3
          print '(3(3es14.6, 4x))', field0%dE(i,:), ff%dE(i,:), field0%dE(i,:) - ff%dE(i,:)
        enddo
      endif

      print *
      print '(2x, a, t49, a, t95, a)', 'dB(em_field)', 'dB(diff)', 'dB(em_field) - dB(diff)'
      do i = 1, 3
        print '(3(3es14.6, 4x))', field0%dB(i,:), ff%dB(i,:), field0%dB(i,:) - ff%dB(i,:)
      enddo

      print *
      print *, ' B(em_field)   B(dA)         B(em_field) - B(dA)'
      do i = 1, 3
        print '(3es14.6)', field0%B(i), ff%B(i), field0%B(i) - ff%B(i)
      enddo

      if (any(field0%E /= 0)) then
        print *
        print *, 'From em_field and taking differences:'
        print '(a, 3es14.6)', 'E:            ', field0%E
        print '(a, 3es14.6)', 'E Grad:       ', ff%dE(1,1),  ff%dE(2,2),  ff%dE(3,3)
        print '(a, 3es14.6)', 'E Curl:       ', ff%dE(2,3) - ff%dE(3,2),  ff%dE(1,3) - ff%dE(3,1), ff%dE(2,1) - ff%dE(1,2)
        if (rf_on) then
          print '(a, 3es14.6)', 'E Curl+dB/dt: ', ff%dE(3,2) - ff%dE(2,3) + dB_dt(1), &
                                                 ff%dE(1,3) - ff%dE(3,1) + dB_dt(2), ff%dE(2,1) - ff%dE(1,2) + dB_dt(3)
        endif
        print '(a, 3es14.6)', 'E Div:  ', ff%dE(1,1) + ff%dE(2,2) + ff%dE(3,3)
      endif

      print *
      print *, 'From em_field and taking differences:'
      print '(a, 3es14.6)', 'B:            ', field0%B
      print '(a, 3es14.6)', 'B Grad:       ', ff%dB(1,1),  ff%dB(2,2),  ff%dB(3,3)
      print '(a, 3es14.6)', 'B Curl:        ', ff%dB(3,2) - ff%dB(2,3),  ff%dB(1,3) - ff%dB(3,1), ff%dB(2,1) - ff%dB(1,2)
      if (rf_on) then
        print '(a, 3es14.6)', 'B Curl-dE/dt: ', ff%dB(3,2) - ff%dB(2,3) - dE_dt(1), &
                                                 ff%dB(1,3) - ff%dB(3,1) - dE_dt(2), ff%dB(2,1) - ff%dB(1,2) - dE_dt(3)
      endif
      print '(a, 3es14.6)', 'B Div:        ', ff%dB(1,1) + ff%dB(2,2) + ff%dB(3,3)
    endif

    write (1, '(3a, 3es16.8)') '"', trim(ele%name), ':B" REL 1E-6', field0%B
    if (any(fp%a /= 0)) write (1, '(3a, 3es16.8)') '"', trim(ele%name), ':B-diff" REL 1E-6', field0%B-ff%B

    write (1, '(3a, 2es16.8)') '"', trim(ele%name), ':Kick" REL 1E-6', dr_ds_kick(2:4:2)
    write (1, '(3a, 2es16.8)') '"', trim(ele%name), ':dKick" REL 1E-6', dr_ds_kick(2:4:2) - dr_ds_track(2:4:2)

    do i = 1, 3
      write (1, '(3a, i0, a, 3es16.8)') '"', trim(ele%name), ':dB Row', i, '" REL 1E-6', field0%dB(i,:)
    enddo

    do i = 1, 3
      write (1, '(3a, i0, a, 3es16.8)') '"', trim(ele%name), ':dB-diff Row', i, '" REL 1E-6', field0%dB(i,:)-ff%dB(i,:)
    enddo

    if (any(field0%E /= 0)) then
      do i = 1, 3
        write (1, '(3a, i0, a, 3es16.8)') '"', trim(ele%name), ':dE Row', i, '" REL 1E-6', field0%dE(i,:)
      enddo

      do i = 1, 3
        write (1, '(3a, i0, a, 3es16.8)') '"', trim(ele%name), ':dE-diff Row', i, '" REL 1E-6', field0%dE(i,:)-ff%dE(i,:)
      enddo
    endif


    write (1, *)

  enddo
enddo

!

close (1)

end program
