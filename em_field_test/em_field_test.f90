!+
! Program em_field_test
!
! This program is part of the Bmad regression testing suite.
!-

program em_field_test

use bmad
use nr

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct) :: orb, dorb, orb2
type (em_field_struct) field0, fp, fm, ff
type (em_potential_struct) p0, pp, pm

real(rp) del, ds, dr_ds_kick(11), dr_ds_track(11) 
integer i, ib, ie, j, nargs
logical err, print_extra
character(100) lat_file

namelist / params / del

!

print_extra = .false.
lat_file = 'em_field_test.bmad'

nargs = cesr_iargc()
if (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit

elseif (nargs > 0)then
  call cesr_getarg(1, lat_file)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.
endif

call bmad_parser (lat_file, lat)

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

    call init_coord (orb, lat%beam_start, ele, inside$)
    orb%vec(2) = 0
    orb%vec(4) = 0
    orb%vec(6) = 0

    call em_field_calc (ele, branch%param, orb%vec(5), 0.0_rp, orb, .false., field0, .true., err, p0)

    ff = em_field_struct()

    do i = 1, 3
      j = 2*i - 1
      dorb = orb
      dorb%vec(j) = orb%vec(j) + del
      call em_field_calc (ele, lat%param, dorb%vec(5), 1.0_rp, dorb, .false., fp, .false., err, pp)
      dorb%vec(j) = orb%vec(j) - del
      call em_field_calc (ele, lat%param, dorb%vec(5), 1.0_rp, dorb, .false., fm, .false., err, pm)
      ff%dE(:,i) = (fp%e - fm%e) / (2 * del)
      ff%dB(:,i) = (fp%b - fm%b) / (2 * del)
      select case (i)
      case (1)
        ff%B(2) = ff%B(2) - (pp%A(3) - pm%A(3)) / (2 * del)
        ff%B(3) = ff%B(3) + (pp%A(2) - pm%A(2)) / (2 * del)
      case (2)
        ff%B(3) = ff%B(3) - (pp%A(1) - pm%A(1)) / (2 * del)
        ff%B(1) = ff%B(1) + (pp%A(3) - pm%A(3)) / (2 * del)
      case (3)
        ff%B(1) = ff%B(1) - (pp%A(2) - pm%A(2)) / (2 * del)
        ff%B(2) = ff%B(2) + (pp%A(1) - pm%A(1)) / (2 * del)
      end select
    enddo

    !

    ds = 1d-6
    call twiss_and_track_intra_ele (ele, branch%param, orb%vec(5), orb%vec(5)+ds, .false., .false., orb, orb2)
    dr_ds_track(2) = (orb2%vec(2) - orb%vec(2)) / ds
    dr_ds_track(4) = (orb2%vec(4) - orb%vec(4)) / ds
    call kick_vector_calc (ele, branch%param, orb%vec(5), 1.0_rp, orb, .false., dr_ds_kick, err)

    !

    if (print_extra) then
      print '(a, 3f10.3)', 'At:', orb%vec(1:5:2)

      print *
      print '(2x, a, t49, a, t95, a)', 'dB(theory)', 'dB(actual)', 'dB(theory) - dB(actual)'
      do i = 1, 3
        print '(3(3es14.6, 4x))', field0%dB(i,:), ff%dB(i,:), field0%dB(i,:) - ff%dB(i,:)
      enddo

      print *
      print *, ' B(actual)     B(theory)     B(actual) - B(theory)'
      do i = 1, 3
        print '(3es14.6)', field0%B(i), ff%B(i), field0%B(i) - ff%B(i)
      enddo

      print *
      print '(a, 3es14.6)', 'B:      ', field0%b
      print '(a, 3es14.6)', 'B Grad: ', ff%dB(1,1),  ff%dB(2,2),  ff%dB(3,3)
      print '(a, 3es14.6)', 'B Curl: ', ff%dB(2,3) - ff%dB(3,2), ff%dB(3,1) - ff%dB(1,3), ff%dB(1,2) - ff%dB(2,1)
      print '(a, 3es14.6)', 'B Div:  ', ff%dB(1,1) + ff%dB(2,2) + ff%dB(3,3)
    endif

    write (1, '(3a, 3es16.8)') '"', trim(ele%name), ':B" REL 1E-6', field0%B
    if (any(pp%a /= 0)) write (1, '(3a, 3es16.8)') '"', trim(ele%name), ':B-diff" REL 1E-6', field0%B-ff%B

    write (1, '(3a, 2es16.8)') '"', trim(ele%name), ':Kick" REL 1E-6', dr_ds_kick(2:4:2)
    write (1, '(3a, 2es16.8)') '"', trim(ele%name), ':dKick" REL 1E-6', dr_ds_kick(2:4:2) - dr_ds_track(2:4:2)

    do i = 1, 3
      write (1, '(3a, i0, a, 3es16.8)') '"', trim(ele%name), ':dB Row', i, '" REL 1E-6', field0%dB(i,:)
    enddo

    do i = 1, 3
      write (1, '(3a, i0, a, 3es16.8)') '"', trim(ele%name), ':dB-diff Row', i, '" REL 1E-6', field0%dB(i,:)-ff%dB(i,:)
    enddo

    write (1, *)

  enddo
enddo

!

close (1)

end program
