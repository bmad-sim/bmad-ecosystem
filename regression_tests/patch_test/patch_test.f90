!+
! Program patch_test
!
! This program is part of the Bmad regression testing suite.
!-

program patch_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, slave, slave2
type (branch_struct), pointer :: branch
type (coord_struct) :: start_orb, start2_orb, end_orb, end2_orb, end_orb_bs

integer ip, ie
character(40) fmt

!

open (1, file = 'output.now')
fmt = '(3a, 3f20.15, 5x, 3f20.15)'

call bmad_parser ('patch_test.bmad', lat)

!

branch => lat%branch(1)
call init_coord (start_orb, lat%particle_start, branch%ele(0), upstream_end$)

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)

  call track1 (start_orb, ele, branch%param, end_orb)
  call reverse_orbit(end_orb, start2_orb)
  call track1 (start2_orb, ele, branch%param, end2_orb)
  call reverse_orbit(end2_orb, end2_orb)
  write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-BS"  ABS 1E-12     ', end_orb%vec
  write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-BS-dif"  ABS 1E-12 ', end2_orb%vec - start_orb%vec

  if (valid_tracking_method(ele, positron$, runge_kutta$)) then
    end_orb_bs = end_orb
    ele%tracking_method = runge_kutta$

    call track1 (start_orb, ele, branch%param, end_orb)
    call reverse_orbit(end_orb, start2_orb)
    call track1 (start2_orb, ele, branch%param, end2_orb)
    call reverse_orbit(end2_orb, end2_orb)
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RK"  ABS 1E-12     ', end_orb%vec
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RK-dif"  ABS 1E-12 ', end2_orb%vec - start_orb%vec

    write (1, '(3a, f12.6)') '"', trim(ele%name), '-L"  ABS 1E-12 ', ele%value(l$)
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RK-BS"  ABS 1E-12  ', end_orb%vec - end_orb_bs%vec
    write (1, *)

    ele%field_calc = custom$
    call track1 (start_orb, ele, branch%param, end_orb)
    call reverse_orbit(end_orb, start2_orb)
    call track1 (start2_orb, ele, branch%param, end2_orb)
    call reverse_orbit(end2_orb, end2_orb)
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RKC"  ABS 1E-12     ', end_orb%vec
    write (1, '(3a, 6es14.6)') '"', trim(ele%name), '-RKC-dif"  ABS 1E-12 ', end2_orb%vec - start_orb%vec
  endif

  write (1, *)
enddo

!

do ip = 1, 3
  ele => lat%ele(ip)
  if (ele%key == marker$) cycle
  call init_coord (start_orb, lat%particle_start, ele, upstream_end$)
  call track1 (start_orb, ele, lat%param, end_orb)
  write (1, '(3a, 6es14.6)') '"', trim(ele%name), '" ABS 0', end_orb%vec
  if (ele%key == patch$) then
    write (1, '(a, f20.14)') '"L" REL 1E-12 ', ele%value(l$)
  endif
enddo

ele => lat%ele(4)
write (1, '(a, 6es10.2)') '"Flexible" REL 1E-15 ', ele%floor%r, ele%floor%theta, ele%floor%phi, ele%floor%psi

!

branch => lat%branch(2)
do ip = 1, 3
  ele => branch%ele(ip)
  write (1, '(a, i0, a, f20.15)') '"L', ip, '-ref" ABS 1E-15 ', ele%value(l$)
enddo

!

close (1)

!------------------------------------------------------
contains

subroutine reverse_orbit(orb_in, orb_out)

type (coord_struct) orb_in, orb_out

orb_out = orb_in
orb_out%vec(2) = -orb_out%vec(2)
orb_out%vec(4) = -orb_out%vec(4)
orb_out%vec(5) = -orb_out%vec(5)
orb_out%direction = -1
orb_out%species = antiparticle(orb_out%species)

end subroutine

end program
