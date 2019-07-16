program hdf5_test

use beam_mod

implicit none

type (lat_struct) lat
type (bunch_struct) bunch1, bunch2
type (beam_struct) beam
type (pmd_header_struct) pmd_header
type (beam_init_struct) beam_init
integer i, j, n_part
logical error, good1(10), good2(10)

!

open (1, file = 'output.now', recl = 200)

beam_init%position_file = 'beam_1.dat'
n_part = 10

call bmad_parser('lat.bmad', lat)
call init_bunch_distribution(lat%ele(0), lat%param, beam_init, 0, bunch1)
call init_bunch_distribution(lat%ele(0), lat%param, beam_init, 0, bunch2)
do i = 1, n_part
  j = i + 2
  if (j > n_part) j = j - n_part
  bunch2%particle(i) = bunch1%particle(j)
  bunch2%particle(i)%species = photon$
  bunch2%particle(i)%spin = 0
  bunch2%particle(i)%field = [i, i+1]
  bunch2%particle(i)%path_len = 100_rp * (i+5.0_rp)
  bunch2%particle(i)%phase = 1d-3 * [j, j-1]
  bunch2%particle(i)%beta = 1
  bunch2%particle(i)%s = j
enddo

call hdf5_write_beam ('bunch.h5', [bunch1, bunch2], .false., error)
call hdf5_read_beam('bunch.h5', beam, pmd_header, error)

do i = 1, n_part
  good1(i) = is_equal(bunch1%particle(i), beam%bunch(1)%particle(i))
  good2(i) = is_equal(bunch2%particle(i), beam%bunch(2)%particle(i))  
enddo

write (1, '(a, 10l1, a)') '"Bunch1-Match" STR "', good1, '"'
write (1, '(a, 10l1, a)') '"Bunch2-Match" STR "', good2, '"'

close (1)

!-----------------------------------------
contains

function is_equal(p1, p2) result (equal)

type (coord_struct) p1, p2
logical equal

equal = is_eqv(p1%vec, p2%vec, 'vec')
equal = equal .and. is_eqv(p1%spin, p2%spin, 'spin')
equal = equal .and. is_eqv(p1%field, p2%field, 'field')
equal = equal .and. is_eqv(p1%phase, p2%phase, 'phase')
equal = equal .and. is_eq(p1%s, p2%s, 's')
equal = equal .and. is_eq(p1%t, p2%t, 't')
equal = equal .and. is_eq(p1%charge, p2%charge, 'charge')
equal = equal .and. is_eq(p1%path_len, p2%path_len, 'path_len')
equal = equal .and. is_eq(p1%p0c, p2%p0c, 'p0c')
equal = equal .and. is_eq(p1%beta, p2%beta, 'beta')

equal = equal .and. is_eqi(p1%ix_ele, p2%ix_ele, 'ix_ele')
equal = equal .and. is_eqi(p1%ix_branch, p2%ix_branch, 'ix_branch')
equal = equal .and. is_eqi(p1%state, p2%state, 'state')
equal = equal .and. is_eqi(p1%direction, p2%direction, 'direction')
equal = equal .and. is_eqi(p1%species, p2%species, 'species')
equal = equal .and. is_eqi(p1%location, p2%location, 'location')

end function is_equal

!-----------------------------------------
! contains

function is_eq(r1, r2, who) result (equal)
real(rp), intent(in) :: r1, r2
character(*) who
logical equal

equal = (abs(r1-r2) <= 1d-15 * (abs(r1) + abs(r2)))
if (.not. equal) print *, who

end function is_eq

!-----------------------------------------
! contains

function is_eqi(r1, r2, who) result (equal)
integer, intent(in) :: r1, r2
character(*) who
logical equal

equal = (r1 == r2)
if (.not. equal) print *, who

end function is_eqi

!-----------------------------------------
! contains

function is_eqv(r1, r2, who) result (equal)
real(rp), intent(in) :: r1(:), r2(:)
character(*) who
logical equal

equal = all(abs(r1-r2) <= 1d-15 * (abs(r1) + abs(r2)))
if (.not. equal) print *, who

end function is_eqv

end program
