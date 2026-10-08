!+
! Program lcavity_test
!
! Tests for bmad_standard lcavity tracking:
!   1) The transfer matrix agrees with a finite difference of the tracking. 
!      This is done for normal and reversed elements and for forward and backward tracking.
!   2) A reversed lcavity with zero voltage tracks like a drift.
!   3) A particle that is decelerated below its rest mass is marked as lost.
!-

program lcavity_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, ele2
type (branch_struct), pointer :: branch
type (coord_struct) orb0, orb1, orb_drift, orb_p, orb_m
type (coord_struct), allocatable :: orbit(:)

real(rp) mat(6,6), mat_fd(6,6), vec0(6), del
integer i, j, ie, idir, track_state
logical err
character(40) name, str

!

open (1, file = 'output.now')

call bmad_parser ('lcavity_test.bmad', lat)
branch => lat%branch(0)
vec0 = lat%particle_start%vec
del = 1e-6_rp

! Transfer matrix vs finite difference tracking.

do idir = 1, -1, -2
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key /= lcavity$) cycle

    if (idir == 1) then
      str = 'Fwd'
    else
      str = 'Bkwd'
    endif
    name = trim(ele%name) // '-' // trim(str)
    if (ele%orientation == -1) name = trim(ele%name) // '-Rev-' // trim(str)

    call set_orbit(orb0, vec0, ele, idir)
    orb1 = orb0
    call mat_make_unit(mat)
    call track1_bmad (orb1, ele, branch%param, mat6 = mat, make_matrix = .true.)

    do j = 1, 6
      call set_orbit(orb_p, vec0 + del * unit_vec(j), ele, idir)
      call set_orbit(orb_m, vec0 - del * unit_vec(j), ele, idir)
      call track1_bmad (orb_p, ele, branch%param)
      call track1_bmad (orb_m, ele, branch%param)
      mat_fd(:,j) = (orb_p%vec - orb_m%vec) / (2 * del)
    enddo

    write (1, '(3a, 6es20.12)') '"', trim(name), ':Orbit"         ABS 1e-12', orb1%vec
    do i = 1, 6
      write (1, '(3a, i0, a, 6es20.12)') '"', trim(name), ':Mat-Row', i, '"     ABS 1e-10', mat(i,:)
    enddo
    write (1, '(3a, es12.2)') '"', trim(name), ':Mat-FD-Diff"    ABS 1e-7 ', maxval(abs(mat - mat_fd))
  enddo
enddo

! Reversed zero voltage lcavity must track like a drift.
! Note: branch%ele(branch%n_ele_track) is the end marker.

ele => branch%ele(branch%n_ele_track-2)       ! Reversed v0
ele2 => branch%ele(branch%n_ele_track-1)      ! Reversed d0
call set_orbit(orb1, vec0, ele, 1)
call set_orbit(orb_drift, vec0, ele2, 1)
call track1 (orb1, ele, branch%param, orb1)
call track1 (orb_drift, ele2, branch%param, orb_drift)
write (1, '(a, 6es12.2)') '"V0-Rev:Drift-Diff" ABS 1e-14', orb1%vec - orb_drift%vec

! Decelerate below the rest mass.

branch => lat%branch(1)
ele => branch%ele(1)

call init_coord (orb0, [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], ele, upstream_end$)
call track1 (orb0, ele, branch%param, orb1)
write (1, '(3a)')          '"Decel-Ref:State"   STR  "', trim(coord_state_name(orb1%state)), '"'

call init_coord (orb0, [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, c_light / (2 * ele%value(rf_frequency$)), 0.0_rp], ele, upstream_end$)
call track1 (orb0, ele, branch%param, orb1)
write (1, '(3a)')          '"Decel-Lost:State"  STR  "', trim(coord_state_name(orb1%state)), '"'
write (1, '(a, l1, a)')    '"Decel-Lost:Has-NaN" STR  "', any(orb1%vec /= orb1%vec), '"'

orb1 = orb0
call mat_make_unit(mat)
call track1_bmad (orb1, ele, branch%param, mat6 = mat, make_matrix = .true.)
write (1, '(3a)')          '"Decel-Lost-Mat:State"  STR  "', trim(coord_state_name(orb1%state)), '"'
write (1, '(a, l1, a)')    '"Decel-Lost-Mat:Has-NaN" STR  "', any(orb1%vec /= orb1%vec) .or. any(mat /= mat), '"'

close (1)

!--------------------------------------------------------------------
contains

subroutine set_orbit(orb, vec, ele, dir)

type (coord_struct) orb
type (ele_struct) ele
real(rp) vec(6)
integer dir

if (dir == 1) then
  call init_coord(orb, vec, ele, upstream_end$, direction = dir)
else
  call init_coord(orb, vec, ele, downstream_end$, direction = dir)
endif

end subroutine set_orbit

!--------------------------------------------------------------------
! contains

function unit_vec(j) result (u)
integer j
real(rp) u(6)
u = 0
u(j) = 1
end function unit_vec

end program
