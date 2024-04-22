!+
! Program coord_test
!
! This program is part of the Bmad regression testing suite.
!-

program coord_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) orbit, orb0
type (branch_struct), pointer :: branch

real(rp) :: vec0(6) = [0.1_rp, 0.02_rp, 0.3_rp, 0.04_rp, 0.5_rp, -0.01_rp]
real(rp) s_out, s_pos, diff_sum(8)

integer ie, nargs
logical debug_mode
character(200) :: lat_file = 'coord_test.bmad'

! 

nargs = command_argument_count()
debug_mode = (nargs > 0)

if (debug_mode) then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
endif

call bmad_parser (lat_file, lat)
bmad_com%spin_tracking_on = .true.

open (1, file = 'output.now')

!

vec0 = lat%particle_start%vec
diff_sum = 0

branch => lat%branch(2)
do ie = 1, branch%n_ele_track-1
  ele => branch%ele(ie)
  call init_coord (orb0, vec0, ele, upstream_end$)
  orb0%t = 0

  call offset_track_drift(orb0, ele,  1,  1, diff_sum)
  call offset_track_drift(orb0, ele, -1,  1, diff_sum)
  call offset_track_drift(orb0, ele,  1, -1, diff_sum)
  call offset_track_drift(orb0, ele, -1, -1, diff_sum)

  call offset_track_static(orb0, ele,  1,  1, diff_sum)
  call offset_track_static(orb0, ele, -1,  1, diff_sum)
  call offset_track_static(orb0, ele,  1, -1, diff_sum)
  call offset_track_static(orb0, ele, -1, -1, diff_sum)
enddo

!

orb0%vec = vec0
orbit = orb0
call canonical_to_angle_coords(orbit)
call angle_to_canonical_coords(orbit)
write (1, '(a, 6es16.8)') '"dAngle" ABS 1E-15', orbit%vec - orb0%vec

!

orb0%vec = lat%particle_start%vec
ele => lat%ele(1)
call init_coord (orbit, orb0%vec, ele, upstream_end$, ele%branch%param%particle)
call offset_particle (ele, set$, orbit)
call offset_particle (ele, unset$, orbit)

write (1, '(a, 6es20.12)') '"orbit1-electron"  ABS 1E-14', orbit%vec
write (1, '(a, 6es20.12)') '"length1-electron" ABS 1E-14', orbit%s, c_light * orbit%t, orbit%dt_ref*c_light


orb0%vec(5) = 0
orb0%vec(6) = sqrt(1 - orb0%vec(2)**2 - orb0%vec(4)**2)

branch => lat%branch(1)
ele => branch%ele(1)
call init_coord (orbit, orb0%vec, ele, upstream_end$, ele%branch%param%particle)
call offset_photon (ele, orbit, set$)
call offset_photon (ele, orbit, unset$)

write (1, '(a, 6es20.12)') '"orbit1-photon"  ABS 1E-14', orbit%vec
write (1, '(a, 6es20.12)') '"length1-photon" ABS 1E-14', orbit%s, c_light * orbit%t, orbit%dt_ref*c_light


ele => branch%ele(2)
orb0%t = 0
call init_coord (orbit, orb0%vec, ele, upstream_end$, ele%branch%param%particle)
call offset_photon (ele, orbit, set$)
call offset_photon (ele, orbit, unset$)

write (1, '(a, 6es20.12)') '"orbit2-photon"  ABS 1E-14', orbit%vec
write (1, '(a, 6es20.12)') '"length2-photon" ABS 1E-14', orbit%s-branch%ele(1)%s, c_light * orbit%t, orbit%dt_ref*c_light

ele => branch%ele(3)
orb0%t = 0
call init_coord (orbit, orb0%vec, ele, upstream_end$, ele%branch%param%particle)
call offset_photon (ele, orbit, set$)
call offset_photon (ele, orbit, unset$)

write (1, '(a, 6es20.12)') '"orbit3-photon"  ABS 1E-14', orbit%vec
write (1, '(a, 6es20.12)') '"length3-photon" ABS 1E-14', orbit%s-branch%ele(1)%s, c_light * orbit%t, orbit%dt_ref*c_light

!

orbit%vec(5) = real_garbage$
orbit%t = 1e-9_rp
call init_coord(orbit, orbit, ele, upstream_end$)
write (1, '(a, 2es20.12)') '"t-set" REL 1E-10', orbit%t, orbit%vec(6)

!

close(1)

!-----------------------------------------------------------------------------
contains

subroutine offset_track_drift (orb0, ele, orient, dir, diff_sum)

type (ele_struct) ele, ele2
type (coord_struct) orb0, orb1, orb2

real(rp) diff_sum(:), s_out
integer orient, dir

character(40) name

!

orb1 = orb0
orb1%direction = dir
orb2 = orb1
ele%orientation = orient
ele2 = ele; call zero_ele_offsets(ele2)
name = trim(ele%name) // ':O:' // int_str(orient) // ':D:' // int_str(dir)
call attribute_bookkeeper(ele2, .true.)

if (debug_mode) print *
call write_orbit(trim(name) // ':Drift-Start', orb1)
call offset_particle (ele, set$, orb1, s_out = s_out)
call write_orbit(trim(name) // ':Drift-In:  ', orb1, s_out)
call track1(orb1, ele2, lat%param, orb1)
call write_orbit(trim(name) // ':Drift-Xfer ', orb1)
call offset_particle (ele, unset$, orb1, s_out = s_out)
call write_orbit(trim(name) // ':Drift-Out  ', orb1, s_out)
if (ele2%key == sbend$) then
  ele2%value(ref_tilt$) = ele%value(ref_tilt$)
  call attribute_bookkeeper(ele2, .true.)
endif
call track1(orb2, ele2, lat%param, orb2)
if (debug_mode) print '(a, 6f13.8, 4x, 2f13.8)', trim(name) // ':NoMis-Track:', orb2%vec, orb2%t*c_light
if (debug_mode) print '(a, 6f13.8, 4x, 2f13.8)', trim(name) // ':Drift-Diff: ', orb1%vec-orb2%vec, (orb1%t-orb2%t)*c_light

diff_sum(1:6) = diff_sum(1:6) + abs(orb1%vec-orb2%vec)
diff_sum(7) = diff_sum(7) + abs(orb1%t-orb2%t)


end subroutine offset_track_drift

!-----------------------------------------------------------------------------
! contains

subroutine offset_track_static (orb0, ele, orient, dir, diff_sum)

type (ele_struct) ele, ele2
type (coord_struct) orb0, orb1, orb2

real(rp) diff_sum(:), s_pos, s_out
integer orient, dir

character(40) name

!

orb1 = orb0
orb1%direction = dir
ele%orientation = orient
ele2 = ele; call zero_ele_offsets(ele2)
name = trim(ele%name) // ':O:' // int_str(orient) // ':D:' // int_str(dir) // ':Static'
s_pos = lat%particle_start%vec(5)

if (debug_mode) print *
call write_orbit(trim(name) // '-Start', orb1, s_pos)
call offset_particle (ele, set$, orb1, drift_to_edge = no$, s_pos = s_pos, s_out = s_out)
call write_orbit(trim(name) // '-In   ', orb1, s_out)
call offset_particle (ele, unset$, orb1, drift_to_edge = no$, s_pos = s_out, s_out = s_out)
call write_orbit(trim(name) // '-Out  ', orb1, s_out)
if (debug_mode) print '(a, 6f13.8, 4x, 2f13.8)', trim(name) // '-Diff: ', orb1%vec-orb0%vec, (orb1%t-orb0%t)*c_light, s_out-s_pos

diff_sum(1:6) = diff_sum(1:6) + abs(orb1%vec-orb0%vec)
diff_sum(7) = diff_sum(7) + abs(orb1%t-orb0%t)*c_light
diff_sum(8) = diff_sum(8) + abs(s_out-s_pos)

end subroutine offset_track_static

!-----------------------------------------------------------------------------
! contains

subroutine write_orbit (name, orbit, s_out)

type (coord_struct) orbit
real(rp), optional :: s_out
character(*) name

!

if (present(s_out)) then
  if (debug_mode) print '(a, 6f13.8, 4x, 2f13.8)', name, orbit%vec, orbit%t*c_light, s_out
  write (1, '(a, t40, a, 6f13.8, 4x, 2f13.8)') quote(trim(name) // '-orb'), 'ABS 1e-8', orbit%vec, orbit%t*c_light, s_out
else
  if (debug_mode) print '(a, 6f13.8, 4x, 2f13.8)', name, orbit%vec, orbit%t*c_light
  write (1, '(a, t40, a, 6f13.8, 4x, 2f13.8)') quote(trim(name) // '-orb'), 'ABS 1e-8', orbit%vec, orbit%t*c_light
endif

write (1, '(a, t40, a, 6f13.8, 4x, 2f13.8)') quote(trim(name) // '-spin'), 'ABS 1e-8', orbit%spin


end subroutine write_orbit

end program
