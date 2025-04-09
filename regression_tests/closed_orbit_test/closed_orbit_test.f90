program co_test

use bmad

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: orb_fwd(:), orb_back(:)

real(rp) :: vec0(6) = 0, df(6), db(6)
integer i, n, n_dim

!

open (1, file = 'output.now')

bmad_com%auto_bookkeeper = .false.
call bmad_parser ('bmad_L9A18A000-_MOVEREC.lat', lat)

n = lat%n_ele_track
call reallocate_coord(orb_fwd, lat)
call reallocate_coord(orb_back, lat)
call init_coord (orb_back(n), vec0, lat%ele(n), downstream_end$, antiparticle(lat%param%particle), -1)

call set_on_off (rfcavity$, lat, off$)

do n_dim = 4, 6
  if (n_dim == 6) then
    bmad_com%radiation_damping_on = .true.
    call set_on_off (rfcavity$, lat, on$)
  endif

  call closed_orbit_calc(lat, orb_fwd, n_dim)
  ! Finding the closed orbit with backwards propagation is problematical since z does not behave 
  ! like it does with forward propagation.
  !! call closed_orbit_calc(lat, orb_back, n_dim, -1)
  df = orb_fwd(n)%vec - orb_fwd(0)%vec;  if (n_dim /= 6) df(5) = 0
  db = orb_fwd(n)%vec - orb_fwd(0)%vec;  if (n_dim /= 6) db(5) = 0

  write (1, '(a, i0, a, 6es14.6, 4x, es10.2)') '"Fwd', n_dim, '"  ABS 1E-10', orb_fwd(0)%vec, maxval(abs(df))
  !! write (1, '(a, i0, a, 6es14.6, 4x, es10.2)') '"Rev', n_dim, '"  ABS 1E-10', orb_back(0)%vec, maxval(abs(db))
enddo

!

call run1 ('ad05f_pos_RF.lat', 1)
call run1 ('ad05f_pos_RF.lat', 2)
call run1 ('ad05f_neg_RF.lat', 1)
call run1 ('ad05f_neg_RF.lat', 2)

!-------------------------------------------------------------------
contains

subroutine run1(lat_file, which)

type (lat_struct) lat
type (coord_struct), allocatable :: orb(:)

character(*) lat_file
integer which

!

call bmad_parser(lat_file, lat)
n = lat%n_ele_track

if (which == 2) then
  call closed_orbit_calc(lat, orb, 6)
  call lat_make_mat6(lat, -1, orb)
  call twiss_at_start(lat)
  write (1, '(3a, i0, a, t35, a, t50, f12.8)') '"', trim(lat_file), '-', which,  '-A Z"', 'ABS 5E-8', orb(0)%vec(5)
!  print '(3a, i0, a, t35, a, t50, 6f12.8)', '"', trim(lat_file), '-', which,  '-A Vec" ', 'ABS 5E-8', orb(0)%vec
!  print '(3a, i0, a, t35, a, t50, 6f12.8)', '"', trim(lat_file), '-', which,  '-A dVec"', 'ABS 5E-8', orb(n)%vec - orb(0)%vec
  deallocate(orb)
endif

call closed_orbit_calc(lat, orb, 6)
write (1, '(3a, i0, a, t35, a, t50, f12.8)') '"', trim(lat_file), '-', which,  '-B Z"', 'ABS 5E-8', orb(0)%vec(5)
!print '(3a, i0, a, t35, a, t50, 6f12.8)', '"', trim(lat_file), '-', which,  '-B Vec" ', 'ABS 5E-8', orb(0)%vec
!print '(3a, i0, a, t35, a, t50, 6f12.8)', '"', trim(lat_file), '-', which,  '-B dVec"', 'ABS 5E-8', orb(n)%vec - orb(0)%vec
call lat_make_mat6(lat, -1, orb)
call twiss_at_start(lat)
call twiss_propagate_all(lat)
call calc_z_tune(lat%branch(0))
write (1, '(3a, i0, a, t35, a, t50, f14.10)') '"', trim(lat_file), '-', which,  '-B Z_tune"', 'ABS 5E-8', lat%z%tune / twopi

end subroutine

end program
