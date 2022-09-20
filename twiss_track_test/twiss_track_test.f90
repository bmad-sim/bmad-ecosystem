program twiss_track_test

use bmad

implicit none

type (lat_struct), target :: lat, lat2
type (ele_struct) ele, ele0, ele1
type (coord_struct), allocatable :: orb(:), orb2(:)
type (coord_struct) orb0, orb1
type (normal_modes_struct) mode, mode2
type (rad_int_all_ele_struct), target :: rad_int, rad_int2, rad_int3

real(rp) chrom_x, chrom_y, delta_e
real(rp) m6(6,6), dorb(6), dt, rf_freq

integer i, j, n, n_lines, version, ix_cache, n_track

character(40) lattice
character(200) lat_file
character(200), allocatable :: lines(:)

!---------------------------------------------------------

open (2, file = 'output.now', recl = 200)

bmad_com%auto_bookkeeper = .false.

lat_file = "bmad_L9A18A000-_MOVEREC.lat"
call bmad_parser (lat_file, lat2)
call write_digested_bmad_file ('digested.file', lat2)
call read_digested_bmad_file ('digested.file', lat2, version)
lat = lat2

allocate (orb(0:lat%n_ele_max))     

bmad_com%radiation_damping_on = .true.

bmad_com%absolute_time_tracking = .true.
call closed_orbit_calc (lat, orb, 6)
write (2, '(a, 6es16.6)') '"Closed Orb 6T Start"  ABS 1e-12', orb(0)%vec
dorb = orb(lat%n_ele_track)%vec - orb(0)%vec
dt = orb(lat%n_ele_track)%t - orb(0)%t
rf_freq = lat%ele(76)%value(rf_frequency$)
dorb(5) = -orb(0)%beta * c_light * (dt - nint(dt * rf_freq) / rf_freq)
write (2, '(a, 6es12.4)') '"Closed Orb 6T Del"    ABS 2e-12', dorb

bmad_com%absolute_time_tracking = .false.
call closed_orbit_calc (lat, orb, 6)
write (2, '(a, 6es16.6)') '"Closed Orb 6 Start"  ABS 1e-12', orb(0)%vec
write (2, '(a, 6es12.4)') '"Closed Orb 6 Del"    ABS 1e-12', orb(lat%n_ele_track)%vec - orb(0)%vec

call set_on_off (rfcavity$, lat, off$)
call closed_orbit_calc (lat, orb, 5)
write (2, '(a, 6es16.6)') '"Closed Orb 5 Start"  ABS 1e-12', orb(0)%vec
write (2, '(a, 6es12.4)') '"Closed Orb 5 Del"    ABS 1e-12', orb(lat%n_ele_track)%vec - orb(0)%vec

call closed_orbit_calc (lat, orb, 4)
write (2, '(a, 6es16.6)') '"Closed Orb 4 Start"  ABS 1e-12', orb(0)%vec
write (2, '(a, 6es12.4)') '"Closed Orb 4 Del"    ABS 1e-12', orb(lat%n_ele_track)%vec - orb(0)%vec
call set_on_off (rfcavity$, lat, on$)


orb(0)%vec = 0
call twiss_at_start (lat)
call closed_orbit_calc (lat, orb, 6)
call lat_make_mat6 (lat, -1, orb)
call twiss_at_start (lat)
call twiss_propagate_all (lat)

!

allocate (lat%ele(96)%descrip)
lat%ele(96)%descrip = 'First'

ele = lat%ele(96)
ele%descrip = 'Second'

write (2, '(3a)') '"First Descrip"      STR  "', trim(lat%ele(96)%descrip), '"'
write (2, '(3a)') '"Second Descrip"     STR  "', trim(ele%descrip), '"'

delta_e = 0.0
call chrom_calc (lat, delta_e, chrom_x, chrom_y)

!

call mat_make_unit (m6)
open (1, file = 'twiss.out')
do n = 1, lat%n_ele_track
  write (1, *) '!------------------------------------'
  write (1, *) 'Index:', n
  call type_ele (lat%ele(n), .false., 0, .false., 0, .false., lines = lines, n_lines = n_lines)
  do i = 1, n_lines
    write (1, '(a)') lines(i)
  enddo
  m6 = matmul (lat%ele(n)%mat6, m6)
  do i = 1, 6
    write (1, '(6f11.5)') (m6(i, j), j = 1, 6)
  enddo
enddo

ix_cache = 0
call radiation_integrals (lat, orb, mode2, ix_cache, 0, rad_int2)
call ri_out(lat, 'rad_int_no_wig_cache.dat', rad_int2)

call radiation_integrals (lat, orb, mode, rad_int_by_ele = rad_int)
call ri_out(lat, 'rad_int_no_wig_no_cache.dat', rad_int)

call ri_diff('no_wig')

!

ele0 = lat%ele(0)
ele1 = lat%ele(lat%n_ele_track)
orb0 = orb(0)
orb1 = orb(lat%n_ele_track)

call data_out (orb1%vec(1) - orb0%vec(1), 1.0d-7, 'Dif:Orb(1)')
call data_out (orb1%vec(2) - orb0%vec(2), 1.0d-7, 'Dif:Orb(2)')
call data_out (orb1%vec(3) - orb0%vec(3), 1.0d-7, 'Dif:Orb(3)')
call data_out (orb1%vec(4) - orb0%vec(4), 1.0d-7, 'Dif:Orb(4)')
call data_out (orb1%vec(5) - orb0%vec(5), 1.0d-7, 'Dif:Orb(5)')
call data_out (orb1%vec(6) - orb0%vec(6), 1.0d-7, 'Dif:Orb(6)')
call data_out (ele1%a%beta - ele0%a%beta, 1.0d-7, 'Dif:Beta_X')
call data_out (ele1%b%beta - ele0%b%beta, 1.0d-7, 'Dif:Beta_Y')
call data_out (ele1%a%alpha - ele0%a%alpha, 1.0d-7, 'Dif:Alpha_X')
call data_out (ele1%b%alpha - ele0%b%alpha, 1.0d-7, 'Dif:Alpha_Y')
call data_out (ele1%a%eta - ele0%a%eta, 1.0d-7, 'Dif:Eta_X ')
call data_out (ele1%b%eta - ele0%b%eta, 1.0d-7, 'Dif:Eta_Y ')
call data_out (ele1%a%etap - ele0%a%etap, 1.0d-6, 'Dif:Etap_X')
call data_out (ele1%b%etap - ele0%b%etap, 1.0d-6, 'Dif:Etap_Y')

!----------------------------------------------------
! Error check.

!! print *
!! print *, 'Non-wiggler lattice check...'

ele = lat%ele(96)

call data_out (ele%a%beta,           1.0D-06, 'Lat1:Beta_a')
call data_out (ele%a%alpha,          1.0D-06, 'Lat1:Alpha_a')
call data_out (ele%a%eta,            1.0D-06, 'Lat1:Eta_a')
call data_out (ele%a%etap,           1.0D-06, 'Lat1:Etap_a')
call data_out (ele%x%eta,            1.0D-06, 'Lat1:Eta_x')
call data_out (ele%x%etap,           1.0D-06, 'Lat1:Etap_x')
call data_out (ele%a%phi,            1.0D-06, 'Lat1:Phi_a')
call data_out (ele%b%beta,           1.0D-06, 'Lat1:Beta_b')
call data_out (ele%b%alpha,          1.0D-06, 'Lat1:Alpha_y')
call data_out (ele%b%eta,            1.0D-06, 'Lat1:Eta_b')
call data_out (ele%b%etap,           1.0D-06, 'Lat1:Etap_y')
call data_out (ele%y%eta,            1.0D-06, 'Lat1:Eta_y')
call data_out (ele%y%etap,           1.0D-06, 'Lat1:Etap_y')
call data_out (ele%b%phi,            1.0D-06, 'Lat1:Phi_y')
call data_out (orb(96)%vec(1),       1.0D-10, 'Lat1:Orb X')
call data_out (orb(96)%vec(2),       1.0D-10, 'Lat1:Orb P_X')
call data_out (orb(96)%vec(3),       1.0D-10, 'Lat1:Orb Y')
call data_out (orb(96)%vec(4),       1.0D-10, 'Lat1:Orb P_Y')
call data_out (orb(96)%vec(5),       1.0D-10, 'Lat1:Orb Z')
call data_out (orb(96)%vec(6),       1.0D-10, 'Lat1:Orb P_Z')
call data_out (chrom_x,              1.0D-05, 'Lat1:Chrom_x')
call data_out (chrom_y,              1.0D-05, 'Lat1:Chrom_y')
call data_out (mode%synch_int(1),    1.0D-06, 'Lat1:Synch_int(1)')
call data_out (mode%synch_int(2),    1.0D-06, 'Lat1:Synch_int(2)')
call data_out (mode%synch_int(3),    1.0D-06, 'Lat1:Synch_int(3)')
call data_out (mode%sige_e,          1.0D-10, 'Lat1:Sige_e')
call data_out (mode%sig_z,           1.0D-08, 'Lat1:Sig_z')
call data_out (mode%e_loss,          1.0D-01, 'Lat1:E_loss')
call data_out (mode%a%emittance,     1.0D-12, 'Lat1:A%Emittance')
call data_out (mode%b%emittance,     1.0D-14, 'Lat1:B%Emittance')
call data_out (mode%z%emittance,     1.0D-11, 'Lat1:Z%Emittance')
call data_out (mode%a%synch_int(4),  1.0D-07, 'Lat1:A%Synch_int(4)')
call data_out (mode%a%synch_int(5),  1.0D-07, 'Lat1:A%Synch_int(5)')
call data_out (mode%b%synch_int(4),  1.0D-07, 'Lat1:B%Synch_int(4)')
call data_out (mode%b%synch_int(5),  1.0D-11, 'Lat1:B%Synch_int(5)')
call data_out (mode%z%synch_int(4),  1.0D-08, 'Lat1:Z%Synch_int(4)')
 
call set_z_tune (lat%branch(0), -0.05 * twopi)

!--------------------------------

write (2, *)

call bmad_parser ('bmad_12wig_20050626.lat', lat)
call set_on_off (rfcavity$, lat, off$)

orb(0)%vec = 0
call twiss_at_start (lat)
call closed_orbit_calc (lat, orb, 4)
call track_all (lat, orb)
call lat_make_mat6 (lat, -1, orb)
call twiss_at_start (lat)
call twiss_propagate_all (lat)

ix_cache = 0
call set_on_off (rfcavity$, lat, on$)
call radiation_integrals (lat, orb, mode2, ix_cache, 0, rad_int2)
call ri_out(lat, 'rad_int_wig_cache.dat', rad_int2)

call radiation_integrals (lat, orb, mode, rad_int_by_ele = rad_int)
call ri_out(lat, 'rad_int_wig_no_cache.dat', rad_int)

call ri_diff('wig')

call chrom_calc (lat, delta_e, chrom_x, chrom_y)

ele = lat%ele(96)
 
call data_out (ele%a%beta,           1.0D-05, 'Lat2:Beta_a')
call data_out (ele%a%alpha,          1.0D-05, 'Lat2:Alpha_a')
call data_out (ele%a%eta,            1.0D-05, 'Lat2:Eta_a')
call data_out (ele%a%etap,           1.0D-05, 'Lat2:Etap_a')
call data_out (ele%x%eta,            1.0D-05, 'Lat2:Eta_x')
call data_out (ele%x%etap,           1.0D-05, 'Lat2:Etap_x')
call data_out (ele%a%phi,            1.0D-05, 'Lat2:Phi_a')
call data_out (ele%b%beta,           1.0D-05, 'Lat2:Beta_b')
call data_out (ele%b%alpha,          1.0D-05, 'Lat2:Alpha_b')
call data_out (ele%b%eta,            1.0D-05, 'Lat2:Eta_b')
call data_out (ele%b%etap,           1.0D-05, 'Lat2:Etap_b')
call data_out (ele%y%eta,            1.0D-05, 'Lat2:Eta_y')
call data_out (ele%y%etap,           1.0D-05, 'Lat2:Etap_y')
call data_out (ele%b%phi,            1.0D-05, 'Lat2:Phi_b')
call data_out (orb(96)%vec(1),       1.0D-07, 'Lat2:Orb X')
call data_out (orb(96)%vec(2),       1.0D-07, 'Lat2:Orb P_X')
call data_out (orb(96)%vec(3),       1.0D-07, 'Lat2:Orb Y')
call data_out (orb(96)%vec(4),       1.0D-07, 'Lat2:Orb P_Y')
call data_out (orb(96)%vec(5),       1.0D-07, 'Lat2:Orb Z')
call data_out (orb(96)%vec(6),       1.0D-07, 'Lat2:Orb P_Z')
call data_out (chrom_x,              1.0D-04, 'Lat2:Chrom_x')
call data_out (chrom_y,              1.0D-04, 'Lat2:Chrom_y')
call data_out (mode%synch_int(1),    1.0D-06, 'Lat2:Synch_int(1)')
call data_out (mode%synch_int(2),    1.0D-06, 'Lat2:Synch_int(2)')
call data_out (mode%synch_int(3),    1.0D-06, 'Lat2:Synch_int(3)')
call data_out (mode%sige_e,          1.0D-10, 'Lat2:Sige_e')
call data_out (mode%sig_z,           1.0D-08, 'Lat2:Sig_z')
call data_out (mode%e_loss,          1.0D-01, 'Lat2:E_loss')
call data_out (mode%a%emittance,     1.0D-12, 'Lat2:A%Emittance')
call data_out (mode%b%emittance,     1.0D-14, 'Lat2:B%Emittance')
call data_out (mode%z%emittance,     1.0D-11, 'Lat2:Z%Emittance')
call data_out (mode%a%synch_int(4),  1.0D-07, 'Lat2:A%Synch_int(4)')
call data_out (mode%a%synch_int(5),  1.0D-07, 'Lat2:A%Synch_int(5)')
call data_out (mode%b%synch_int(4),  1.0D-07, 'Lat2:B%Synch_int(4)')
call data_out (mode%b%synch_int(5),  1.0D-11, 'Lat2:B%Synch_int(5)')
call data_out (mode%z%synch_int(4),  1.0D-08, 'Lat2:Z%Synch_int(4)')

write (2, '(a, l1, a)') '"Lat2:Lat"      STR  "', associated(lat2%ele(100)%branch, lat2%branch(0)), '"'

!-----------------------------------------------------------------------------------
contains

subroutine data_out (now, err_tol, what)

implicit none

real(rp) now, theory, err_tol
integer ix
character(*) what
character(200) line

!

write (2, '(3a, t30, a, es10.1, es25.14)') '"', what, '" ', 'ABS', err_tol, now

end subroutine

!----------------------------------------------------------------------------------
! contains

subroutine ri_out (lat, file_name, rad_int)

type (lat_struct) lat
type (rad_int_all_ele_struct), target :: rad_int
type (rad_int1_struct), pointer :: r1
character(*) file_name
integer i

!

open (1, file = file_name)

write (1, '(45x, a)') 'I1          I2          I3          I4a         I4b         I5a         I5b'
do i = 1, lat%n_ele_max
  r1 => rad_int%branch(0)%ele(i)
  if (all([r1%i1, r1%i2, r1%i3, r1%i4a, r1%i4b, r1%i5a, r1%i5b] == 0)) cycle
  write (1, '(i4, 2x, a30, 7es12.3)') i, lat%ele(i)%name, r1%i1, r1%i2, r1%i3, r1%i4a, r1%i4b, r1%i5a, r1%i5b
enddo
write (1, '(45x, a)') 'I1          I2          I3          I4a         I4b         I5a         I5b'

close (1)

end subroutine

!------------------------------------------------------------------------------
! contains

subroutine ri_diff(str)

character(*) str

call ri_diff1('Cache Diff: I1-' // str, rad_int%branch(0)%ele%i1,  rad_int2%branch(0)%ele%i1)
call ri_diff1('Cache Diff: I2-' // str, rad_int%branch(0)%ele%i2,  rad_int2%branch(0)%ele%i2)
call ri_diff1('Cache Diff: I3-' // str, rad_int%branch(0)%ele%i3,  rad_int2%branch(0)%ele%i3)
call ri_diff1('Cache Diff: I4a-' // str, rad_int%branch(0)%ele%i4a, rad_int2%branch(0)%ele%i4a)
call ri_diff1('Cache Diff: I4b-' // str, rad_int%branch(0)%ele%i4b, rad_int2%branch(0)%ele%i4b)
call ri_diff1('Cache Diff: I5a-' // str, rad_int%branch(0)%ele%i5a, rad_int2%branch(0)%ele%i5a)
call ri_diff1('Cache Diff: I5b-' // str, rad_int%branch(0)%ele%i5b, rad_int2%branch(0)%ele%i5b)

end subroutine

!---------------------------------------------------------------------------
! contains

subroutine ri_diff1 (str, vec1, vec2)

character(*) str
real(rp) vec1(0:), vec2(0:), mdiff, mrdiff
real(rp) dvec, vmax
integer i, im, imr

!

vmax = maxval(abs(vec1))
mdiff = -1
mrdiff = -1

do i = 0, ubound(vec1, 1)

  dvec = abs(vec1(i) - vec2(i))
  if (dvec > mdiff) then
    mdiff = dvec
    im = i
  endif

  if (vec1(i) == 0 .and. vec2(i) == 0) cycle
  if (abs(vec1(i)) < 1e-10 * vmax .and. abs(vec1(i)) < 1e-10 * vmax) cycle

  dvec = 2 * dvec / (abs(vec1(i)) + abs(vec2(i)))
  if (dvec > mrdiff) then
    mrdiff = dvec
    imr = i
  endif
enddo

write (2, '(3a, t30, a, i6, 3es20.9)') '"', str, '"', 'ABS   1e-8', im, mdiff/vmax!, vec1(im), vec2(im)

end subroutine

end program


