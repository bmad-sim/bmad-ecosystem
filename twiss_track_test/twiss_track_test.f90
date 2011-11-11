program twiss_track_test

use bmad
use z_tune_mod

implicit none

type (lat_struct)  lat, lat2
type (ele_struct) ele, ele0, ele1
type (coord_struct), allocatable :: orb(:)   
type (coord_struct) orb0, orb1
type (normal_modes_struct) mode
type (rad_int_common_struct) rad_int, rad_int2

real(rp) chrom_x, chrom_y, delta_e
real(rp) m6(6,6)

integer i, j, n, n_lines, version, ix_cache

character(40) lattice
character(200) lat_file
character(100), pointer :: lines(:)

!

bmad_com%auto_bookkeeper = .false.

lat_file = "bmad_L9A18A000-_MOVEREC.lat"
call bmad_parser (lat_file, lat2)
call write_digested_bmad_file ('digested.file', lat2, 0)
call read_digested_bmad_file ('digested.file', lat2, version)
lat = lat2

print *, 'Number of elements:', lat%n_ele_max
print *, 'Precision:', rp

allocate (orb(0:lat%n_ele_max))     

bmad_com%rel_tolerance = 1e-7
bmad_com%abs_tolerance = 1e-10

orb(0)%vec = 0
!! call set_on_off (rfcavity$, lat, off$)  
call twiss_at_start (lat)
call closed_orbit_calc (lat, orb, 6)
call lat_make_mat6 (lat, -1, orb)
call twiss_at_start (lat)
call twiss_propagate_all (lat)

call type_ele (lat%ele(96), .false., 6, .true., radians$, .true., lat)
print *
call type_coord (orb(96))

!

allocate (lat%ele(96)%descrip)
lat%ele(96)%descrip = 'First'

ele = lat%ele(96)
ele%descrip = 'Second'

print *, 'Descriptions Should be: First, Second.'
print *, '    ', trim(lat%ele(96)%descrip)
print *, '    ', trim(ele%descrip)

delta_e = 0.0
call chrom_calc (lat, delta_e, chrom_x, chrom_y)
print *
print *, 'Chrom:', chrom_x, chrom_y

!

call mat_make_unit (m6)
open (1, file = 'twiss.out')
do n = 1, lat%n_ele_track
  write (1, *) '!------------------------------------'
  write (1, *) 'Index:', n
  call type2_ele (lat%ele(n), lines, n_lines, .false., 0, .false., 0, .false., lat)  
  do i = 1, n_lines
    write (1, '(a)') lines(i)
  enddo
  m6 = matmul (lat%ele(n)%mat6, m6)
  do i = 1, 6
    write (1, '(6f11.5)') (m6(i, j), j = 1, 6)
  enddo
enddo

ix_cache = 0
call radiation_integrals (lat, orb, mode, ix_cache, rad_int2)
call ri_out('rad_int_no_wig_cache.dat', rad_int2)

call radiation_integrals (lat, orb, mode, rad_int_by_ele = rad_int)
call ri_out('rad_int_no_wig_no_cache.dat', rad_int)

call ri_diff

open (1, file = 'rad_int.out')
do n = 1, lat%n_ele_track
  if (rad_int%i4a(n) == 0) cycle
  write (1, '(i5, 2x, a16, 2f12.8)') n, lat%ele(n)%name, &
                                    rad_int%i4a(n), sum(rad_int%i4a(1:n))
enddo
close (1)

print *
print '(a,1p,6e18.10)', 'Synch_int1/3:', mode%synch_int
print '(a,1p,6e18.10)', 'Sig_e/z:     ', mode%sige_e, mode%sig_z, mode%e_loss
print '(a,1p,6e18.10)', 'a', mode%a%emittance, mode%a%synch_int(4:5)
print '(a,1p,6e18.10)', 'b', mode%b%emittance, mode%b%synch_int(4:5)
print '(a,1p,6e18.10)', 'z', mode%z%emittance, mode%z%synch_int(4)

!

ele0 = lat%ele(0)
ele1 = lat%ele(lat%n_ele_track)
orb0 = orb(0)
orb1 = orb(lat%n_ele_track)
print *
print *, 'End Check:'
print '(a, 2f13.8)', 'Orb(1): ', orb1%vec(1), orb1%vec(1) - orb0%vec(1)
print '(a, 2f13.8)', 'Orb(2): ', orb1%vec(2), orb1%vec(2) - orb0%vec(2)
print '(a, 2f13.8)', 'Orb(3): ', orb1%vec(3), orb1%vec(3) - orb0%vec(3)
print '(a, 2f13.8)', 'Orb(4): ', orb1%vec(4), orb1%vec(4) - orb0%vec(4)
print '(a, 2f13.8)', 'Orb(5): ', orb1%vec(5), orb1%vec(5) - orb0%vec(5)
print '(a, 2f13.8)', 'Orb(6): ', orb1%vec(6), orb1%vec(6) - orb0%vec(6)
print '(a, 2f13.8)', 'Beta_X: ', ele1%a%beta, ele1%a%beta - ele0%a%beta
print '(a, 2f13.8)', 'Beta_Y: ', ele1%b%beta, ele1%b%beta - ele0%b%beta
print '(a, 2f13.8)', 'Alpha_X:', ele1%a%alpha, ele1%a%alpha - ele0%a%alpha
print '(a, 2f13.8)', 'Alpha_Y:', ele1%b%alpha, ele1%b%alpha - ele0%b%alpha
print '(a, 2f13.8)', 'Eta_X:  ', ele1%a%eta, ele1%a%eta - ele0%a%eta
print '(a, 2f13.8)', 'Eta_Y:  ', ele1%b%eta, ele1%b%eta - ele0%b%eta
print '(a, 2f13.8)', 'Etap_X: ', ele1%a%etap, ele1%a%etap - ele0%a%etap
print '(a, 2f13.8)', 'Etap_Y: ', ele1%b%etap, ele1%b%etap - ele0%b%etap

!----------------------------------------------------
! Error check.

print *
print *, 'Non-wiggler lattice check...'

ele = lat%ele(96)

open (2, file = 'output.now', recl = 200)

call check (ele%a%beta,           1.0D-06, 'Lat1:Beta_a')
call check (ele%a%alpha,          1.0D-06, 'Lat1:Alpha_a')
call check (ele%a%eta,            1.0D-06, 'Lat1:Eta_a')
call check (ele%a%etap,           1.0D-06, 'Lat1:Etap_a')
call check (ele%x%eta,            1.0D-06, 'Lat1:Eta_x')
call check (ele%x%etap,           1.0D-06, 'Lat1:Etap_x')
call check (ele%a%phi,            1.0D-06, 'Lat1:Phi_a')
call check (ele%b%beta,           1.0D-06, 'Lat1:Beta_b')
call check (ele%b%alpha,          1.0D-06, 'Lat1:Alpha_y')
call check (ele%b%eta,            1.0D-06, 'Lat1:Eta_b')
call check (ele%b%etap,           1.0D-06, 'Lat1:Etap_y')
call check (ele%y%eta,            1.0D-06, 'Lat1:Eta_y')
call check (ele%y%etap,           1.0D-06, 'Lat1:Etap_y')
call check (ele%b%phi,            1.0D-06, 'Lat1:Phi_y')
call check (orb(96)%vec(1),       1.0D-10, 'Lat1:Orb X')
call check (orb(96)%vec(2),       1.0D-10, 'Lat1:Orb P_X')
call check (orb(96)%vec(3),       1.0D-10, 'Lat1:Orb Y')
call check (orb(96)%vec(4),       1.0D-10, 'Lat1:Orb P_Y')
call check (orb(96)%vec(5),       1.0D-10, 'Lat1:Orb Z')
call check (orb(96)%vec(6),       1.0D-10, 'Lat1:Orb P_Z')
call check (chrom_x,              1.0D-05, 'Lat1:Chrom_x')
call check (chrom_y,              1.0D-05, 'Lat1:Chrom_y')
call check (mode%synch_int(1),    1.0D-06, 'Lat1:Synch_int(1)')
call check (mode%synch_int(2),    1.0D-06, 'Lat1:Synch_int(2)')
call check (mode%synch_int(3),    1.0D-06, 'Lat1:Synch_int(3)')
call check (mode%sige_e,          1.0D-10, 'Lat1:Sige_e')
call check (mode%sig_z,           1.0D-08, 'Lat1:Sig_z')
call check (mode%e_loss,          1.0D-01, 'Lat1:E_loss')
call check (mode%a%emittance,     1.0D-12, 'Lat1:A%Emittance')
call check (mode%b%emittance,     1.0D-14, 'Lat1:B%Emittance')
call check (mode%z%emittance,     1.0D-11, 'Lat1:Z%Emittance')
call check (mode%a%synch_int(4),  1.0D-07, 'Lat1:A%Synch_int(4)')
call check (mode%a%synch_int(5),  1.0D-07, 'Lat1:A%Synch_int(5)')
call check (mode%b%synch_int(4),  1.0D-07, 'Lat1:B%Synch_int(4)')
call check (mode%b%synch_int(5),  1.0D-11, 'Lat1:B%Synch_int(5)')
call check (mode%z%synch_int(4),  1.0D-08, 'Lat1:Z%Synch_int(4)')
 
call set_z_tune (lat, -0.05 * twopi)

!--------------------------------

write (2, *)
write (2, *)

print *
print *, 'Wiggler lattice check...'

call bmad_parser ('bmad_12wig_20050626.lat', lat)

orb(0)%vec = 0
call twiss_at_start (lat)
call closed_orbit_calc (lat, orb, 4)
call track_all (lat, orb)
call lat_make_mat6 (lat, -1, orb)
call twiss_at_start (lat)
call twiss_propagate_all (lat)

ix_cache = 0
call radiation_integrals (lat, orb, mode, ix_cache, rad_int2)
call ri_out('rad_int_wig_cache.dat', rad_int2)

call radiation_integrals (lat, orb, mode, rad_int_by_ele = rad_int)
call ri_out('rad_int_wig_no_cache.dat', rad_int)

call ri_diff

call chrom_calc (lat, delta_e, chrom_x, chrom_y)

ele = lat%ele(96)
 
call check (ele%a%beta,           1.0D-05, 'Lat2:Beta_a')
call check (ele%a%alpha,          1.0D-05, 'Lat2:Alpha_a')
call check (ele%a%eta,            1.0D-05, 'Lat2:Eta_a')
call check (ele%a%etap,           1.0D-05, 'Lat2:Etap_a')
call check (ele%x%eta,            1.0D-05, 'Lat2:Eta_x')
call check (ele%x%etap,           1.0D-05, 'Lat2:Etap_x')
call check (ele%a%phi,            1.0D-05, 'Lat2:Phi_a')
call check (ele%b%beta,           1.0D-05, 'Lat2:Beta_b')
call check (ele%b%alpha,          1.0D-05, 'Lat2:Alpha_b')
call check (ele%b%eta,            1.0D-05, 'Lat2:Eta_b')
call check (ele%b%etap,           1.0D-05, 'Lat2:Etap_b')
call check (ele%y%eta,            1.0D-05, 'Lat2:Eta_y')
call check (ele%y%etap,           1.0D-05, 'Lat2:Etap_y')
call check (ele%b%phi,            1.0D-05, 'Lat2:Phi_b')
call check (orb(96)%vec(1),       1.0D-07, 'Lat2:Orb X')
call check (orb(96)%vec(2),       1.0D-07, 'Lat2:Orb P_X')
call check (orb(96)%vec(3),       1.0D-07, 'Lat2:Orb Y')
call check (orb(96)%vec(4),       1.0D-07, 'Lat2:Orb P_Y')
call check (orb(96)%vec(5),       1.0D-07, 'Lat2:Orb Z')
call check (orb(96)%vec(6),       1.0D-07, 'Lat2:Orb P_Z')
call check (chrom_x,              1.0D-04, 'Lat2:Chrom_x')
call check (chrom_y,              1.0D-04, 'Lat2:Chrom_y')
call check (mode%synch_int(1),    1.0D-06, 'Lat2:Synch_int(1)')
call check (mode%synch_int(2),    1.0D-06, 'Lat2:Synch_int(2)')
call check (mode%synch_int(3),    1.0D-06, 'Lat2:Synch_int(3)')
call check (mode%sige_e,          1.0D-10, 'Lat2:Sige_e')
call check (mode%sig_z,           1.0D-08, 'Lat2:Sig_z')
call check (mode%e_loss,          1.0D-01, 'Lat2:E_loss')
call check (mode%a%emittance,     1.0D-12, 'Lat2:A%Emittance')
call check (mode%b%emittance,     1.0D-14, 'Lat2:B%Emittance')
call check (mode%z%emittance,     1.0D-11, 'Lat2:Z%Emittance')
call check (mode%a%synch_int(4),  1.0D-07, 'Lat2:A%Synch_int(4)')
call check (mode%a%synch_int(5),  1.0D-07, 'Lat2:A%Synch_int(5)')
call check (mode%b%synch_int(4),  1.0D-07, 'Lat2:B%Synch_int(4)')
call check (mode%b%synch_int(5),  1.0D-11, 'Lat2:B%Synch_int(5)')
call check (mode%z%synch_int(4),  1.0D-08, 'Lat2:Z%Synch_int(4)')

!--------------------------------------------------------------------
contains

subroutine check (now, err_tol, what)

implicit none

real(rp) now, theory, err_tol
integer ix
character(*) what
character(200) line

!

write (2, '(3a, t25, es9.1, a, es22.14)') '"', what, '" : REL :', err_tol, ' : ', now

end subroutine

!--------------------------------------------------------------------
! contains

subroutine ri_out (file_name, rad_int)

type (rad_int_common_struct) rad_int
character(*) file_name
integer i

!

open (1, file = file_name)

write (1, '(a)') '             I1          I2          I3          I4a         I4b         I5a         I5b'
do i = 1, lat%n_ele_track
  if (all([rad_int%i1(i), rad_int%i2(i), rad_int%i3(i), &
                        rad_int%i4a(i), rad_int%i4b(i), rad_int%i5a(i), rad_int%i5b(i)] == 0)) cycle
  write (1, '(i4, 7es12.3)') i, rad_int%i1(i), rad_int%i2(i), rad_int%i3(i), &
                        rad_int%i4a(i), rad_int%i4b(i), rad_int%i5a(i), rad_int%i5b(i)
enddo
write (1, '(a)') '             I1          I2          I3          I4a         I4b         I5a         I5b'

close (1)

end subroutine

!--------------------------------------------------------------------
! contains

subroutine ri_diff

print *, 'Max radiation integral diffs between caching and no caching:'
call ri_diff1('I1 ', rad_int%i1,  rad_int2%i1)
call ri_diff1('I2 ', rad_int%i2,  rad_int2%i2)
call ri_diff1('I3 ', rad_int%i3,  rad_int2%i3)
call ri_diff1('I4a', rad_int%i4a, rad_int2%i4a)
call ri_diff1('I4b', rad_int%i4b, rad_int2%i4b)
call ri_diff1('I5a', rad_int%i5a, rad_int2%i5a)
call ri_diff1('I5b', rad_int%i5b, rad_int2%i5b)

end subroutine

!--------------------------------------------------------------------
! contains

subroutine ri_diff1 (str, vec1, vec2)

character(*) str
real(rp) vec1(0:), vec2(0:), mdiff, mrdiff
real(rp) dvec, vmax
integer i, im, imr

!

vmax = maxval(abs(vec1))
mdiff = 0
mrdiff = 0

do i = 0, ubound(vec1, 1)
  if (vec1(i) == 0 .and. vec2(i) == 0) cycle
  if (abs(vec1(i)) < 1e-10*vmax .and. abs(vec2(i)) < 1e-10*vmax) cycle

  dvec = abs(vec1(i) - vec2(i))
  if (dvec > mdiff) then
    mdiff = dvec
    im = i
  endif

  dvec = 2 * dvec / (abs(vec1(i)) + abs(vec2(i)))
  if (dvec > mrdiff) then
    mrdiff = dvec
    imr = i
  endif
enddo

print '(a, i6, es12.3, i6, f8.4)', str, im, mdiff/vmax, imr, mrdiff

end subroutine

end program


