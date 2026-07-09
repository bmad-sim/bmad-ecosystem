program gg_smoke
use bmad
implicit none
type (lat_struct) lat
type (ele_struct), pointer :: ele
type (em_field_struct) field
type (coord_struct) orb
real(rp) s

call bmad_parser ('gg_smoke.bmad', lat)
ele => lat%ele(1)

call init_coord (orb, lat%particle_start, ele, upstream_end$)
orb%vec(1) = 0.011_rp
orb%vec(3) = -0.007_rp

s = 0.0_rp   ! Base plane iz0 (z_rel = 0).
call em_field_calc (ele, lat%param, s, orb, .true., field, calc_potential = .true.)

print '(a, 3es22.14)', 'B = ', field%B
print '(a, 3es22.14)', 'A = ', field%A
end program
