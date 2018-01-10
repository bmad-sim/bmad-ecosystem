program wall3d_test

use wall3d_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (wall3d_struct), pointer :: wall

real(rp) r0(6), r1(6), perp(3), d_radius, r_wall, r_particle
real(rp) d1, d2, theta, max_d, origin(3), offset

integer i, ix

!

call bmad_parser ('wall3d_test.bmad', lat)

!

branch => lat%branch(1)
ix = element_at_s (lat, 3.0_rp, .true., 1)
ele => branch%ele(ix)

max_d = 0
offset = 1e-6

do i = 0, 19
  r0 = 0
  theta = i * twopi / 20.0
  r0(1:5:2) = [10*cos(theta), 10*sin(theta), 3.0_rp]

  d_radius = wall3d_d_radius (r0, ele, 1, perp, origin = origin)
  r_particle = norm2(r0(1:5:2) - origin)
  r_wall = r_particle - d_radius
  r0(1:5:2) = origin + (r0(1:5:2) - origin) * r_wall / r_particle

  d_radius = wall3d_d_radius (r0, ele, 1, perp)

  r1 = r0
  r1(1:5:2) = r1(1:5:2) + offset * [-perp(3), 0.0_rp, perp(1)]
  d1 = wall3d_d_radius (r1, ele) / offset

  r1 = r0
  r1(1:5:2) = r1(1:5:2) + offset * [0.0_rp, -perp(3), perp(2)]
  d2 = wall3d_d_radius (r1, ele) / offset

  print '(a, i4, 3es12.4, 3x, 3f10.5)', 'd_radius:', i, d_radius, d1, d2, perp
  max_d = max(max_d, abs(d_radius), abs(d1), abs(d2))
enddo

print *, 'Max:', max_d

!

wall => lat%branch(0)%wall3d(1)

open (1, file = 'output.now')


write (1, '(a, i0)') '"n_wall"  ABS 0  ', size(wall%section)

do i = 1, size(wall%section)
  write (1, '(a, i3, a, i3, f10.3)') '"Section', i, '" ABS 0  ', &
                                wall%section(i)%ix_ele, wall%section(i)%s
enddo

write (1, '(a, i0)') '"N_section" ABS 0  ', size(lat%branch(2)%wall3d(1)%section) 

end program
