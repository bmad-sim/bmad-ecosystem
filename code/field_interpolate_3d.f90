!+
! Function field_interpolate_3d (position, field_mesh, deltas, position0)
!
! Function to interpolate a 3d field.
! The interpolation is such that the derivative is continuous.
!
! Note: For "interpolation" outside of the region covered by the field_mesh
! it is assumed that the field is constant, Equal to the field at the
! boundary.
!
! Modules needed:
!   use bmad_interface
!
! Input:
!   position(3)       -- Real: (x, y, z) position.
!   field_mesh(:,:,:) -- Real: Grid of field points.
!   deltas(3)         -- Real: (dx, dy, dz) distances between mesh points.
!   position0(3)      -- Real, optional:  position at (ix0, iy0, iz0) where 
!                            (ix0, iy0, iz0) is the lower bound of the 
!                            filed_mesh(i, j, k) array. If not present then
!                            position0 is taken to be (0.0, 0.0, 0.0)
! Output:
!   field_interpolate_3d -- Real: interpolated field.
!-

real function field_interpolate_3d (position, field_mesh, deltas, position0)

  implicit none

  real, optional, intent(in) :: position0(3)
  real, intent(in) :: position(3), field_mesh(0:,0:,0:), deltas(3)

  real r(3), f(-1:2), g(-1:2), h(-1:2), r_frac(3)

  integer i0(3), ix, iy, iz, iix, iiy, iiz

!

  if (present(position0)) then
    r = (position - position0) / deltas 
  else
    r = position / deltas 
  endif

  i0 = int(r) 
  r_frac = r - i0

  do ix = -1, 2
   iix = min(max(ix + i0(1), 0), ubound(field_mesh, 1))
   do iy = -1, 2
      iiy = min(max(iy + i0(2), 0), ubound(field_mesh, 2))
      do iz = -1, 2
        iiz = min(max(iz + i0(3), 0), ubound(field_mesh, 3))
        f(iz) = field_mesh(iix, iiy, iiz)
      enddo
      g(iy) = interpolate_1d (r_frac(3), f)
    enddo
    h(ix) = interpolate_1d (r_frac(2), g)
  enddo
  field_interpolate_3d = interpolate_1d (r_frac(1), h)

!---------------------------------------------------------------

contains

! interpolation in 1 dimension using 4 equally spaced points: P1, P2, P3, P4.
!   x = interpolation point. 
!           x = 0 -> point is at P2.
!           x = 1 -> point is at P3.
! Interpolation is done so that the derivative is continuous.
! The interpolation uses a cubic polynomial

real function interpolate_1d (x, field_1d)

  implicit none

  real x, field_1d(4), df_2, df_3
  real c0, c1, c2, c3

!

  df_2 = (field_1d(3) - field_1d(1)) / 2   ! derivative at P2
  df_3 = (field_1d(4) - field_1d(2)) / 2   ! derivative at P3

  c0 = field_1d(2)
  c1 = df_2
  c2 = 3 * field_1d(3) - df_3 - 3 * field_1d(2) - 2 * df_2 
  c3 = df_3 - 2 * field_1d(3) + 2 * field_1d(2) + df_2

  interpolate_1d = c0 + c1 * x + c2 * x**2 + c3 * x**3

end function

end function
