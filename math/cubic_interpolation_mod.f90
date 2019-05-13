!+
! Module cubic_interpolation_mod
!
! Module implementing cubic interpolation of a field that has been measured on a grid.
!
! See the Wikipedia article on bicubic interpolation on a 2D grid.
!
! The tricubic interpolation on a 3D grid is based upon:
!   "Tricubic Interpolation in Three Dimensions", 
!   F. Lekien, J. Marsden, Journal of Numerical Methods and Engineering, (2005).
!
! Note: Cubic interpolation is C^1. That is, it is continuous and has a continuous first derivative.
!-

module cubic_interpolation_mod

use sim_utils_struct

implicit none

!--------------------------------------------
! Real bicubic interpolation structures

! Field and derivatives at a grid point needed for bicubic interpolation.
! The derivatives are normalized by the distance between grid points dx, dy.
! For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy

type field1_at_2D_pt_struct
  real(rp) f                                 ! Field
  real(rp) df_dx, df_dy                      ! Normalized field 1st derivatives
  real(rp) d2f_dxdy                          ! Normalized field 2nd derivative
end type

! Field and derivatives at the 4 grid points defining the box for bicubic interpolation.

type field_at_2D_box_struct
  type (field1_at_2D_pt_struct) pt(0:1, 0:1)
end type

! Coefficients needed to evaluate bicubic interpolation.

type bicubic_coef_struct
  ! (ix, iy) = index of grid element coefs are for.
  integer :: ix = int_garbage$, iy = int_garbage$
  real(rp) :: coef(0:3,0:3) = 0 ! Coefs
end type

!--------------------------------------------
! Complex bicubic interpolation structures

! Field and derivatives at a grid point needed for bicubic interpolation.
! The derivatives are normalized by the distance between grid points dx, dy.
! For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy

type cmplx_field1_at_2D_pt_struct
  complex(rp) f                                 ! Field
  complex(rp) df_dx, df_dy                      ! Normalized field 1st derivatives
  complex(rp) d2f_dxdy                          ! Normalized field 2nd derivative
end type

! Field and derivatives at the 4 grid points defining the box for bicubic interpolation.

type cmplx_field_at_2D_box_struct
  type (cmplx_field1_at_2D_pt_struct) pt(0:1, 0:1)
end type

! Coefficients needed to evaluate bicubic interpolation.

type bicubic_cmplx_coef_struct
  ! (ix, iy) = index of grid element coefs are for.
  integer :: ix = int_garbage$, iy = int_garbage$
  complex(rp) :: coef(0:3,0:3) = 0 ! Coefs
end type

!--------------------------------------------
! Real tricubic interpolation structures

! Field and derivatives at a grid point needed for tricubic interpolation.
! The derivatives are normalized by the distance between grid points dx, dy, dz.
! For example: %d2f_dxdz (structure component) = d^2f/dxdz * dx * dz

type field1_at_3D_pt_struct
  real(rp) f                                 ! Field
  real(rp) df_dx, df_dy, df_dz               ! Normalized field 1st derivatives
  real(rp) d2f_dxdy, d2f_dxdz, d2f_dydz      ! Normalized field 2nd derivatives
  real(rp) d3f_dxdydz                        ! Normalized field 3rd derivative
end type

! Field and derivatives at the 8 grid points defining the box for tricubic interpolation.

type field_at_3D_box_struct
  type (field1_at_3D_pt_struct) pt(0:1, 0:1, 0:1)
end type

! Coefficients needed to evaluate tricubic interpolation.

type tricubic_coef_struct
  ! (ix, iy, iz) = index of grid element coefs are for.
  integer :: ix = int_garbage$, iy = int_garbage$, iz = int_garbage$
  real(rp) :: coef(0:3,0:3,0:3) = 0 ! Coefs
end type

!--------------------------------------------
! Complex tricubic interpolation structures

! Field and derivatives at a grid point needed for tricubic interpolation.
! The derivatives are normalized by the distance between grid points dx, dy, dz.
! For example: %d2f_dxdz (structure component) = d^2f/dxdz * dx * dz

type cmplx_field1_at_3D_pt_struct
  complex(rp) f                                 ! Field
  complex(rp) df_dx, df_dy, df_dz               ! Normalized field 1st derivatives
  complex(rp) d2f_dxdy, d2f_dxdz, d2f_dydz      ! Normalized field 2nd derivatives
  complex(rp) d3f_dxdydz                        ! Normalized field 3rd derivative
end type

! Field and derivatives at the 8 grid points defining the box for tricubic interpolation.

type cmplx_field_at_3D_box_struct
  type (cmplx_field1_at_3D_pt_struct) pt(0:1, 0:1, 0:1)
end type

! Coefficients needed to evaluate tricubic interpolation.

type tricubic_cmplx_coef_struct
  ! (ix, iy, iz) = index of grid element coefs are for.
  integer :: ix = int_garbage$, iy = int_garbage$, iz = int_garbage$
  complex(rp) :: coef(0:3,0:3,0:3) = 0 ! Coefs
end type

!---------------------------------------------------------

type m16_row_struct
  integer pt(16)
end type

type m16_matrix_struct
  type (m16_row_struct) row(16)
end type

private m16_row_struct, m16_matrix_struct

type (m16_matrix_struct), parameter, private :: m16_matrix = m16_matrix_struct( [ &
 m16_row_struct([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
 m16_row_struct([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
 m16_row_struct([-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
 m16_row_struct([2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
 m16_row_struct([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]), &
 m16_row_struct([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]), &
 m16_row_struct([0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0]), &
 m16_row_struct([0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0]), &
 m16_row_struct([-3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0]), &
 m16_row_struct([0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0]), &
 m16_row_struct([9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1]), &
 m16_row_struct([-6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1]), &
 m16_row_struct([2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0]), &
 m16_row_struct([0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0]), &
 m16_row_struct([-6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1]), &
 m16_row_struct([4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1]) ])

type m64_row_struct
  integer pt(64)
end type

type m64_matrix_struct
  type (m64_row_struct) row(64)
end type

private m64_row_struct, m64_matrix_struct

type (m64_matrix_struct), parameter, private :: m64_matrix = m64_matrix_struct( [ &
  m64_row_struct([ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 9,-9,-9, 9, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([-6, 6, 6,-6, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([-6, 6, 6,-6, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 4,-4,-4, 4, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0]), &
  m64_row_struct([-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 9,-9, 0, 0,-9, 9, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([-6, 6, 0, 0, 6,-6, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9, 0, 0,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0]), &
  m64_row_struct([ 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0]), &
  m64_row_struct([-27,27,27,-27,27,-27,-27,27,-18,-9,18, 9,18, 9,-18,-9,-18,18,-9, 9,18,-18, 9,-9,-18,18,18,-18,-9, 9, 9,-9,-12,-6,-6,-3,12, 6, 6, 3,-12,-6,12, 6,-6,-3, 6, 3,-12,12,-6, 6,-6, 6,-3, 3,-8,-4,-4,-2,-4,-2,-2,-1]), &
  m64_row_struct([18,-18,-18,18,-18,18,18,-18, 9, 9,-9,-9,-9,-9, 9, 9,12,-12, 6,-6,-12,12,-6, 6,12,-12,-12,12, 6,-6,-6, 6, 6, 6, 3, 3,-6,-6,-3,-3, 6, 6,-6,-6, 3, 3,-3,-3, 8,-8, 4,-4, 4,-4, 2,-2, 4, 4, 2, 2, 2, 2, 1, 1]), &
  m64_row_struct([-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0]), &
  m64_row_struct([18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1]), &
  m64_row_struct([-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1]), &
  m64_row_struct([ 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([-6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0]), &
  m64_row_struct([-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0]), &
  m64_row_struct([18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1]), &
  m64_row_struct([-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1]), &
  m64_row_struct([ 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]), &
  m64_row_struct([ 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0]), &
  m64_row_struct([-12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1]), &
  m64_row_struct([ 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1]) ])

contains

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine bicubic_compute_field_at_2D_box (field, ig0, ix, iy, extrapolation, field_at_box, err_flag)
!
! Routine to compute the field and normalized derivatives at 4 points forming a box on a 2D grid.
! The box spans from (ix, iy) to (ix+1, iy+1).
!
! Note: The grid size must be at least 2 x 2.
!
! The grid box is "just outside" the grid if some of the box points are within the grid region and some are outside.
! In this case, the extrapolation argument is used to determine how the field is calculated at outside box points.
!
! The output of this routine can be used in the routine bicubic_interpolation_coefs.
!
! Input:
!   field(:,:)      -- real(rp): Field at grid points.
!   ig0(2)          -- integer: grid lower bounds.
!   ix, iy          -- integer: Grid box lower bounds
!   extrapolation   -- character(*): Used if the grid box is just outside of the grid.
!                       = 'NONE'      ! No extrapolation and err_flag is set to True.
!                       = 'ZERO'      ! Assume field and derivatives are zero at the outside box points.
!                       = 'LINEAR'    ! Linear extrapolation to the outside box points.
!                       = 'CONSTANT'  ! Constant field is assumed outside the grid.
!
! Output:
!   field_at_box    -- field_at_2D_box_struct: Field and derivatives.
!   err_flag        -- logical: Set True if there is an error such as the box is not in or just outside the grid.
!-

subroutine bicubic_compute_field_at_2D_box (field, ig0, ix, iy, extrapolation, field_at_box, err_flag)

type (field_at_2D_box_struct) field_at_box

integer ig0(2), ig1(2)
real(rp) field(ig0(1):, ig0(2):)
real(rp) f(-1:2, -1:2)

integer i, j, ix, iy, ixx, iyy, ix0, ix1, iy0, iy1
logical err_flag

character(*) extrapolation

! Init

ig1 = ubound(field)
err_flag = .true.

if (any(ig0 >= ig1)) return

if (extrapolation == 'NONE') then
  if (ix < ig0(1) .or. ix >= ig1(1)) return
  if (iy < ig0(2) .or. iy >= ig1(2)) return
else
  if (ix < ig0(1) - 1 .or. ix > ig1(1)) return
  if (iy < ig0(2) - 1 .or. iy > ig1(2)) return
endif

err_flag = .false.

! Setup g mini-grid

ix0 = ig0(1); ix1 = ig1(1)
iy0 = ig0(2); iy1 = ig1(2)

do i = -1, 2
do j = -1, 2
  ixx = ix + i
  iyy = iy + j
  if (ix >= ix0 .and. ix <= ix1 .and. iy >= iy0 .and. iy <= iy1) then
    f(i,j) = field(ixx, iyy)
    cycle
  endif

  select case (extrapolation)
  case ('ZERO')
    f(i,j) = 0

  case ('LINEAR')
    if (ixx == ix0-1 .and. iyy == iy0-1) then
      f(i,j) = 3 * field(ixx+1, iyy+1) - field(ixx+2, iyy+1) - field(ixx+1, iyy+2)
    elseif (ixx == ix0-1 .and. iyy == iy1+1) then
      f(i,j) = 3 * field(ixx+1, iyy-1) - field(ixx+2, iyy-1) - field(ixx+1, iyy-2)
    elseif (ixx == ix1+1 .and. iyy == iy0-1) then
      f(i,j) = 3 * field(ixx-1, iyy+1) - field(ixx-2, iyy+1) - field(ixx-1, iyy+2)
    elseif (ixx == ix1+1 .and. iyy == iy1+1) then
      f(i,j) = 3 * field(ixx-1, iyy-1) - field(ixx-2, iyy-1) - field(ixx-1, iyy-2)
    elseif (ixx == ix0) then
      f(i,j) = 2 * field(ixx+1, iyy) - field(ixx+2, iyy)
    elseif (ixx == ix1) then
      f(i,j) = 2 * field(ixx-1, iyy) - field(ixx-2, iyy)
    elseif (iyy == iy0) then
      f(i,j) = 2 * field(ixx, iyy+1) - field(ixx, iyy+2)
    elseif (iyy == iy1) then
      f(i,j) = 2 * field(ixx, iyy-1) - field(ixx, iyy-2)
    endif

  case ('CONSTANT')
    if (ixx == ix0-1) ixx = ix0
    if (ixx == ix1+1) ixx = ix1
    if (iyy == iy0-1) iyy = iy0
    if (iyy == iy1+1) iyy = iy1
    f(i,j) = field(ixx, iyy)
  end select
enddo
enddo

! Evaluate field and derivatives.

do i = 0, 1
do j = 0, 1
  field_at_box%pt(i,j)%f = f(i, j)
  field_at_box%pt(i,j)%df_dx = (f(i+1, j) - f(i-1, j)) / 2
  field_at_box%pt(i,j)%df_dy = (f(i, j+1) - f(i, j-1)) / 2
  field_at_box%pt(i,j)%d2f_dxdy = (f(i+1, j+1) - f(i+1, j-1) - f(i-1, j+1) + f(i-1, j-1)) / 4
enddo
enddo

end subroutine bicubic_compute_field_at_2D_box

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! subroutine bicubic_interpolation_coefs (field_at_box, bi_coef)
! 
! Routine to compute the coefficients for bicubic interpolation.
!
! Use the routine bicubic_eval to evaluate the interpolation function.
!
! Note: The derivatives in field_at_box are normalized by the distance between grid points dx, dy.
! For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy
!
! Input:
!   field_at_box   -- field_at_2D_box_struct: Field and normalized derivatives at the 4 grid points.
!
! Output:
!   bi_coef         -- bicubic_coef_struct: Coefficients.
!-

subroutine bicubic_interpolation_coefs (field_at_box, bi_coef)

type (field_at_2D_box_struct), target :: field_at_box
type (bicubic_coef_struct) bi_coef
type (field1_at_2D_pt_struct), pointer :: fg

real(rp) field_array(16)
integer i, j, ij

!

do i = 0, 1
do j = 0, 1
  ij = i + 2*j + 1
  fg => field_at_box%pt(i,j)
  field_array(ij:ij+12:4) = [fg%f, fg%df_dx, fg%df_dy, fg%d2f_dxdy]
enddo
enddo

do i = 0, 3
do j = 0, 3
  ij = i + 4*j + 1
  bi_coef%coef(i,j) = sum(m16_matrix%row(ij)%pt(:) * field_array(:))
enddo
enddo

end subroutine bicubic_interpolation_coefs

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Function bicubic_eval (x_norm, y_norm, bi_coef, df_dx, df_dy) result (f_val)
!
! Routine to evaluate a bicubic interpolating function.
!
! Use the routine bicubic_interpolation_coefs to generate bi_coef.
!
! Note: In the equations below, the four points of the grid box being interpolated range 
! from (x0, y0) to (x0+dx, y0+dy).
!
! Input:
!   x_norm    -- real(rp): x_norm = (x - x0) / dx 
!   y_norm    -- real(rp): y_norm = (y - y0) / dy
!   bi_coef   -- bicubic_coef_struct: Coefficients.
!
! Output:
!   f_val     -- real(rp): Value of f.
!   df_dx     -- real(rp): Normalized first derivative: True df/dx = df_dx * dx 
!   df_dy     -- real(rp): Normalized first derivative: True df/dy = df_dy * dy
!- 

function bicubic_eval (x_norm, y_norm, bi_coef, df_dx, df_dy) result (f_val)

type (bicubic_coef_struct) bi_coef

real(rp) x_norm, y_norm, z_norm, f_val
real(rp), optional :: df_dx, df_dy
real(rp) x_power(0:3), y_power(0:3), z_power(0:3)

integer i, j

! Init

x_power(0) = 1
y_power(0) = 1

do i = 1, 3
  x_power(i) = x_power(i-1) * x_norm
  y_power(i) = y_power(i-1) * y_norm
enddo

! Evaluate f_val

f_val = 0
do i = 0, 3
do j = 0, 3
  f_val = f_val + x_power(i) * y_power(j) * bi_coef%coef(i,j)
enddo
enddo

! Evaluate derivatives

if (present(df_dx)) then
  df_dx = 0
  do i = 1, 3
  do j = 0, 3
    df_dx = df_dx + i * x_power(i-1) * y_power(j) * bi_coef%coef(i,j)
  enddo
  enddo
endif

if (present(df_dy)) then
  df_dy = 0
  do i = 0, 3
  do j = 1, 3
    df_dy = df_dy + j * x_power(i) * y_power(j-1) * bi_coef%coef(i,j)
  enddo
  enddo
endif

end function bicubic_eval

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine tricubic_compute_field_at_3D_box (field, ig0, ix, iy, iz, extrapolation, field_at_box, err_flag)
!
! Routine to compute the field and normalized derivatives at 8 points forming a box on a 2D grid.
! The box spans from (ix, iy, iz) to (ix+1, iy+1, iz+1).
!
! Note: The grid size must be at least 2 x 2 x 2.
!
! The grid box is "just outside" the grid if some of the box points are within the grid region and some are outside.
! In this case, the extrapolation argument is used to determine how the field is calculated at outside box points.
!
! The output of this routine can be used in the routine tricubic_interpolation_coefs.
!
! Input:
!   field(:,:,:)    -- real(rp): Field at grid points.
!   ig0(3)          -- integer: grid lower bounds.
!   ix, iy, iz      -- integer: Grid box lower bounds
!   extrapolation   -- character(*): Used if the grid box is just outside of the grid.
!                       = 'NONE'      ! No extrapolation and err_flag is set to True.
!                       = 'ZERO'      ! Assume field and derivatives are zero at the outside box points.
!                       = 'LINEAR'    ! Linear extrapolation to the outside box points.
!                       = 'CONSTANT'  ! Constant field is assumed outside the grid.
!
! Output:
!   field_at_box    -- field_at_3D_box_struct: Field and derivatives.
!   err_flag        -- logical: Set True if there is an error such as the box is not in or just outside the grid.
!-

subroutine tricubic_compute_field_at_3D_box (field, ig0, ix, iy, iz, extrapolation, field_at_box, err_flag)

type (field_at_3D_box_struct) field_at_box

integer ig0(3), ig1(3)
real(rp) field(ig0(1):, ig0(2):, ig0(3):)
real(rp) f(-1:2, -1:2, -1:2)

integer i, j, k, ix, iy, iz, ixx, iyy, izz, ix0, ix1, iy0, iy1, iz0, iz1
logical err_flag

character(*) extrapolation

! Init

ig1 = ubound(field)
err_flag = .true.

if (any(ig0 >= ig1)) return

if (extrapolation == 'NONE') then
  if (ix < ig0(1) .or. ix >= ig1(1)) return
  if (iy < ig0(2) .or. iy >= ig1(2)) return
  if (iz < ig0(2) .or. iz >= ig1(2)) return
else
  if (ix < ig0(1) - 1 .or. ix > ig1(1)) return
  if (iy < ig0(2) - 1 .or. iy > ig1(2)) return
  if (iz < ig0(2) - 1 .or. iz > ig1(2)) return
endif

err_flag = .false.

! Setup g mini-grid

ix0 = ig0(1); ix1 = ig1(1)
iy0 = ig0(2); iy1 = ig1(2)
iz0 = ig0(3); iz1 = ig1(3)

do i = -1, 2
do j = -1, 2
do k = -1, 2
  ixx = ix + i
  iyy = iy + j
  izz = iz + k
  if (ix >= ix0 .and. ix <= ix1 .and. iy >= iy0 .and. iy <= iy1 .and. iz >= iz0 .and. iz <= iz1) then
    f(i,j,k) = field(ixx, iyy, izz)
    cycle
  endif

  select case (extrapolation)
  case ('ZERO')
    f(i,j,k) = 0

  case ('LINEAR')
    if (ixx == ix0-1 .and. iyy == iy0-1 .and. izz == iz0-1) then
      f(i,j,k) = 4 * field(ixx+1, iyy+1, izz+1) - field(ixx+2, iyy+1, izz+1) - field(ixx+1, iyy+2, izz+1) - field(ixx+1, iyy+1, izz+2)
    elseif (ixx == ix0-1 .and. iyy == iy0-1 .and. izz == iz1+1) then
      f(i,j,k) = 4 * field(ixx+1, iyy+1, izz-1) - field(ixx+2, iyy+1, izz-1) - field(ixx+1, iyy+2, izz-1) - field(ixx+1, iyy+1, izz-2)
    elseif (ixx == ix0-1 .and. iyy == iy1+1 .and. izz == iz0-1) then
      f(i,j,k) = 4 * field(ixx+1, iyy-1, izz+1) - field(ixx+2, iyy-1, izz+1) - field(ixx+1, iyy-2, izz+1) - field(ixx+1, iyy-1, izz+2)
    elseif (ixx == ix0-1 .and. iyy == iy1+1 .and. izz == iz1+1) then
      f(i,j,k) = 4 * field(ixx+1, iyy-1, izz-1) - field(ixx+2, iyy-1, izz-1) - field(ixx+1, iyy-2, izz-1) - field(ixx+1, iyy-1, izz-2)
    elseif (ixx == ix1+1 .and. iyy == iy0-1 .and. izz == iz0-1) then
      f(i,j,k) = 4 * field(ixx-1, iyy+1, izz+1) - field(ixx-2, iyy+1, izz+1) - field(ixx-1, iyy+2, izz+1) - field(ixx-1, iyy+1, izz+2)
    elseif (ixx == ix1+1 .and. iyy == iy0-1 .and. izz == iz1+1) then
      f(i,j,k) = 4 * field(ixx-1, iyy+1, izz-1) - field(ixx-2, iyy+1, izz-1) - field(ixx-1, iyy+2, izz-1) - field(ixx-1, iyy+1, izz-2)
    elseif (ixx == ix1+1 .and. iyy == iy1+1 .and. izz == iz0-1) then
      f(i,j,k) = 4 * field(ixx-1, iyy-1, izz+1) - field(ixx-2, iyy-1, izz+1) - field(ixx-1, iyy-2, izz+1) - field(ixx-1, iyy-1, izz+2)
    elseif (ixx == ix1+1 .and. iyy == iy1+1 .and. izz == iz1+1) then
      f(i,j,k) = 4 * field(ixx-1, iyy-1, izz-1) - field(ixx-2, iyy-1, izz-1) - field(ixx-1, iyy-2, izz-1) - field(ixx-1, iyy-1, izz-2)
    elseif (ixx == ix0) then
      f(i,j,k) = 2 * field(ixx+1, iyy, izz) - field(ixx+2, iyy, izz)
    elseif (ixx == ix1) then
      f(i,j,k) = 2 * field(ixx-1, iyy, izz) - field(ixx-2, iyy, izz)
    elseif (iyy == iy0) then
      f(i,j,k) = 2 * field(ixx, iyy+1, izz) - field(ixx, iyy+2, izz)
    elseif (iyy == iy1) then
      f(i,j,k) = 2 * field(ixx, iyy-1, izz) - field(ixx, iyy-2, izz)
    elseif (izz == iz0) then
      f(i,j,k) = 2 * field(ixx, iyy, izz+1) - field(ixx, iyy, izz+2)
    elseif (izz == iz1) then
      f(i,j,k) = 2 * field(ixx, iyy, izz-1) - field(ixx, iyy, izz-2)
    endif

  case ('CONSTANT')
    if (ixx == ix0-1) ixx = ix0
    if (ixx == ix1+1) ixx = ix1
    if (iyy == iy0-1) iyy = iy0
    if (iyy == iy1+1) iyy = iy1
    if (izz == iz0-1) izz = iz0
    if (izz == iz1+1) izz = iz1
    f(i,j,k) = field(ixx, iyy, izz)
  end select
enddo
enddo
enddo

! Evaluate field and derivatives.

do i = 0, 1
do j = 0, 1
do k = 0, 1
  field_at_box%pt(i,j,k)%f = f(i, j, k)
  field_at_box%pt(i,j,k)%df_dx = (f(i+1, j, k) - f(i-1, j, k)) / 2
  field_at_box%pt(i,j,k)%df_dy = (f(i, j+1, k) - f(i, j-1, k)) / 2
  field_at_box%pt(i,j,k)%df_dz = (f(i, j, k+1) - f(i, j, k+1)) / 2
  field_at_box%pt(i,j,k)%d2f_dxdy = (f(i+1, j+1, k) - f(i+1, j-1, k) - f(i-1, j+1, k) + f(i-1, j-1, k)) / 4
  field_at_box%pt(i,j,k)%d2f_dxdz = (f(i+1, j, k+1) - f(i+1, j, k-1) - f(i-1, j, k+1) + f(i-1, j, k-1)) / 4
  field_at_box%pt(i,j,k)%d2f_dydz = (f(i, j+1, k+1) - f(i, j+1, k-1) - f(i, j-1, k+1) + f(i, j-1, k-1)) / 4
  field_at_box%pt(i,j,k)%d3f_dxdydz = (f(i+1, j+1, k+1) - f(i-1, j+1, k+1) - f(i+1, j-1, k+1) - f(i+1, j+1, k-1) + &
                                       f(i-1, j-1, k+1) + f(i-1, j+1, k-1) + f(i+1, j-1, k-1) - f(i-1, j-1, k-1)) / 8
enddo
enddo
enddo

end subroutine tricubic_compute_field_at_3D_box

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine tricubic_interpolation_coefs (field_at_box, tri_coef)
! 
! Routine to compute the coefficients for tricubic interpolation.
!
! Use the routine tricubic_eval to evaluate the interpolation function.
!
! Note: The derivatives in field_at_box are normalized by the distance between grid points dx, dy.
! For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy
!
! Input:
!   field_at_box   -- field_at_2D_box_struct: Field and normalized derivatives at the 4 grid points.
!
! Output:
!   tri_coef        -- tricubic_coef_struct: Coefficients.
!-

subroutine tricubic_interpolation_coefs (field_at_box, tri_coef)

type (field_at_3D_box_struct), target :: field_at_box
type (tricubic_coef_struct) tri_coef
type (field1_at_3D_pt_struct), pointer :: fg

real(rp) field_array(64)
integer i, j, k, ijk

!

do i = 0, 1
do j = 0, 1
do k = 0, 1
  ijk = i + 2*j + 4*k + 1
  fg => field_at_box%pt(i,j,k)
  field_array(ijk:ijk+56:8) = [fg%f, fg%df_dx, fg%df_dy, fg%df_dz, fg%d2f_dxdy, fg%d2f_dxdz, fg%d2f_dydz, fg%d3f_dxdydz]
enddo
enddo
enddo

do i = 0, 3
do j = 0, 3
do k = 0, 3
  ijk = i + 4*j + 16*k + 1
  tri_coef%coef(i,j,k) = sum(m64_matrix%row(ijk)%pt(:) * field_array(:))
enddo
enddo
enddo

end subroutine tricubic_interpolation_coefs

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Function tricubic_eval (x_norm, y_norm, z_norm, tri_coef, df_dx, df_dy, df_dz) result (f_val)
!
! Routine to evaluate a tricubic interpolating function.
!
! Use the routine tricubic_interpolation_coefs to generate tri_coef.
!
! Note: In the equations below, the eight points of the grid box being interpolated range 
! from (x0, y0, z0) to (x0+dx, y0+dy, z0+dz).
!
! Input:
!   x_norm    -- real(rp): x_norm = (x - x0) / dx 
!   y_norm    -- real(rp): y_norm = (y - y0) / dy
!   z_norm    -- real(rp): z_norm = (z - z0) / dz
!   tri_coef  -- tricubic_coef_struct: Coefficients.
!
! Output:
!   f_val     -- real(rp): Value of f.
!   df_dx     -- real(rp): Normalized first derivative: True df/dx = df_dx * dx
!   df_dy     -- real(rp): Normalized first derivative: True df/dy = df_dy * dy
!   df_dz     -- real(rp): Normalized first derivative: True df/dz = df_dz * dz
!- 

function tricubic_eval (x_norm, y_norm, z_norm, tri_coef, df_dx, df_dy, df_dz) result (f_val)

type (tricubic_coef_struct) tri_coef

real(rp) x_norm, y_norm, z_norm, f_val
real(rp), optional :: df_dx, df_dy, df_dz
real(rp) x_power(0:3), y_power(0:3), z_power(0:3)

integer i, j, k

! Init

x_power(0) = 1
y_power(0) = 1
z_power(0) = 1

do i = 1, 3
  x_power(i) = x_power(i-1) * x_norm
  y_power(i) = y_power(i-1) * y_norm
  z_power(i) = z_power(i-1) * z_norm
enddo

! Evaluate f_val

f_val = 0
do i = 0, 3
do j = 0, 3
do k = 0, 3
  f_val = f_val + x_power(i) * y_power(j) * z_power(k) * tri_coef%coef(i,j,k)
enddo
enddo
enddo

! Evaluate derivatives

if (present(df_dx)) then
  df_dx = 0
  do i = 1, 3
  do j = 0, 3
  do k = 0, 3
    df_dx = df_dx + i * x_power(i-1) * y_power(j) * z_power(k) * tri_coef%coef(i,j,k)
  enddo
  enddo
  enddo
endif

if (present(df_dy)) then
  df_dy = 0
  do i = 0, 3
  do j = 1, 3
  do k = 0, 3
    df_dy = df_dy + j * x_power(i) * y_power(j-1) * z_power(k) * tri_coef%coef(i,j,k)
  enddo
  enddo
  enddo
endif

if (present(df_dz)) then
  df_dz = 0
  do i = 0, 3
  do j = 0, 3
  do k = 1, 3
    df_dz = df_dz + k * x_power(i) * y_power(j) * z_power(k-1) * tri_coef%coef(i,j,k)
  enddo
  enddo
  enddo
endif

end function tricubic_eval

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine bicubic_compute_cmplx_field_at_2D_box (field, ig0, ix, iy, extrapolation, field_at_box, err_flag)
!
! Routine to compute the field and normalized derivatives at 4 points forming a box on a 2D grid.
! The box spans from (ix, iy) to (ix+1, iy+1).
!
! Note: The grid size must be at least 2 x 2.
!
! The grid box is "just outside" the grid if some of the box points are within the grid region and some are outside.
! In this case, the extrapolation argument is used to determine how the field is calculated at outside box points.
!
! The output of this routine can be used in the routine bicubic_interpolation_coefs.
!
! Input:
!   field(:,:)      -- complex(rp): Field at grid points.
!   ig0(2)          -- integer: grid lower bounds.
!   ix, iy          -- integer: Grid box lower bounds
!   extrapolation   -- character(*): Used if the grid box is just outside of the grid.
!                       = 'NONE'      ! No extrapolation and err_flag is set to True.
!                       = 'ZERO'      ! Assume field and derivatives are zero at the outside box points.
!                       = 'LINEAR'    ! Linear extrapolation to the outside box points.
!                       = 'CONSTANT'  ! Constant field is assumed outside the grid.
!
! Output:
!   field_at_box    -- cmplx_field_at_2D_box_struct: Field and derivatives.
!   err_flag        -- logical: Set True if there is an error such as the box is not in or just outside the grid.
!-

subroutine bicubic_compute_cmplx_field_at_2D_box (field, ig0, ix, iy, extrapolation, field_at_box, err_flag)

type (cmplx_field_at_2D_box_struct) field_at_box

integer ig0(2), ig1(2)
complex(rp) field(ig0(1):, ig0(2):)
complex(rp) f(-1:2, -1:2)

integer i, j, ix, iy, ixx, iyy, ix0, ix1, iy0, iy1
logical err_flag

character(*) extrapolation

! Init

ig1 = ubound(field)
err_flag = .true.

if (any(ig0 >= ig1)) return

if (extrapolation == 'NONE') then
  if (ix < ig0(1) .or. ix >= ig1(1)) return
  if (iy < ig0(2) .or. iy >= ig1(2)) return
else
  if (ix < ig0(1) - 1 .or. ix > ig1(1)) return
  if (iy < ig0(2) - 1 .or. iy > ig1(2)) return
endif

err_flag = .false.

! Setup g mini-grid

ix0 = ig0(1); ix1 = ig1(1)
iy0 = ig0(2); iy1 = ig1(2)

do i = -1, 2
do j = -1, 2
  ixx = ix + i
  iyy = iy + j
  if (ix >= ix0 .and. ix <= ix1 .and. iy >= iy0 .and. iy <= iy1) then
    f(i,j) = field(ixx, iyy)
    cycle
  endif

  select case (extrapolation)
  case ('ZERO')
    f(i,j) = 0

  case ('LINEAR')
    if (ixx == ix0-1 .and. iyy == iy0-1) then
      f(i,j) = 3 * field(ixx+1, iyy+1) - field(ixx+2, iyy+1) - field(ixx+1, iyy+2)
    elseif (ixx == ix0-1 .and. iyy == iy1+1) then
      f(i,j) = 3 * field(ixx+1, iyy-1) - field(ixx+2, iyy-1) - field(ixx+1, iyy-2)
    elseif (ixx == ix1+1 .and. iyy == iy0-1) then
      f(i,j) = 3 * field(ixx-1, iyy+1) - field(ixx-2, iyy+1) - field(ixx-1, iyy+2)
    elseif (ixx == ix1+1 .and. iyy == iy1+1) then
      f(i,j) = 3 * field(ixx-1, iyy-1) - field(ixx-2, iyy-1) - field(ixx-1, iyy-2)
    elseif (ixx == ix0) then
      f(i,j) = 2 * field(ixx+1, iyy) - field(ixx+2, iyy)
    elseif (ixx == ix1) then
      f(i,j) = 2 * field(ixx-1, iyy) - field(ixx-2, iyy)
    elseif (iyy == iy0) then
      f(i,j) = 2 * field(ixx, iyy+1) - field(ixx, iyy+2)
    elseif (iyy == iy1) then
      f(i,j) = 2 * field(ixx, iyy-1) - field(ixx, iyy-2)
    endif

  case ('CONSTANT')
    if (ixx == ix0-1) ixx = ix0
    if (ixx == ix1+1) ixx = ix1
    if (iyy == iy0-1) iyy = iy0
    if (iyy == iy1+1) iyy = iy1
    f(i,j) = field(ixx, iyy)
  end select
enddo
enddo

! Evaluate field and derivatives.

do i = 0, 1
do j = 0, 1
  field_at_box%pt(i,j)%f = f(i, j)
  field_at_box%pt(i,j)%df_dx = (f(i+1, j) - f(i-1, j)) / 2
  field_at_box%pt(i,j)%df_dy = (f(i, j+1) - f(i, j-1)) / 2
  field_at_box%pt(i,j)%d2f_dxdy = (f(i+1, j+1) - f(i+1, j-1) - f(i-1, j+1) + f(i-1, j-1)) / 4
enddo
enddo

end subroutine bicubic_compute_cmplx_field_at_2D_box

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine bicubic_interpolation_cmplx_coefs (field_at_box, bi_coef)
! 
! Routine to compute the complex coefficients for bicubic interpolation.
!
! Use the routine bicubic_cmplx_eval to evaluate the interpolation function.
!
! Note: The derivatives in field_at_box are normalized by the distance between grid points dx, dy.
! For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy
!
! Input:
!   field_at_box   -- field_at_2D_box_struct: Field and normalized derivatives at the 4 grid points.
!
! Output:
!   bi_coef         -- bicubic_cmplx_coef_struct: Coefficients.
!-

subroutine bicubic_interpolation_cmplx_coefs (field_at_box, bi_coef)

type (cmplx_field_at_2D_box_struct), target :: field_at_box
type (bicubic_cmplx_coef_struct) bi_coef
type (cmplx_field1_at_2D_pt_struct), pointer :: fg

complex(rp) field_array(16)
integer i, j, ij

!

do i = 0, 1
do j = 0, 1
  ij = i + 2*j + 1
  fg => field_at_box%pt(i,j)
  field_array(ij:ij+12:4) = [fg%f, fg%df_dx, fg%df_dy, fg%d2f_dxdy]
enddo
enddo

do i = 0, 3
do j = 0, 3
  ij = i + 4*j + 1
  bi_coef%coef(i,j) = sum(m16_matrix%row(ij)%pt(:) * field_array(:))
enddo
enddo

end subroutine bicubic_interpolation_cmplx_coefs

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Function bicubic_cmplx_eval (x_norm, y_norm, bi_coef, df_dx, df_dy) result (f_val)
!
! Routine to evaluate a bicubic interpolating complex function.
!
! Use the routine bicubic_interpolation_cmplx_coefs to generate bi_coef.
!
! Note: In the equations below, the four points of the grid box being interpolated range 
! from (x0, y0) to (x0+dx, y0+dy).
!
! Input:
!   x_norm    -- real(rp): x_norm = (x - x0) / dx 
!   y_norm    -- real(rp): y_norm = (y - y0) / dy
!   bi_coef   -- bicubic_cmplx_coef_struct: Coefficients.
!
! Output:
!   f_val     -- complex(rp): Value of f.
!   df_dx     -- complex(rp): Normalized first derivative: True df/dx = df_dx * dx 
!   df_dy     -- complex(rp): Normalized first derivative: True df/dy = df_dy * dy
!- 

function bicubic_cmplx_eval (x_norm, y_norm, bi_coef, df_dx, df_dy) result (f_val)

type (bicubic_cmplx_coef_struct) bi_coef

real(rp) x_norm, y_norm, z_norm
real(rp) x_power(0:3), y_power(0:3), z_power(0:3)

complex(rp) f_val
complex(rp), optional :: df_dx, df_dy

integer i, j

! Init

x_power(0) = 1
y_power(0) = 1

do i = 1, 3
  x_power(i) = x_power(i-1) * x_norm
  y_power(i) = y_power(i-1) * y_norm
enddo

! Evaluate f_val

f_val = 0
do i = 0, 3
do j = 0, 3
  f_val = f_val + x_power(i) * y_power(j) * bi_coef%coef(i,j)
enddo
enddo

! Evaluate derivatives

if (present(df_dx)) then
  df_dx = 0
  do i = 1, 3
  do j = 0, 3
    df_dx = df_dx + i * x_power(i-1) * y_power(j) * bi_coef%coef(i,j)
  enddo
  enddo
endif

if (present(df_dy)) then
  df_dy = 0
  do i = 0, 3
  do j = 1, 3
    df_dy = df_dy + j * x_power(i) * y_power(j-1) * bi_coef%coef(i,j)
  enddo
  enddo
endif

end function bicubic_cmplx_eval

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine tricubic_compute_cmplx_field_at_3D_box (field, ig0, ix, iy, iz, extrapolation, field_at_box, err_flag)
!
! Routine to compute the field and normalized derivatives at 8 points forming a box on a 2D grid.
! The box spans from (ix, iy, iz) to (ix+1, iy+1, iz+1).
!
! Note: The grid size must be at least 2 x 2 x 2.
!
! The grid box is "just outside" the grid if some of the box points are within the grid region and some are outside.
! In this case, the extrapolation argument is used to determine how the field is calculated at outside box points.
!
! The output of this routine can be used in the routine tricubic_interpolation_coefs.
!
! Input:
!   field(:,:,:)    -- complex(rp): Field at grid points.
!   ig0(3)          -- integer: grid lower bounds.
!   ix, iy, iz      -- integer: Grid box lower bounds
!   extrapolation   -- character(*): Used if the grid box is just outside of the grid.
!                       = 'NONE'      ! No extrapolation and err_flag is set to True.
!                       = 'ZERO'      ! Assume field and derivatives are zero at the outside box points.
!                       = 'LINEAR'    ! Linear extrapolation to the outside box points.
!                       = 'CONSTANT'  ! Constant field is assumed outside the grid.
!
! Output:
!   field_at_box    -- cmplx_field_at_3D_box_struct: Field and derivatives.
!   err_flag        -- logical: Set True if there is an error such as the box is not in or just outside the grid.
!-

subroutine tricubic_compute_cmplx_field_at_3D_box (field, ig0, ix, iy, iz, extrapolation, field_at_box, err_flag)

type (cmplx_field_at_3D_box_struct) field_at_box

integer ig0(3), ig1(3)
complex(rp) field(ig0(1):, ig0(2):, ig0(3):)
complex(rp) f(-1:2, -1:2, -1:2)

integer i, j, k, ix, iy, iz, ixx, iyy, izz, ix0, ix1, iy0, iy1, iz0, iz1
logical err_flag

character(*) extrapolation

! Init

ig1 = ubound(field)
err_flag = .true.

if (any(ig0 >= ig1)) return

if (extrapolation == 'NONE') then
  if (ix < ig0(1) .or. ix >= ig1(1)) return
  if (iy < ig0(2) .or. iy >= ig1(2)) return
  if (iz < ig0(2) .or. iz >= ig1(2)) return
else
  if (ix < ig0(1) - 1 .or. ix > ig1(1)) return
  if (iy < ig0(2) - 1 .or. iy > ig1(2)) return
  if (iz < ig0(2) - 1 .or. iz > ig1(2)) return
endif

err_flag = .false.

! Setup g mini-grid

ix0 = ig0(1); ix1 = ig1(1)
iy0 = ig0(2); iy1 = ig1(2)
iz0 = ig0(3); iz1 = ig1(3)

do i = -1, 2
do j = -1, 2
do k = -1, 2
  ixx = ix + i
  iyy = iy + j
  izz = iz + k
  if (ix >= ix0 .and. ix <= ix1 .and. iy >= iy0 .and. iy <= iy1 .and. iz >= iz0 .and. iz <= iz1) then
    f(i,j,k) = field(ixx, iyy, izz)
    cycle
  endif

  select case (extrapolation)
  case ('ZERO')
    f(i,j,k) = 0

  case ('LINEAR')
    if (ixx == ix0-1 .and. iyy == iy0-1 .and. izz == iz0-1) then
      f(i,j,k) = 4 * field(ixx+1, iyy+1, izz+1) - field(ixx+2, iyy+1, izz+1) - field(ixx+1, iyy+2, izz+1) - field(ixx+1, iyy+1, izz+2)
    elseif (ixx == ix0-1 .and. iyy == iy0-1 .and. izz == iz1+1) then
      f(i,j,k) = 4 * field(ixx+1, iyy+1, izz-1) - field(ixx+2, iyy+1, izz-1) - field(ixx+1, iyy+2, izz-1) - field(ixx+1, iyy+1, izz-2)
    elseif (ixx == ix0-1 .and. iyy == iy1+1 .and. izz == iz0-1) then
      f(i,j,k) = 4 * field(ixx+1, iyy-1, izz+1) - field(ixx+2, iyy-1, izz+1) - field(ixx+1, iyy-2, izz+1) - field(ixx+1, iyy-1, izz+2)
    elseif (ixx == ix0-1 .and. iyy == iy1+1 .and. izz == iz1+1) then
      f(i,j,k) = 4 * field(ixx+1, iyy-1, izz-1) - field(ixx+2, iyy-1, izz-1) - field(ixx+1, iyy-2, izz-1) - field(ixx+1, iyy-1, izz-2)
    elseif (ixx == ix1+1 .and. iyy == iy0-1 .and. izz == iz0-1) then
      f(i,j,k) = 4 * field(ixx-1, iyy+1, izz+1) - field(ixx-2, iyy+1, izz+1) - field(ixx-1, iyy+2, izz+1) - field(ixx-1, iyy+1, izz+2)
    elseif (ixx == ix1+1 .and. iyy == iy0-1 .and. izz == iz1+1) then
      f(i,j,k) = 4 * field(ixx-1, iyy+1, izz-1) - field(ixx-2, iyy+1, izz-1) - field(ixx-1, iyy+2, izz-1) - field(ixx-1, iyy+1, izz-2)
    elseif (ixx == ix1+1 .and. iyy == iy1+1 .and. izz == iz0-1) then
      f(i,j,k) = 4 * field(ixx-1, iyy-1, izz+1) - field(ixx-2, iyy-1, izz+1) - field(ixx-1, iyy-2, izz+1) - field(ixx-1, iyy-1, izz+2)
    elseif (ixx == ix1+1 .and. iyy == iy1+1 .and. izz == iz1+1) then
      f(i,j,k) = 4 * field(ixx-1, iyy-1, izz-1) - field(ixx-2, iyy-1, izz-1) - field(ixx-1, iyy-2, izz-1) - field(ixx-1, iyy-1, izz-2)
    elseif (ixx == ix0) then
      f(i,j,k) = 2 * field(ixx+1, iyy, izz) - field(ixx+2, iyy, izz)
    elseif (ixx == ix1) then
      f(i,j,k) = 2 * field(ixx-1, iyy, izz) - field(ixx-2, iyy, izz)
    elseif (iyy == iy0) then
      f(i,j,k) = 2 * field(ixx, iyy+1, izz) - field(ixx, iyy+2, izz)
    elseif (iyy == iy1) then
      f(i,j,k) = 2 * field(ixx, iyy-1, izz) - field(ixx, iyy-2, izz)
    elseif (izz == iz0) then
      f(i,j,k) = 2 * field(ixx, iyy, izz+1) - field(ixx, iyy, izz+2)
    elseif (izz == iz1) then
      f(i,j,k) = 2 * field(ixx, iyy, izz-1) - field(ixx, iyy, izz-2)
    endif

  case ('CONSTANT')
    if (ixx == ix0-1) ixx = ix0
    if (ixx == ix1+1) ixx = ix1
    if (iyy == iy0-1) iyy = iy0
    if (iyy == iy1+1) iyy = iy1
    if (izz == iz0-1) izz = iz0
    if (izz == iz1+1) izz = iz1
    f(i,j,k) = field(ixx, iyy, izz)
  end select
enddo
enddo
enddo

! Evaluate field and derivatives.

do i = 0, 1
do j = 0, 1
do k = 0, 1
  field_at_box%pt(i,j,k)%f = f(i, j, k)
  field_at_box%pt(i,j,k)%df_dx = (f(i+1, j, k) - f(i-1, j, k)) / 2
  field_at_box%pt(i,j,k)%df_dy = (f(i, j+1, k) - f(i, j-1, k)) / 2
  field_at_box%pt(i,j,k)%df_dz = (f(i, j, k+1) - f(i, j, k+1)) / 2
  field_at_box%pt(i,j,k)%d2f_dxdy = (f(i+1, j+1, k) - f(i+1, j-1, k) - f(i-1, j+1, k) + f(i-1, j-1, k)) / 4
  field_at_box%pt(i,j,k)%d2f_dxdz = (f(i+1, j, k+1) - f(i+1, j, k-1) - f(i-1, j, k+1) + f(i-1, j, k-1)) / 4
  field_at_box%pt(i,j,k)%d2f_dydz = (f(i, j+1, k+1) - f(i, j+1, k-1) - f(i, j-1, k+1) + f(i, j-1, k-1)) / 4
  field_at_box%pt(i,j,k)%d3f_dxdydz = (f(i+1, j+1, k+1) - f(i-1, j+1, k+1) - f(i+1, j-1, k+1) - f(i+1, j+1, k-1) + &
                                       f(i-1, j-1, k+1) + f(i-1, j+1, k-1) + f(i+1, j-1, k-1) - f(i-1, j-1, k-1)) / 8
enddo
enddo
enddo

end subroutine tricubic_compute_cmplx_field_at_3D_box

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine tricubic_interpolation_cmplx_coefs (field_at_box, tri_coef)
! 
! Routine to compute the complex coefficients for tricubic interpolation.
!
! Use the routine tricubic_cmplx_eval to evaluate the interpolation function.
!
! Note: The derivatives in field_at_box are normalized by the distance between grid points dx, dy.
! For example: %d2f_dxdy (structure component) = d^2f/dxdz * dx * dy
!
! Input:
!   field_at_box   -- cmplx_field_at_2D_box_struct: Field and normalized derivatives at the 4 grid points.
!
! Output:
!   tri_coef        -- tricubic_cmplx_coef_struct: Coefficients.
!-

subroutine tricubic_interpolation_cmplx_coefs (field_at_box, tri_coef)

type (cmplx_field_at_3D_box_struct), target :: field_at_box
type (tricubic_cmplx_coef_struct) tri_coef
type (cmplx_field1_at_3D_pt_struct), pointer :: fg

complex(rp) field_array(64)
integer i, j, k, ijk

!

do i = 0, 1
do j = 0, 1
do k = 0, 1
  ijk = i + 2*j + 4*k + 1
  fg => field_at_box%pt(i,j,k)
  field_array(ijk:ijk+56:8) = [fg%f, fg%df_dx, fg%df_dy, fg%df_dz, fg%d2f_dxdy, fg%d2f_dxdz, fg%d2f_dydz, fg%d3f_dxdydz]
enddo
enddo
enddo

do i = 0, 3
do j = 0, 3
do k = 0, 3
  ijk = i + 4*j + 16*k + 1
  tri_coef%coef(i,j,k) = sum(m64_matrix%row(ijk)%pt(:) * field_array(:))
enddo
enddo
enddo

end subroutine tricubic_interpolation_cmplx_coefs

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Function tricubic_cmplx_eval (x_norm, y_norm, z_norm, tri_coef, df_dx, df_dy, df_dz) result (f_val)
!
! Routine to evaluate a tricubic interpolating complex function.
!
! Use the routine tricubic_interpolation_cmplx_coefs to generate tri_coef.
!
! Note: In the equations below, the eight points of the grid box being interpolated range 
! from (x0, y0, z0) to (x0+dx, y0+dy, z0+dz).
!
! Input:
!   x_norm    -- real(rp): x_norm = (x - x0) / dx 
!   y_norm    -- real(rp): y_norm = (y - y0) / dy
!   z_norm    -- real(rp): z_norm = (z - z0) / dz
!   tri_coef  -- tricubic_cmplx_coef_struct: Coefficients.
!
! Output:
!   f_val     -- complex(rp): Value of f.
!   df_dx     -- complex(rp): Normalized first derivative: True df/dx = df_dx * dx
!   df_dy     -- complex(rp): Normalized first derivative: True df/dy = df_dy * dy
!   df_dz     -- complex(rp): Normalized first derivative: True df/dz = df_dz * dz
!- 

function tricubic_cmplx_eval (x_norm, y_norm, z_norm, tri_coef, df_dx, df_dy, df_dz) result (f_val)

type (tricubic_cmplx_coef_struct) tri_coef

real(rp) x_norm, y_norm, z_norm
real(rp) x_power(0:3), y_power(0:3), z_power(0:3)

complex(rp) f_val
complex(rp), optional :: df_dx, df_dy, df_dz

integer i, j, k

! Init

x_power(0) = 1
y_power(0) = 1
z_power(0) = 1

do i = 1, 3
  x_power(i) = x_power(i-1) * x_norm
  y_power(i) = y_power(i-1) * y_norm
  z_power(i) = z_power(i-1) * z_norm
enddo

! Evaluate f_val

f_val = 0
do i = 0, 3
do j = 0, 3
do k = 0, 3
  f_val = f_val + x_power(i) * y_power(j) * z_power(k) * tri_coef%coef(i,j,k)
enddo
enddo
enddo

! Evaluate derivatives

if (present(df_dx)) then
  df_dx = 0
  do i = 1, 3
  do j = 0, 3
  do k = 0, 3
    df_dx = df_dx + i * x_power(i-1) * y_power(j) * z_power(k) * tri_coef%coef(i,j,k)
  enddo
  enddo
  enddo
endif

if (present(df_dy)) then
  df_dy = 0
  do i = 0, 3
  do j = 1, 3
  do k = 0, 3
    df_dy = df_dy + j * x_power(i) * y_power(j-1) * z_power(k) * tri_coef%coef(i,j,k)
  enddo
  enddo
  enddo
endif

if (present(df_dz)) then
  df_dz = 0
  do i = 0, 3
  do j = 0, 3
  do k = 1, 3
    df_dz = df_dz + k * x_power(i) * y_power(j) * z_power(k-1) * tri_coef%coef(i,j,k)
  enddo
  enddo
  enddo
endif

end function tricubic_cmplx_eval

end module 
