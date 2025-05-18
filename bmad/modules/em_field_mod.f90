!+
! Module em_field_mod
!
! Module to define the electric and magnetic fields for an elemet.
!-

module em_field_mod

use spline_mod
use taylor_mod

implicit none

contains

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!+
! Function g_bend_from_em_field (B, E, orbit) result (g_bend)
!
! Routine to calculate the bending strength (1/bending_radius) for a given particle for a given field.
! This will include the dipole bending field of an sbend.
!
! Input:
!   B(3)  -- real(rp): Magnetic field.
!   E(3)  -- real(rp): Electric field
!   orbit -- coord_struct: particle orbit
!
! Output:
!   g_bend(3) -- real(rp): bending strength vector.
!-

function g_bend_from_em_field (B, E, orbit) result (g_bend)

type (coord_struct) orbit
real(rp) b(3), e(3), g_bend(3)
real(rp) vel(3), rel_pc, force(3)

! vel is normalized velocity

rel_pc = 1 + orbit%vec(6)
vel(1:2) = [orbit%vec(2), orbit%vec(4)] / rel_pc
vel(3) = sqrt(1 - vel(1)**2 - vel(2)**2) * orbit%direction * orbit%time_dir

force = (E + cross_product(vel, B) * orbit%beta * c_light) * charge_of(orbit%species)
g_bend = -(force - vel * (dot_product(force, vel))) / (orbit%p0c * rel_pc)

end function g_bend_from_em_field

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine to_fieldmap_coords (ele, local_orb, s_body, ele_anchor_pt, r0, curved_ref_frame,
!                                                               x, y, z, cos_ang, sin_ang, err_flag)
!
! Routine to return the (x,y,s) position relative to a field map.
!
! Input:
!   ele               -- ele_struct: Element being tracked through.
!   local_orb         -- coord_struct: Particle orbit. Must be in local element coordinates.
!   s_body            -- real(rp): Longitudinal position relative to the entrance end of the element.
!   ele_anchor_pt     -- integer: anchor point of the field map (anchor_beginning$, anchor_center$, or anchor_end$).
!   r0(3)             -- real(rp): origin point of the fieldmap.
!   curved_ref_frame  -- logical: If the element is a bend: Does the field map follow the bend reference coords?
!
! Outpt:
!   x, y, z           -- real(rp): Coords relative to the field map.
!   cos_ang, sin_ang  -- real(rp): cos and sin of coordinate rotation angle.
!   err_flag          -- logical: Set True if there is an error. False otherwise.

subroutine to_fieldmap_coords (ele, local_orb, s_body, ele_anchor_pt, r0, curved_ref_frame, &
                                                                      x, y, z, cos_ang, sin_ang, err_flag)

type (ele_struct) ele
type (coord_struct) local_orb

real(rp) :: s_body, r0(3), x, y, z, x_save, s0, cos_ang, sin_ang
integer ele_anchor_pt
logical curved_ref_frame
logical :: err_flag

character(*), parameter :: r_name = 'to_fieldmap_coords'

!

err_flag = .false.

select case (ele_anchor_pt)
case (anchor_beginning$); s0 = 0
case (anchor_center$);    s0 = ele%value(l$) / 2
case (anchor_end$);       s0 = ele%value(l$)
case default
  call out_io (s_fatal$, r_name, 'BAD ELE_ANCHOR_PT FOR FIELD GRID IN ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
  err_flag = .true.
  return
end select

!

x = local_orb%vec(1)
z = s_body - s0

!

if ((ele%key == sbend$ .or. ele%key == rf_bend$) .and. ele%value(g$) /= 0 .and. .not. curved_ref_frame) then
  cos_ang = cos(z*ele%value(g$))
  sin_ang = sin(z*ele%value(g$))

  x_save = x
  x = (x_save + ele%value(rho$)) * cos_ang - ele%value(rho$)
  z = (x_save + ele%value(rho$)) * sin_ang 
endif

!

x = x - r0(1)
y = local_orb%vec(3) - r0(2)
z = z - r0(3)

end subroutine to_fieldmap_coords

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine rotate_em_field (field, w_mat, w_inv, calc_dfield, calc_potential)
!
! Routine to transform the fields using the given rotation matrices.
!
! Input:
!   field           -- em_field_struct: E and B fields and derivatives.
!   w_mat(3,3)      -- real(rp): rotation matrix.
!   w_inv(3,3)      -- real(rp): rotation matrix inverse = transpose(w_mat)
!   calc_dfield     -- Logical, optional: If present and True then rotate the field derivatives.
!   calc_potential  -- logical, optional: Rotate the magnetic vector potential? Default is false. 
!
! Output:
!   field           -- em_field_struct: E and B fields and derivatives.
!-

subroutine rotate_em_field (field, w_mat, w_inv, calc_dfield, calc_potential)

type (em_field_struct) field

real(rp) w_mat(3,3), w_inv(3,3)
logical, optional :: calc_dfield, calc_potential

!

field%B = matmul(w_mat, field%B)
field%E = matmul(w_mat, field%E)

if (logic_option(.false., calc_potential)) then
  field%A = matmul(w_mat, field%A)
endif

if (logic_option (.false., calc_dfield)) then
  field%dB = matmul(w_mat, matmul(field%dB, w_inv))
  field%dE = matmul(w_mat, matmul(field%dE, w_inv))
endif

end subroutine rotate_em_field

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine grid_field_interpolate (ele, orbit, grid, field, err_flag, x1, x2, x3, &
!                                                              allow_s_out_of_bounds, print_err)
!
! Subroutine to interpolate the E and B fields on a rectilinear grid.
!
! Input:
!   ele       -- ele_struct: Element containing the grid.
!   orbit     -- coord_struct: Used for constructing an error message if the particle is out of bounds.
!   grid      -- grid_field_struct: Grid to interpolate.
!   err_flag  -- Logical: Set to true if there is an error. False otherwise.
!   x1        -- real(rp): dimension 1 interpolation point.
!   x2        -- real(rp), optional: dimension 2 interpolation point.
!   x3        -- real(rp), optional: dimension 3 interpolation point.
!   allow_s_out_of_bounds -- logical, optional: allow s-coordinate grossly out of bounds to return
!                 zero field without an error. This is used when the field of one element overlaps
!                 the field of another. Default is False.
!   print_err -- logical, optional: print an error message if the particle is out of bounds? Default is True.
!
! Output:
!   field     -- grid_field_pt_struct: Interpolated field (complex)
!-

subroutine grid_field_interpolate (ele, orbit, grid, g_field, err_flag, x1, x2, x3, allow_s_out_of_bounds, print_err)

type (ele_struct) ele
type (coord_struct) orbit
type (grid_field_struct), target :: grid
type (grid_field_pt1_struct), intent(out) :: g_field
type (cmplx_field_at_2D_box_struct) field2_at_box
type (cmplx_field_at_3D_box_struct) field3_at_box
type (bicubic_cmplx_coef_struct), pointer :: bi_coef(:,:)    ! Save computed coefs for faster tracking
type (tricubic_cmplx_coef_struct), pointer :: tri_coef(:,:)  ! Save computed coefs for faster tracking

real(rp) :: x1
real(rp), optional :: x2, x3
real(rp) rel_x1, rel_x2, rel_x3, r2_x1

integer i, n, i1, i2, i3, n2, n3, grid_dim, allow_what, lbnd, ubnd, nn
integer, parameter :: allow_tiny$ = 1, allow_some$ = 2, allow_all$ = 3

logical err_flag, allow_out
logical, optional :: allow_s_out_of_bounds, print_err

character(*), parameter :: r_name = 'grid_field_interpolate'
character(40) extrapolation

! Pick appropriate dimension 

err_flag = .false.

allow_out = logic_option(.false., allow_s_out_of_bounds)
allow_what = allow_some$
if (allow_out) allow_what = allow_all$

grid_dim = grid_field_dimension(grid%geometry)

! xz grid

select case(grid_dim)
case (2)

  lbnd = lbound(grid%ptr%pt, 2); ubnd = ubound(grid%ptr%pt, 2)
  n3 = lbound(grid%ptr%pt, 3)

  call get_this_index(x2, 2, i2, rel_x2, err_flag, allow_what, .false.); if (err_flag) return
  ! If grossly out of longitudinal bounds just return zero field. Do not test transverse position in this case.
  if (i2 < lbnd - 1 .or. i2 > ubnd) return 

  call get_this_index(x1, 1, i1, rel_x1, err_flag, allow_tiny$, (.not. allow_out)); if (err_flag) return

  ! BiCubic interpolation

  if (grid%interpolation_order == 3) then
    ! Look for coefs already calculated
    n = size(grid%bi_coef, 1)
    do i = 1, n
      if (any(grid%bi_coef(i,1,1)%i_box /= [i1, i2])) cycle
      bi_coef => grid%bi_coef(i,:,:)
      exit
    enddo

    if (i == n+1) then
      if (i1 == 1) then
        extrapolation = 'SYMMETRIC:ZERO'
      else
        extrapolation = 'LINEAR:ZERO'
      endif

      grid%bi_coef(1:n-1,:,:) = grid%bi_coef(2:n,:,:)
      bi_coef => grid%bi_coef(4,:,:)

      do i = 1, 3
        call bicubic_compute_cmplx_field_at_2D_box(grid%ptr%pt(:,:,n3)%B(i), lbound(grid%ptr%pt), i1, i2, extrapolation, field2_at_box, err_flag)
        call bicubic_interpolation_cmplx_coefs (field2_at_box, bi_coef(1,i))
        call bicubic_compute_cmplx_field_at_2D_box(grid%ptr%pt(:,:,n3)%E(i), lbound(grid%ptr%pt), i1, i2, extrapolation, field2_at_box, err_flag)
        call bicubic_interpolation_cmplx_coefs (field2_at_box, bi_coef(2,i))
      enddo
    endif

    do i = 1, 3
      g_field%B(i) = bicubic_cmplx_eval(rel_x1, rel_x2, bi_coef(1,i))
      g_field%E(i) = bicubic_cmplx_eval(rel_x1, rel_x2, bi_coef(2,i))
    enddo

    return
  endif

  ! Do bilinear interpolation. If just outside longitudinally, interpolate between grid edge and zero.
  ! If using rotationally_symmetric_rz$ then the z component of the fields are even in r.
  ! In this case interpolate the z component using r^2 and not r.

  if (grid%geometry == rotationally_symmetric_rz$) then
    nn = 2
    r2_x1 = (2*i1*rel_x1 + rel_x1**2) / (2*i1 + 1)  ! = ((i1+r1)^2 - i1^2) / ((i1+1)^2 - i1^2)
  else
    nn = 3
  endif

  if (i2 == lbnd - 1 .or. i2 == ubnd) then  ! Just outside entrance end or just outside exit end
    if (i2 == lbnd - 1) then
      i2 = lbnd
      rel_x2 = 1 - rel_x2
    endif

    g_field%E(1:nn) = (1-rel_x1)*(1-rel_x2)   * grid%ptr%pt(i1,   i2, n3)%E(1:nn) &
                    + (rel_x1)*(1-rel_x2)     * grid%ptr%pt(i1+1, i2, n3)%E(1:nn) 

    g_field%B(1:nn) = (1-rel_x1)*(1-rel_x2)   * grid%ptr%pt(i1,   i2, n3)%B(1:nn) &
                    + (rel_x1)*(1-rel_x2)     * grid%ptr%pt(i1+1, i2, n3)%B(1:nn)  

    if (grid%geometry == rotationally_symmetric_rz$) then
      g_field%E(3) = (1-r2_x1)*(1-rel_x2)   * grid%ptr%pt(i1,   i2, n3)%E(3) &
                   + (r2_x1)*(1-rel_x2)     * grid%ptr%pt(i1+1, i2, n3)%E(3) 

      g_field%B(3) = (1-r2_x1)*(1-rel_x2)   * grid%ptr%pt(i1,   i2, n3)%B(3) &
                   + (r2_x1)*(1-rel_x2)     * grid%ptr%pt(i1+1, i2, n3)%B(3)  
    endif

  else  ! Inside
    g_field%E(1:nn) = (1-rel_x1)*(1-rel_x2) * grid%ptr%pt(i1,   i2,   n3)%E(1:nn) &
                    + (1-rel_x1)*(rel_x2)   * grid%ptr%pt(i1,   i2+1, n3)%E(1:nn) &
                    + (rel_x1)*(1-rel_x2)   * grid%ptr%pt(i1+1, i2,   n3)%E(1:nn) &
                    + (rel_x1)*(rel_x2)     * grid%ptr%pt(i1+1, i2+1, n3)%E(1:nn) 

    g_field%B(1:nn) = (1-rel_x1)*(1-rel_x2) * grid%ptr%pt(i1,   i2,   n3)%B(1:nn) &
                    + (1-rel_x1)*(rel_x2)   * grid%ptr%pt(i1,   i2+1, n3)%B(1:nn) &
                    + (rel_x1)*(1-rel_x2)   * grid%ptr%pt(i1+1, i2,   n3)%B(1:nn) &
                    + (rel_x1)*(rel_x2)     * grid%ptr%pt(i1+1, i2+1, n3)%B(1:nn)  

    if (grid%geometry == rotationally_symmetric_rz$) then
      g_field%E(3) = (1-r2_x1)*(1-rel_x2) * grid%ptr%pt(i1,   i2,   n3)%E(3) &
                   + (1-r2_x1)*(rel_x2)   * grid%ptr%pt(i1,   i2+1, n3)%E(3) &
                   + (r2_x1)*(1-rel_x2)   * grid%ptr%pt(i1+1, i2,   n3)%E(3) &
                   + (r2_x1)*(rel_x2)     * grid%ptr%pt(i1+1, i2+1, n3)%E(3) 

      g_field%B(3) = (1-r2_x1)*(1-rel_x2) * grid%ptr%pt(i1,   i2,   n3)%B(3) &
                   + (1-r2_x1)*(rel_x2)   * grid%ptr%pt(i1,   i2+1, n3)%B(3) &
                   + (r2_x1)*(1-rel_x2)   * grid%ptr%pt(i1+1, i2,   n3)%B(3) &
                   + (r2_x1)*(rel_x2)     * grid%ptr%pt(i1+1, i2+1, n3)%B(3)  
    endif
  endif

! xyz grid

case (3)

  lbnd = lbound(grid%ptr%pt, 3); ubnd = ubound(grid%ptr%pt, 3)
  n2 = lbound(grid%ptr%pt, 2);   n3 = lbound(grid%ptr%pt, 3)

  call get_this_index(x3, 3, i3, rel_x3, err_flag, allow_what, .false.); if (err_flag) return
  ! If grossly out of longitudinal bounds just return zero field. Do not test transverse position in this case.
  if (i3 < lbnd - 1 .or. i3 > ubnd) return 

  call get_this_index(x1, 1, i1, rel_x1, err_flag, allow_tiny$, (.not. allow_out)); if (err_flag) return
  call get_this_index(x2, 2, i2, rel_x2, err_flag, allow_tiny$, (.not. allow_out)); if (err_flag) return

  ! TriCubic interpolation

  if (grid%interpolation_order == 3) then
    ! Look for coefs already calculated
    n = size(grid%tri_coef, 1)
    do i = 1, n
      if (any(grid%tri_coef(i,n2,n3)%i_box /= [i1, i2, i3])) cycle
      tri_coef => grid%tri_coef(i,:,:)
      exit
    enddo

    if (i == n+1) then
      extrapolation = 'LINEAR:LINEAR:ZERO'

      grid%tri_coef(1:n-1,:,:) = grid%tri_coef(2:n,:,:)
      tri_coef => grid%tri_coef(4,:,:)

      do i = 1, 3
        call tricubic_compute_cmplx_field_at_3D_box(grid%ptr%pt%B(i), lbound(grid%ptr%pt), i1, i2, i3, extrapolation, field3_at_box, err_flag)
        call tricubic_interpolation_cmplx_coefs (field3_at_box, tri_coef(1,i))
        call tricubic_compute_cmplx_field_at_3D_box(grid%ptr%pt%E(i), lbound(grid%ptr%pt), i1, i2, i3, extrapolation, field3_at_box, err_flag)
        call tricubic_interpolation_cmplx_coefs (field3_at_box, tri_coef(2,i))
      enddo
    endif

    do i = 1, 3
      g_field%B(i) = tricubic_cmplx_eval(rel_x1, rel_x2, rel_x3, tri_coef(1,i))
      g_field%E(i) = tricubic_cmplx_eval(rel_x1, rel_x2, rel_x3, tri_coef(2,i))
    enddo

    return
  endif

  ! Do trilinear interpolation. If just outside longitudinally, interpolate between grid edge and zero.

  if (i3 == lbnd - 1 .or. i3 == ubnd) then  ! Just outside entrance end or just outside exit end
    if (i3 == lbnd - 1) then
      i3 = lbnd
      rel_x3 = 1 - rel_x3
    endif

    g_field%E(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3)   * grid%ptr%pt(i1,   i2,   i3)%E(:) &
                 + (1-rel_x1)*(rel_x2)  *(1-rel_x3)   * grid%ptr%pt(i1,   i2+1, i3)%E(:) &
                 + (rel_x1)  *(1-rel_x2)*(1-rel_x3)   * grid%ptr%pt(i1+1, i2,   i3)%E(:) &
                 + (rel_x1)  *(rel_x2)  *(1-rel_x3)   * grid%ptr%pt(i1+1, i2+1, i3)%E(:)               
               
    g_field%B(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3)   * grid%ptr%pt(i1,   i2,   i3)%B(:) &
                 + (1-rel_x1)*(rel_x2)  *(1-rel_x3)   * grid%ptr%pt(i1,   i2+1, i3)%B(:) &
                 + (rel_x1)  *(1-rel_x2)*(1-rel_x3)   * grid%ptr%pt(i1+1, i2,   i3)%B(:) &
                 + (rel_x1)  *(rel_x2)  *(1-rel_x3)   * grid%ptr%pt(i1+1, i2+1, i3)%B(:)

  else    ! Inside
    g_field%E(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1,   i2,   i3  )%E(:) &
                 + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1,   i2+1, i3  )%E(:) &
                 + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1+1, i2,   i3  )%E(:) &
                 + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1+1, i2+1, i3  )%E(:) &
                 + (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1,   i2,   i3+1)%E(:) &
                 + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1,   i2+1, i3+1)%E(:) &
                 + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1+1, i2,   i3+1)%E(:) &
                 + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1+1, i2+1, i3+1)%E(:)               
               
    g_field%B(:) = (1-rel_x1)*(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1,   i2,   i3  )%B(:) &
                 + (1-rel_x1)*(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1,   i2+1, i3  )%B(:) &
                 + (rel_x1)  *(1-rel_x2)*(1-rel_x3) * grid%ptr%pt(i1+1, i2,   i3  )%B(:) &
                 + (rel_x1)  *(rel_x2)  *(1-rel_x3) * grid%ptr%pt(i1+1, i2+1, i3  )%B(:) &
                 + (1-rel_x1)*(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1,   i2,   i3+1)%B(:) &
                 + (1-rel_x1)*(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1,   i2+1, i3+1)%B(:) &
                 + (rel_x1)  *(1-rel_x2)*(rel_x3)   * grid%ptr%pt(i1+1, i2,   i3+1)%B(:) &
                 + (rel_x1)  *(rel_x2)  *(rel_x3)   * grid%ptr%pt(i1+1, i2+1, i3+1)%B(:) 
  endif

case default
  call out_io (s_fatal$, r_name, 'BAD DIMENSION: \i0\ ', grid_field_dimension(grid%geometry))
  if (global_com%exit_on_error) call err_exit
  err_flag = .true.
  return
end select

!-------------------------------------------------------------------------------------
contains

subroutine get_this_index (x, ix_x, i0, rel_x0, err_flag, allow_what, lost_if_out_of_bounds)

type (coord_struct) orb2
real(rp) x, rel_x0, x_norm, x_diff, x_ave
integer ix_x, i0, ig0, ig1, allow_what
logical err_flag, lost_if_out_of_bounds
character(40) str

!

ig0 = lbound(grid%ptr%pt, ix_x)
ig1 = ubound(grid%ptr%pt, ix_x)

x_norm = x / grid%dr(ix_x)  ! Note that to_fieldmap_coords has already been called.
i0 = floor(x_norm)          ! index of lower 1 data point
rel_x0 = x_norm - i0        ! Relative distance from lower x1 grid point

! Out of bounds?

if (i0 < ig0 .or. i0 >= ig1) then
  g_field%E = 0
  g_field%B = 0

  select case (allow_what)
  case (allow_tiny$)
    ! Here do extrapolation is the point is within one dr/2 of the grid boundary.
    ! Otherwise it is an error.
    if (i0 == ig0 - 1 .and. rel_x0 > 0.5) then
      i0 = ig0
      rel_x0 = rel_x0 - 1
      return
    elseif (i0 == ig1 .and. rel_x0 < 0.5) then
      i0 = ig1 - 1
      rel_x0 = rel_x0 + 1
      return
    endif

  case (allow_some$)
    ! Here only generate an error message if the particle is grossly outside of the grid region.
    ! Here "gross" is defined as dOut > L_grid/2 where dOut is the distance between the
    ! particle and the grid edge and L_grid is the length of the grid.
    x_diff = (ig1 - ig0) * grid%dr(ix_x)
    x_ave = (ig1 + ig0) * grid%dr(ix_x) / 2
    if (abs(x - x_ave) < x_diff .or. i0 == ig0-1 .or. i0 == ig1) return

  case (allow_all$)
    return
  end select

  ! Avoid nedless error messages if the particle is outside the aperture.

  err_flag = .true.

  orb2%state = alive$
  if (ele%aperture_at == continuous$) then
    orb2 = orbit
    call check_aperture_limit(orb2, ele, in_between$, ele%branch%param)
  endif

  if (orb2%state == alive$ .and. logic_option(.true., print_err)) then
    if (lost_if_out_of_bounds) then
      str = 'PARTICLE WILL BE MARKED AS LOST'
    else
      str = 'SETTING FIELD TO ZERO'
    endif
    call out_io (s_error$, r_name, '\i0\D GRID_FIELD INTERPOLATION INDEX OUT OF BOUNDS: I\i0\ = \i0\ (POSITION = \f12.6\)', &
                                 'FOR ELEMENT: ' // ele%name // '  ' // trim(ele_loc_name(ele, parens = '()')), &
                                 'PARTICLE POSITION: \3F12.6\ ', str, &
                                 i_array = [grid_dim, ix_x, i0], r_array = [x, orbit%vec(1), orbit%vec(3), orbit%s-ele%s_start])
  endif

  if (lost_if_out_of_bounds .and. orbit%state == alive$) orbit%state = lost$
endif

end subroutine get_this_index 

end subroutine grid_field_interpolate

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Function field_interpolate_3d (position, field_mesh, deltas, position0) result (field)
!
! Function to interpolate a 3d field.
! The interpolation is such that the derivative is continuous.
!
! Note: For "interpolation" outside of the region covered by the field_mesh
! it is assumed that the field is constant, Equal to the field at the
! boundary.
!
! Input:
!   position(3)       -- Real(rp): (x, y, z) position.
!   field_mesh(:,:,:) -- Real(rp): Grid of field points.
!   deltas(3)         -- Real(rp): (dx, dy, dz) distances between mesh points.
!   position0(3)      -- Real(rp), optional:  position at (ix0, iy0, iz0) where
!                            (ix0, iy0, iz0) is the lower bound of the
!                            filed_mesh(i, j, k) array. If not present then
!                            position0 is taken to be (0.0, 0.0, 0.0)
! Output:
!   field -- Real(rp): interpolated field.
!-

function field_interpolate_3d (position, field_mesh, deltas, position0) result (field)

real(rp), optional, intent(in) :: position0(3)
real(rp), intent(in) :: position(3), field_mesh(0:,0:,0:), deltas(3)
real(rp) field

real(rp) r(3), f(-1:2), g(-1:2), h(-1:2), r_frac(3)

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
field = interpolate_1d (r_frac(1), h)

!---------------------------------------------------------------
contains

! interpolation in 1 dimension using 4 equally spaced points: P1, P2, P3, P4.
!   x = interpolation point.
!           x = 0 -> point is at P2.
!           x = 1 -> point is at P3.
! Interpolation is done so that the derivative is continuous.
! The interpolation uses a cubic polynomial

function interpolate_1d (x, field1_in) result (field1)

real(rp) field1, x, field1_in(4), df_2, df_3
real(rp) c0, c1, c2, c3

!

df_2 = (field1_in(3) - field1_in(1)) / 2   ! derivative at P2
df_3 = (field1_in(4) - field1_in(2)) / 2   ! derivative at P3

c0 = field1_in(2)
c1 = df_2
c2 = 3 * field1_in(3) - df_3 - 3 * field1_in(2) - 2 * df_2
c3 = df_3 - 2 * field1_in(3) + 2 * field1_in(2) + df_2

field1 = c0 + c1 * x + c2 * x**2 + c3 * x**3

end function interpolate_1d

end function field_interpolate_3d 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine em_field_derivatives (ele, param, s_pos, orbit, local_ref_frame, dfield, grid_allow_s_out_of_bounds, rf_time)
!
! Routine to calculate field derivatives.
! In theory this should be handled by em_filed_calc. In practice, em_field_calc is currently incomplete.
!
! Input
!   ele             -- Ele_struct: Element
!   param           -- lat_param_struct: Lattice parameters.
!   s_pos           -- Real(rp): Longitudinal position relative to the upstream edge of the element.
!   time            -- Real(rp): Particle time.
!                       For absolute time tracking this is the absolute time.
!                       For relative time tracking this is relative to the reference particle entering the element.
!   orbit           -- Coord_struct: Transverse coordinates.
!     %vec(1), %vec(3)  -- Transverse coords. These are the only components used in the calculation.
!   local_ref_frame     -- Logical, If True then take the input coordinates and output fields 
!                                   as being with respect to the frame of referene of the element (ignore misalignments). 
!   grid_allow_s_out_of_bounds 
!                    -- logical, optional: For grids, allow s-coordinate to be grossly out of bounds 
!                         and return zero instead of an error? Default: False. Used internally for overlapping fields.
!   rf_time          -- real(rp), optional: RF clock time. If not present then the time will be calculated using the standard algorithm.
!
! Output:
!   dfield       -- em_field_struct: E and B field derivatives. dfield%E and dfield%B are not touched.
!-

subroutine em_field_derivatives (ele, param, s_pos, orbit, local_ref_frame, dfield, grid_allow_s_out_of_bounds, rf_time)

type (ele_struct), target :: ele
type (lat_param_struct) param
type (em_field_struct) :: dfield, f0, f1
type (coord_struct) :: orbit, orb, local_orb

real(rp), optional :: rf_time
real(rp) s_pos, s0, s1, del, s_body
logical local_ref_frame
logical, optional :: grid_allow_s_out_of_bounds

!

local_orb = orbit
if (local_ref_frame) then
  s_body = s_pos
else
  call offset_particle (ele, set$, local_orb, set_hvkicks = .false., s_pos = s_pos, s_out = s_body)
endif

!

orb = local_orb
del = bmad_com%d_orb(1)

orb%vec(1) = local_orb%vec(1) - del
call em_field_calc (ele, param, s_pos, orb, .true., f0, &
                            grid_allow_s_out_of_bounds = grid_allow_s_out_of_bounds, rf_time = rf_time)

orb%vec(1) = local_orb%vec(1) + del
call em_field_calc (ele, param, s_pos, orb, .true., f1, &
                            grid_allow_s_out_of_bounds = grid_allow_s_out_of_bounds, rf_time = rf_time)

dfield%dB(:,1) = (f1%B - f0%B) / (2 * del)
dfield%dE(:,1) = (f1%E - f0%E) / (2 * del)

!

orb = local_orb
del = bmad_com%d_orb(3)

orb%vec(3) = local_orb%vec(3) - del
call em_field_calc (ele, param, s_pos, orb, .true., f0, &
                            grid_allow_s_out_of_bounds = grid_allow_s_out_of_bounds, rf_time = rf_time)

orb%vec(3) = local_orb%vec(3) + del
call em_field_calc (ele, param, s_pos, orb, .true., f1, &
                            grid_allow_s_out_of_bounds = grid_allow_s_out_of_bounds, rf_time = rf_time)

dfield%dB(:,2) = (f1%B - f0%B) / (2 * del)
dfield%dE(:,2) = (f1%E - f0%E) / (2 * del)

!

orb = local_orb
del = bmad_com%d_orb(5)

s0 = max(0.0_rp, s_pos-del)
s1 = min(ele%value(l$), s_pos+del)
if (s1 == s0) return  ! Cannot calc if zero length

call em_field_calc (ele, param, s0, local_orb, .true., f0, &
                            grid_allow_s_out_of_bounds = grid_allow_s_out_of_bounds, rf_time = rf_time)

call em_field_calc (ele, param, s1, local_orb, .true., f1, &
                            grid_allow_s_out_of_bounds = grid_allow_s_out_of_bounds, rf_time = rf_time)

dfield%dB(:,3) = (f1%B - f0%B) / (s1 - s0)
dfield%dE(:,3) = (f1%E - f0%E) / (s1 - s0)

!

if (.not. local_ref_frame) call convert_field_ele_to_lab (ele, s_body, .true., dfield, .true.)

end subroutine em_field_derivatives

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Fuction gen_grad_field (deriv, gg, rho, theta) result (field)
!
! Routine to calculate the field from a generalized gradient.
!
! Input:
!   deriv(0:)     -- real(rp): Array of derivatives.
!   gg            -- gen_grad1_struct: Contains m, sincos, and n_deriv_max parameters. 
!   rho, theta    -- real(rp): Particle transverse position in cylindrical coords.
!
! Output:
!   field(3)      -- real(rp): Field.
!-

function gen_grad_field (deriv, gg, rho, theta) result (field)

type (gen_grad1_struct) :: gg

real(rp) deriv(0:), field(3), rho, theta
real(rp) cd, f, ff, F_rho, F_theta, p_rho

integer id, nn
logical is_even

!

field = 0

do id = 0, gg%n_deriv_max
  cd = deriv(id)
  if (cd == 0) cycle

  is_even = (modulo(id,2) == 0)
  if (is_even) then
    nn = id / 2
  else
    nn = (id - 1) / 2
  endif

  f = (-0.25_rp)**nn * factorial(gg%m) / (factorial(nn) * factorial(nn+gg%m))

  if (id+gg%m-1 <= 0) then  ! Covers case where rho = 0
    p_rho = 1
  else
    p_rho = rho**(id+gg%m-1)
  endif

  ff = f * p_rho * cd

  if (is_even) then
    if (gg%sincos == sin$) then
      F_rho   = ff * (2*nn+gg%m) * sin(gg%m*theta)
      F_theta = ff * gg%m * cos(gg%m*theta)
    else
      F_rho   = ff * (2*nn+gg%m) * cos(gg%m*theta)
      F_theta = -ff * gg%m * sin(gg%m*theta)
    endif
    field(1) = field(1) + F_rho * cos(theta) - F_theta * sin(theta)
    field(2) = field(2) + F_rho * sin(theta) + F_theta * cos(theta)

  else
    if (gg%sincos == sin$) then
      field(3) = field(3) + ff * sin(gg%m*theta)
    else
      field(3) = field(3) + ff * cos(gg%m*theta)
    endif
  endif

enddo

end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine convert_field_ele_to_lab (ele, s_here, forward_transform, field, calc_dfield, calc_potential)
!
! Convert fields: ele to lab coords
!
! Input:
!   ele                 -- Ele_struct: Lattice element.
!   s_here              -- real(rp) S-position.
!   forward_transform   -- Transform foward (to lab) or reverse.
!   calc_dfield         -- Logical, optional: If present and True then calculate the field derivatives.
!   calc_potential      -- logical, optional: Calc electric and magnetic potentials? Default is false. 
!                         This is experimental and only implemented for wigglers at present.
!
! Output:
!   field               -- em_field_struct: EM field.
!-

subroutine convert_field_ele_to_lab (ele, s_here, forward_transform, field, calc_dfield, calc_potential)

type (ele_struct) ele
type (em_field_struct) field

real(rp) s_here, w_mat(3,3), w_inv(3,3), w_s(3,3), w_rt(3,3), w_rt_inv(3,3)
real(rp) theta
logical forward_transform
logical, optional :: calc_dfield, calc_potential

!

if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(roll$), w_mat)
  theta = ele%value(g$) * s_here - ele%value(angle$)/2
  w_s = w_mat_for_x_pitch (theta)
  if (ele%value(ref_tilt_tot$) == 0) then
    w_mat = matmul(matmul(w_s, w_mat), transpose(w_s))
  else
    w_rt = w_mat_for_tilt (ele%value(ref_tilt_tot$))
    w_rt_inv = w_mat_for_tilt (ele%value(ref_tilt_tot$), .true.)
    w_mat = matmul(matmul(matmul(matmul(matmul(w_rt, w_s), w_rt_inv), w_mat), w_rt), transpose(w_s))
  endif
  w_inv = transpose(w_mat)
else
  call floor_angles_to_w_mat (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%value(tilt_tot$), w_mat, w_inv)
endif

if (forward_transform) then
  call rotate_em_field (field, w_mat, w_inv, calc_dfield, calc_potential)
else
  call rotate_em_field (field, w_inv, w_mat, calc_dfield, calc_potential)
endif

end subroutine convert_field_ele_to_lab

end module
