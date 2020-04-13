!+
! Subroutine track_a_converter (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an converter element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: converter element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_converter (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_converter
use random_mod
use super_recipes_mod

implicit none

! This common struct is needed to get around ifort bug where debug info (when a program is compiled with debug) is 
! not generated for variables that are used in both the main and a contained routine.
type outgoing_particle_struct
  real(rp) pc_out, r
  real(rp) dx_ds, dy_ds
  logical :: lost = .false.
end type

type this_common_struct
  type (ele_struct), pointer :: ele
  real(rp) pc_out, r
  real(rp) pc_out_min, pc_out_max
  real(rp) dxy_ds_max
  real(rp) A_dir
  real(rp) beta, alpha_x, alpha_y, c_x
end type

type (this_common_struct) com
type (outgoing_particle_struct) out1, out2
type (coord_struct) :: orbit, orb0
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (converter_struct), pointer :: conv

real(rp), optional :: mat6(6,6)
real(rp) r_ran(3), pc_in, thickness, drd, drd2, azimuth_angle

integer i, nd, ix

logical, optional :: make_matrix
logical thickness_interpolate, err_flag

character(*), parameter :: r_name = 'track_a_converter'

!

conv => ele%converter
nd = size(conv%dist)
if (nd == 0) then
  call out_io (s_error$, r_name, 'CONVERTER ELEMENT DOES NOT HAVE ANY ASSOCIATED DISTRIBUTION(S) FOR: ' // ele%name, &
                                 'PARTICLE WILL BE MARKED AS LOST.')
  orbit%state = lost$
  return
endif

if (.not. allocated(ele%converter%dist(1)%sub_dist(1)%prob_pc_r%p_norm)) call probability_tables_setup(ele, err_flag)

call offset_particle (ele, param, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

if (ele%value(l$) == 0) then
  if (size(conv%dist) > 1) then
    call out_io (s_error$, r_name, 'CONVERTER ELEMENT HAS MULTIPLE DISTRIBUTIONS WITH MULTIPLE THICKNESSES BUT', &
                                   'NO LENGTH HAS BEEN SET. FOR ELEMENT: ' // ele%name, &
                                   'PARTICLE WILL BE MARKED AS LOST.')
    orbit%state = lost$
    return
  endif
endif

! Find dists that straddle thickness

pc_in = orbit%p0c * (1 + orbit%vec(6))
call ran_uniform(r_ran)
thickness_interpolate = .false.

if (size(conv%dist) == 1) then
  call calc_out_coords(ele, conv%dist(1), pc_in, r_ran, out1)
  thickness = conv%dist(1)%thickness

else
  if (ele%value(l$) < conv%dist(1)%thickness .or. ele%value(l$) > conv%dist(nd)%thickness) then
    call out_io (s_error$, r_name, 'CONVERTER LENGTH OUT OF RANGE OF THICKNESSES SET IN THE DISTRIBUTIONS FOR: '// ele%name, &
                                   'PARTICLE WILL BE MARKED AS LOST.')
    orbit%state = lost$
    return
  endif

  thickness = ele%value(l$)
  ix = bracket_index(ele%value(l$), conv%dist%thickness, 1, dr = drd)
  call calc_out_coords(ele, conv%dist(ix), pc_in, r_ran, out1)
  if (ix /= nd) then
    thickness_interpolate = .true.
    call calc_out_coords(ele, conv%dist(ix+1), pc_in, r_ran, out2)
  endif
endif

!

if (thickness_interpolate) then
  out1%lost = (out1%lost .or. out2%lost)
  drd2 = 1 - drd
  out1%pc_out = drd * out1%pc_out + drd2 * out2%pc_out
  out1%r      = drd * out1%r      + drd2 * out2%r
  out1%dx_ds  = drd * out1%dx_ds  + drd2 * out2%dx_ds
  out1%dy_ds  = drd * out1%dy_ds  + drd2 * out2%dy_ds
endif

if (out1%lost) then
  orbit%state = lost$
  return
endif

azimuth_angle = twopi * r_ran(3)

orb0 = orbit
orbit%species = conv%species_out
orbit%p0c = ele%value(p0c$)
orbit%vec(6) = out1%pc_out / ele%value(p0c$) - 1
call convert_pc_to (out1%pc_out, orbit%species, beta = orbit%beta)

orbit%vec(1) = orbit%vec(1) + out1%r * cos(azimuth_angle)
orbit%vec(2) = orbit%vec(2) + (out1%dx_ds * cos(azimuth_angle) - out1%dy_ds * sin(azimuth_angle)) * (1 + orbit%vec(6))
orbit%vec(3) = orbit%vec(3) + out1%r * sin(azimuth_angle)
orbit%vec(4) = orbit%vec(4) + (out1%dx_ds * sin(azimuth_angle) + out1%dy_ds * cos(azimuth_angle)) * (1 + orbit%vec(6))
orbit%vec(5) = orbit%vec(5) * orb0%beta / orbit%beta

call offset_particle (ele, param, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

!------------------------------------------------
contains

! Calculate output coords for a given thickness.

subroutine calc_out_coords(ele, dist, pc_in, r_ran, out)

type (ele_struct) ele
type (converter_distribution_struct) dist
type (outgoing_particle_struct) out, out1, out2

real(rp) pc_in, r_ran(:), rsd, rsd2
integer nsd, ix_sd

! Interpolate between energies.

nsd = size(dist%sub_dist)
if (pc_in < dist%sub_dist(1)%pc_in .or. pc_in >= dist%sub_dist(nsd)%pc_in) then
  call out_io (s_fatal$, r_name, 'INCOMING PARTICLE ENERGY OUT OF RANGE: \es12.4\ FOR CONVERTER: ' // ele%name, &
                                 'PARTICLE WILL BE MARKED AS LOST.')
  out%lost = .true.
  return
endif

ix_sd = bracket_index(pc_in, dist%sub_dist%pc_in, 1, rsd)
call calc_out_coords2 (ele, dist%sub_dist(ix_sd), r_ran, out1)
call calc_out_coords2 (ele, dist%sub_dist(ix_sd+1), r_ran, out2)

rsd2 = 1 - rsd
out%pc_out = rsd * out1%pc_out + rsd2 * out2%pc_out
out%r     = rsd * out1%r     + rsd2 * out2%r
out%dx_ds = rsd * out1%dx_ds + rsd2 * out2%dx_ds
out%dy_ds = rsd * out1%dy_ds + rsd2 * out2%dy_ds

end subroutine calc_out_coords

!------------------------------------------------
! contains

! Calculate output coords for a given thickness and a given incoming particle energy.

subroutine calc_out_coords2(ele, sub_dist, r_ran, out)

type (ele_struct) ele
type (converter_sub_distribution_struct), target :: sub_dist
type (outgoing_particle_struct) out
type (converter_prob_pc_r_struct), pointer :: ppcr

real(rp) r_ran(:), dr
integer ix
logical err_flag

!

ppcr => sub_dist%prob_pc_r
ix = bracket_index(r_ran(1), ppcr%integ_pc_out, 1, dr)


end subroutine calc_out_coords2

!------------------------------------------------
! contains

subroutine probability_tables_setup(ele, err_flag)

type (ele_struct), target :: ele
type (converter_struct), pointer :: conv
type (converter_distribution_struct) d_temp
type (converter_distribution_struct), pointer :: dist
type (converter_direction_out_struct), pointer :: dir_out
type (converter_prob_pc_r_struct), pointer :: ppcr

real(rp) dxy_ds_max, angle_max, prob
real(rp) pc_out, r, A_dir
integer id, isd, ie, ir, ne, nr
logical err_flag, ordered

! Order distributions in thickness

conv => ele%converter

ordered  = .true.
do
  do id = 2, size(conv%dist)
    if (conv%dist(id-1)%thickness > conv%dist(id)%thickness) then
      ordered = .false.
      d_temp = conv%dist(id-1)
      conv%dist(id-1) = conv%dist(id)
      conv%dist(id) = d_temp
    endif
  enddo
  if (ordered) exit
enddo

!

do id = 1, size(conv%dist)
  dist =>conv%dist(id)
  dxy_ds_max = dist%dxy_ds_max
  if (ele%value(angle_out_max$) /= 0) angle_max = min(ele%value(angle_out_max$), angle_max)
  do isd = 1, size(dist%sub_dist)
    dir_out => dist%sub_dist(isd)%dir_out

    ppcr => dist%sub_dist(isd)%prob_pc_r
    ne = size(ppcr%pc_out)
    nr = size(ppcr%r)

    if (.not. allocated(ppcr%p_norm)) then
      allocate (ppcr%p_norm(ne,nr))
    endif

    ppcr%pc_out_max = ppcr%pc_out(ne)
    if (ele%value(pc_out_max$) /= 0) ppcr%pc_out_max = min(ppcr%pc_out_max, ele%value(pc_out_max$))
    ppcr%pc_out_min = max(ppcr%pc_out(1), ele%value(pc_out_min$))
    ppcr%dxy_ds_max = min(dist%dxy_ds_max, atan(ele%value(angle_out_max$)))
    if (ppcr%pc_out_max <= ppcr%pc_out_min) then
      call out_io(s_fatal$, r_name, 'PC_OUT_MAX IS LESS THAN OR EQUAL TO PC_OUT_MIN. FOR ELEMENT: ' // ele%name, &
                                    'PARTICLE WILL BE MARKED AS LOST.')
      err_flag = .true.
    endif

    do ie = 1, ne
      pc_out = ppcr%pc_out(ie)
      do ir = 1, nr
        r = ppcr%r(ir)
        !!! com = calc_angle_coefs(ele, pc_out, r, 1.0_rp)
        !!! prob = super_qromb(p1_angle, -dxy_ds_max, dxy_ds_max, 1e-5_rp, 0.0_rp, 5, err_flag)
        A_dir = 1 / prob   ! A_dir is not normalized by angle_max.
        !!! com = this_common_struct(ele, pc_out, r, A_dir)
        !!! prob = super_qromb(p1_angle_func, -com%dxy_ds_max, com%dxy_ds_max, 1e-5_rp, 0.0_rp, 5, err_flag)
        ppcr%p_norm(ie,ir) = ppcr%prob(ie,ir) * prob  ! p_norm is now normalized by angle_max
      enddo
    enddo

    ppcr%integ_pc_out(1) = 0
    do ie = 2, ne
      !!! ppcr%integ_pc_out(ie) = ppcr%integ_pc_out(ie-1) + super_qromb_2D(p2_norm_func, ppcr%pc_out(ie-1), ppcr%pc_out(ie), &
      !!!                                                               0.0_rp, ppcr%r_out(nr), 1e-4, 0.0_rp, 2, err_flag)
    enddo
    ppcr%integ_pc_out = ppcr%integ_pc_out(ne)
    ppcr%integ_pc_out = ppcr%integ_pc_out / ppcr%integrated_prob

    ppcr%integ_r(:,1) = 0
    do ir = 2, nr
      do ie = 1, ne
        !!! ppcr%integ_r(ie,ir) = ppcr%integ_r(ie,ir-1) + super_qromb(p1_norm_func, ppcr%r(ir-1), ppcr%r(ir), &
        !!!                                                                             1e-4, 0.0_rp, 2, err_flag)
      enddo
    enddo

    do ie = 1, ne
      ppcr%integ_r(ie,:) = ppcr%integ_r(ie,:) / ppcr%integ_r(ie,nr)
    enddo

  enddo
enddo

end subroutine probability_tables_setup

!------------------------------------------------
! contains

subroutine p1_angle (x, value)

real(rp),intent(in) :: x(:)
real(rp) value(:)

value = 1

end subroutine p1_angle

end subroutine
