!+
! Subroutine track_a_converter (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an converter element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: converter element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_converter (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_converter
use random_mod
use super_recipes_mod
use spline_mod

implicit none

integer, parameter :: n_pt$ = 100

type converter_param_storage_struct
  real(rp) pc_out, r
  real(rp) rel_spin_z
  real(rp) dxds, dyds
  real(rp) weight
  logical :: lost = .false.
  real(rp) beta, alpha_x, alpha_y, c_x, integ_prob_tot
  real(rp) dxds_min, dxds_max, dyds_max
end type

type converter_common_struct
  type (converter_prob_pc_r_struct), pointer :: ppcr
  type (converter_distribution_struct), pointer :: dist
  type (spline_struct) dxds_spline(0:n_pt$)
  real(rp) dxds_integ(0:n_pt$)
  real(rp) r_ran, integ_prob_tot
  integer ipc
end type

type (converter_common_struct), target :: com
type (converter_common_struct), pointer :: com_ptr  ! Used to get around ifort problem with debugging code.
type (converter_param_storage_struct) out1, out2
type (coord_struct) :: orbit, orb0
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (converter_struct), pointer :: conv

real(rp), optional :: mat6(6,6)
real(rp) r_ran(5), pc_in, thickness, drd, drd2, azimuth_angle, ps

integer i, nd, ix

logical, optional :: make_matrix
logical thickness_interpolate, err_flag

character(*), parameter :: r_name = 'track_a_converter'

!

com_ptr => com
conv => ele%converter
nd = size(conv%dist)
if (nd == 0) then
  call out_io (s_error$, r_name, 'CONVERTER ELEMENT DOES NOT HAVE ANY ASSOCIATED DISTRIBUTION(S) FOR: ' // ele%name, &
                                 'PARTICLE WILL BE MARKED AS LOST.')
  orbit%state = lost$
  return
endif

if (ele%value(l$) == 0) then
  if (size(conv%dist) > 1) then
    call out_io (s_error$, r_name, 'CONVERTER ELEMENT HAS MULTIPLE DISTRIBUTIONS WITH MULTIPLE THICKNESSES BUT', &
                                   'NO LENGTH HAS BEEN SET. FOR ELEMENT: ' // ele%name, &
                                   'PARTICLE WILL BE MARKED AS LOST.')
    orbit%state = lost$
    return
  endif
endif

if (logic_option(.false., make_matrix)) call mat_make_unit(mat6)

call offset_particle (ele, set$, orbit)

! Find dists that straddle thickness

if (.not. allocated(ele%converter%dist(1)%sub_dist(1)%prob_pc_r%p_norm)) then
  call probability_tables_setup(ele, err_flag)
  if (err_flag) goto 9000
endif

pc_in = orbit%p0c * (1 + orbit%vec(6))
call ran_uniform(r_ran)
thickness_interpolate = .false.

if (size(conv%dist) == 1) then
  call calc_out_coords(ele, conv%dist(1), pc_in, r_ran, out1, err_flag);  if (err_flag) goto 9000
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
  call calc_out_coords(ele, conv%dist(ix), pc_in, r_ran, out1, err_flag);  if (err_flag) goto 9000
  if (ix /= nd) then
    thickness_interpolate = .true.
    call calc_out_coords(ele, conv%dist(ix+1), pc_in, r_ran, out2, err_flag);  if (err_flag) goto 9000
  endif
endif

!

if (thickness_interpolate) then
  out1%lost = (out1%lost .or. out2%lost)
  drd2 = 1 - drd
  out1%pc_out = drd2 * out1%pc_out + drd * out2%pc_out
  out1%r      = drd2 * out1%r      + drd * out2%r
  out1%dxds   = drd2 * out1%dxds   + drd * out2%dxds
  out1%dyds   = drd2 * out1%dyds   + drd * out2%dyds
  out1%weight = drd2 * out1%weight + drd * out2%weight
  out1%rel_spin_z = drd2 * out1%rel_spin_z + drd * out2%rel_spin_z
endif

if (out1%lost) then
  orbit%state = lost$
  return
endif

call track_a_drift (orbit, ele%value(l$))

azimuth_angle = twopi * r_ran(5)

orb0 = orbit
orbit%species = conv%species_out
orbit%p0c = ele%value(p0c$)
orbit%vec(6) = out1%pc_out / ele%value(p0c$) - 1
call convert_pc_to (out1%pc_out, orbit%species, beta = orbit%beta)

ps = sqrt((1 + orbit%vec(6))**2 / (1 + out1%dxds**2 + out1%dyds**2))
orbit%charge = out1%weight
orbit%vec(1) = orbit%vec(1) + out1%r * cos(azimuth_angle)
orbit%vec(2) = orbit%vec(2) + (out1%dxds * cos(azimuth_angle) - out1%dyds * sin(azimuth_angle)) * ps
orbit%vec(3) = orbit%vec(3) + out1%r * sin(azimuth_angle)
orbit%vec(4) = orbit%vec(4) + (out1%dxds * sin(azimuth_angle) + out1%dyds * cos(azimuth_angle)) * ps
orbit%vec(5) = orbit%vec(5) * orb0%beta / orbit%beta
orbit%spin   = [0, 0, 1] * orbit%spin(3) * out1%rel_spin_z

call offset_particle (ele, unset$, orbit)
return

! Error detected

9000 continue
orbit%state = lost$
return

!------------------------------------------------
contains

! Calculate output coords for a given thickness.

subroutine calc_out_coords(ele, dist, pc_in, r_ran, out, err_flag)

type (ele_struct) ele
type (converter_distribution_struct) dist
type (converter_param_storage_struct) out, out1, out2

real(rp) pc_in, r_ran(:), rsd, rsd2
integer nsd, ix_sd
logical err_flag

! Interpolate between energies.

nsd = size(dist%sub_dist)
if (pc_in < dist%sub_dist(1)%pc_in .or. pc_in >= dist%sub_dist(nsd)%pc_in) then
  call out_io (s_fatal$, r_name, 'INCOMING PARTICLE MOMENTUM \es12.4\ OUT OF THE DISTRIBUTION RANGE: [\es10.4\, \es10.4\]', &
                                 'FOR CONVERTER: ' // ele%name, &
                                 'PARTICLE WILL BE MARKED AS LOST.', r_array = [pc_in, dist%sub_dist(1)%pc_in, dist%sub_dist(nsd)%pc_in])
  out%lost = .true.
  return
endif

ix_sd = bracket_index(pc_in, dist%sub_dist%pc_in, 1, rsd, restrict = .true.)
call calc_out_coords2 (ele, dist, dist%sub_dist(ix_sd), r_ran, out1, err_flag);  if (err_flag) return
call calc_out_coords2 (ele, dist, dist%sub_dist(ix_sd+1), r_ran, out2, err_flag);  if (err_flag) return

rsd2 = 1 - rsd
out%pc_out = rsd2 * out1%pc_out + rsd * out2%pc_out
out%r      = rsd2 * out1%r      + rsd * out2%r
out%dxds   = rsd2 * out1%dxds   + rsd * out2%dxds
out%dyds   = rsd2 * out1%dyds   + rsd * out2%dyds
out%weight = rsd2 * dist%sub_dist(ix_sd)%prob_pc_r%integrated_prob + rsd * dist%sub_dist(ix_sd+1)%prob_pc_r%integrated_prob 

end subroutine calc_out_coords

!------------------------------------------------
! contains

! Calculate output coords for a given thickness and a given incoming particle energy.

subroutine calc_out_coords2(ele, dist, sub_dist, r_ran, out, err_flag)

type (ele_struct) ele
type (converter_distribution_struct), target :: dist
type (converter_sub_distribution_struct), target :: sub_dist
type (converter_param_storage_struct) out
type (converter_prob_pc_r_struct), pointer :: ppcr

real(rp) r_ran(:), dpc, dpc2, dr, rx, delta, dx, k_const, b
integer ix, ix_pc, ix_r, n, status

logical err_flag

! pc_out calc

ppcr => sub_dist%prob_pc_r
ix_pc = bracket_index(r_ran(1), ppcr%integ_pc_out, 1, dpc, restrict = .true.)
out%pc_out = (1-dpc) * ppcr%pc_out(ix_pc) + dpc * ppcr%pc_out(ix_pc+1)

! r calc

ppcr%integ_r_ave = (1-dpc) * ppcr%integ_r(ix_pc,:) + dpc * ppcr%integ_r(ix_pc+1,:)
ix_r = bracket_index(r_ran(2), ppcr%integ_r_ave(:), 1, dr, restrict = .true.)
out%r = (1-dr) * ppcr%r(ix_r) + dr * ppcr%r(ix_r+1)

! spin calc

out%rel_spin_z = (1-dpc) * ((1-dr) * ppcr%spin_z(ix_pc, ix_r) + dr * ppcr%spin_z(ix_pc, ix_r+1)) + &
                    dpc * ((1-dr) * ppcr%spin_z(ix_pc+1, ix_r) + dr * ppcr%spin_z(ix_pc+1, ix_r+1))

! dx/ds calc

call calc_dir_out_params(ele, dist, sub_dist, .true., out, err_flag);  if (err_flag) return

com%r_ran = r_ran(3)
com%integ_prob_tot = out%integ_prob_tot
out%dxds = super_zbrent(dxds_func, out%dxds_min, out%dxds_max, 0.0_rp, 1e-4_rp, status)

! dy/ds calc

b = sqrt(1 + (out%alpha_x * (out%dxds - out%c_x))**2)
k_const = atan(out%alpha_y * out%dyds_max / b)
out%dyds = b * tan((2 * r_ran(4) - 1) * k_const) / out%alpha_y

end subroutine calc_out_coords2

!------------------------------------------------
! contains

function dxds_func (x, status) result (value)

real(rp), intent(in) :: x
real(rp) value
integer status, ix

!

ix = bracket_index(x, com%dxds_spline%x0, 0, restrict = .true.)
value = (com%dxds_integ(ix) + spline1(com%dxds_spline(ix), x, -1)) - com%integ_prob_tot * com%r_ran

end function dxds_func

!------------------------------------------------
! contains

subroutine probability_tables_setup (ele, err_flag)

type (ele_struct), target :: ele
type (converter_struct), pointer :: conv
type (converter_distribution_struct) d_temp
type (converter_distribution_struct), pointer :: dist
type (converter_sub_distribution_struct), pointer :: sub_dist
type (converter_direction_out_struct), pointer :: dir_out
type (converter_prob_pc_r_struct), pointer :: ppcr
type (converter_param_storage_struct) out

real(rp) prob, unnorm_integ_prob_tot
integer id, isd, ipc, ir, npc, nr
logical err_flag, ordered

! Order distributions in thickness

err_flag = .true.
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

  do isd = 1, size(dist%sub_dist)
    sub_dist => dist%sub_dist(isd)
    dir_out => sub_dist%dir_out

    com%dist => dist
    com%ppcr => sub_dist%prob_pc_r
    ppcr => sub_dist%prob_pc_r
    npc = size(ppcr%pc_out)
    nr = size(ppcr%r)

    if (.not. allocated(ppcr%p_norm)) then
      allocate (ppcr%p_norm(npc,nr), ppcr%integ_pc_out(npc), ppcr%integ_r(npc,nr), ppcr%integ_r_ave(nr))
    endif

    ppcr%pc_out_max = ppcr%pc_out(npc)
    if (ele%value(pc_out_max$) /= 0) ppcr%pc_out_max = min(ppcr%pc_out_max, ele%value(pc_out_max$))
    ppcr%pc_out_min = max(ppcr%pc_out(1), ele%value(pc_out_min$))
    if (ppcr%pc_out_max <= ppcr%pc_out_min) then
      call out_io(s_fatal$, r_name, 'PC_OUT_MAX IS LESS THAN OR EQUAL TO PC_OUT_MIN. FOR ELEMENT: ' // ele%name, &
                                    'PARTICLE WILL BE MARKED AS LOST.')
      return
    endif

    do ipc = 1, npc
      do ir = 1, nr
        out%pc_out = ppcr%pc_out(ipc)
        out%r = ppcr%r(ir)
        call calc_dir_out_params (ele, dist, sub_dist, .false., out, err_flag);  if (err_flag) return
        unnorm_integ_prob_tot = out%integ_prob_tot   ! integ_prob_tot is not normalized by angle_max.
        call calc_dir_out_params (ele, dist, sub_dist, .true., out, err_flag);  if (err_flag) return
        ppcr%p_norm(ipc,ir) = ppcr%prob(ipc,ir) *  unnorm_integ_prob_tot / out%integ_prob_tot  ! p_norm is normalized by angle_max
      enddo
    enddo

    ppcr%integ_pc_out(1) = 0
    do ipc = 2, npc
      ppcr%integ_pc_out(ipc) = ppcr%integ_pc_out(ipc-1) + super_qromb_2D(p2_norm_func, ppcr%pc_out(ipc-1), ppcr%pc_out(ipc), &
                                                                               0.0_rp, ppcr%r(nr), 1e-4_rp, 0.0_rp, 2, err_flag)
    enddo
    ppcr%integrated_prob = ppcr%integ_pc_out(npc)
    ppcr%integ_pc_out = ppcr%integ_pc_out / ppcr%integrated_prob

    ppcr%integ_r(:,1) = 0
    do ir = 2, nr
      do ipc = 1, npc
        com%ipc = ipc
        ppcr%integ_r(ipc,ir) = ppcr%integ_r(ipc,ir-1) + super_qromb(p1_norm_func, ppcr%r(ir-1), ppcr%r(ir), &
                                                                                       1e-4_rp, 0.0_rp, 2, err_flag)
      enddo
    enddo

    do ipc = 1, npc
      ppcr%integ_r(ipc,:) = ppcr%integ_r(ipc,:) / ppcr%integ_r(ipc,nr)
    enddo

  enddo
enddo

end subroutine probability_tables_setup

!------------------------------------------------
! contains

function p2_norm_func (x,y) result (value)

real(rp) x, y, value
real(rp) dx, dy
integer ix, iy, nx, ny

!

nx = size(com%ppcr%pc_out)
ny = size(com%ppcr%r)
ix = bracket_index(x, com%ppcr%pc_out, 1, dx, restrict = .true.)
iy = bracket_index(y, com%ppcr%r, 1, dy, restrict = .true.)

if (x < com%ppcr%pc_out_min .or. x > com%ppcr%pc_out_max) then
  value = 0
  return
endif

value = (1-dx)*(1-dy)*com%ppcr%p_norm(ix,iy) + dx*(1-dy)*com%ppcr%p_norm(ix+1,iy) + &
        (1-dx)*dy*com%ppcr%p_norm(ix,iy+1) + dx*dy*com%ppcr%p_norm(ix+1,iy+1) 

end function p2_norm_func

!------------------------------------------------
! contains

function p1_norm_func (r) result (value)

real(rp), intent(in) :: r(:)
real(rp) :: value(size(r)), dr
integer i, ix

!

do i = 1, size(r)
  ix = bracket_index(r(i), com%ppcr%r, 1, dr, restrict = .true.)
  value(i) = (1-dr)*com%ppcr%p_norm(com%ipc, ix) + dr*com%ppcr%p_norm(com%ipc, ix+1)
enddo


end function p1_norm_func

!------------------------------------------------
! contains

subroutine calc_dir_out_params (ele, dist, sub_dist, restrict_angle, out, err_flag)

type (ele_struct) ele
type (converter_distribution_struct), target :: dist
type (converter_sub_distribution_struct), target :: sub_dist
type (converter_param_storage_struct) out
type (spline_struct), pointer :: spn(:)

real(rp) dr, dx, x_min, x, rad, a_tan, b1, drad, arg, dxy_ds_max
real(rp), pointer :: integ(:)

integer i, n, ix
logical restrict_angle, err_flag

! 

err_flag = .true.

if (.not. dir_out_calc (out%beta, out, sub_dist%dir_out%beta)) return
if (.not. dir_out_calc (out%c_x, out, sub_dist%dir_out%c_x)) return
if (.not. dir_out_calc (out%alpha_x, out, sub_dist%dir_out%alpha_x)) return
if (.not. dir_out_calc (out%alpha_y, out, sub_dist%dir_out%alpha_y)) return
if (.not. dir_out_calc (out%dxds_min, out, sub_dist%dir_out%dxds_min)) return
if (.not. dir_out_calc (out%dxds_max, out, sub_dist%dir_out%dxds_max)) return
if (.not. dir_out_calc (out%dyds_max, out, sub_dist%dir_out%dyds_max)) return

dxy_ds_max = atan(ele%value(angle_out_max$))
if (restrict_angle .and. dxy_ds_max > 0) then
  out%dxds_max = min(out%dxds_max, dxy_ds_max)
  out%dxds_min = min(max(out%dxds_min, -dxy_ds_max), out%dxds_max)
  out%dyds_max = min(out%dyds_max, dxy_ds_max)
endif

if (1 + out%beta * out%dxds_min < 0 .or. 1 + out%beta * out%dxds_max < 0) then
  call out_io (s_error$, r_name, 'BETA VALUE TOO LARGE: \f12.4\ ', &
                                 '  AT (PC_OUT, R) = (\2es12.4\)', &
                                 '  PARTICLE WILL BE MARKED AS LOST.', r_array = [out%beta, out%pc_out, out%r])
  return
endif

if (out%alpha_x < 0) then
  call out_io (s_error$, r_name, 'ALPHA_X VALUE IS NEGATIVE: \f12.4\ ', &
                                 '  AT (PC_OUT, R) = (\2es12.4\)', &
                                 '  PARTICLE WILL BE MARKED AS LOST.', r_array = [out%alpha_x, out%pc_out, out%r])
  return
endif

if (out%alpha_y < 0) then
  call out_io (s_error$, r_name, 'ALPHA_Y VALUE IS NEGATIVE: \f12.4\ ', &
                                 '  AT (PC_OUT, R) = (\2es12.4\)', &
                                 '  PARTICLE WILL BE MARKED AS LOST.', r_array = [out%alpha_y, out%pc_out, out%r])
  return
endif

! For dx/ds use a spine fit.

spn => com%dxds_spline
integ => com%dxds_integ
integ(0) = 0
dx = (out%dxds_max - out%dxds_min) / n_pt$
x_min = out%dxds_min

do i = 0, n_pt$
  x = x_min + i * dx
  rad = 1.0_rp / sqrt(1 + (out%alpha_x * (x - out%c_x))**2)
  drad = -out%alpha_x**2 * (x - out%c_x) * rad**3
  arg = out%alpha_y * out%dyds_max * rad
  a_tan = atan(arg)
  b1 = 1 + out%beta * x

  spn(i)%x0 = x
  spn(i)%y0 = b1 * rad * a_tan
  spn(i)%coef(1) = out%beta * rad * a_tan + b1 * drad * (rad * out%alpha_y * out%dyds_max /(1 + arg**2) + a_tan)

  if (i > 0) then
    spn(i-1) = create_a_spline([spn(i-1)%x0, spn(i-1)%y0], [spn(i)%x0, spn(i)%y0], spn(i-1)%coef(1), spn(i)%coef(1))
    integ(i) = integ(i-1) + spline1(spn(i-1), x, -1)
  endif
enddo

out%integ_prob_tot = integ(n_pt$)
err_flag = .false.

end subroutine calc_dir_out_params

!------------------------------------------------
! contains

function dir_out_calc (value, out, dc) result (is_ok)

type (converter_dir_coef_struct) dc
type (converter_param_storage_struct) out
real(rp) value, dr
integer n, ix
logical is_ok

!

n = size(dc%fit_1d_r)

if (out%pc_out >= dc%fit_1d_r(n)%pc_out) then
  value = poly_eval(dc%fit_2d_pc%poly, out%pc_out) * poly_eval(dc%fit_2d_r%poly, out%r) * &
                 exp(-dc%fit_2d_pc%k * out%pc_out - dc%fit_2d_r%k * out%r) + dc%c0
else
  ix = bracket_index(out%pc_out, dc%fit_1d_r%pc_out, 1, dr, restrict = .true.)
  value = (1 - dr) * poly_eval(dc%fit_1d_r(ix)%poly, out%r) + dr * poly_eval(dc%fit_1d_r(ix+1)%poly, out%r) 
endif

is_ok = .true.

end function dir_out_calc

end subroutine track_a_converter
