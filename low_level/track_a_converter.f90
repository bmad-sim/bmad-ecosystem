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
use spline_mod

implicit none

integer, parameter :: n_pt = 100

type converter_param_storage_struct
  real(rp) pc_out, r
  real(rp) dx_ds, dy_ds
  real(rp) weight
  logical :: lost = .false.
  real(rp) beta, alpha_x, alpha_y, c_x, integ_prob_tot
end type

type converter_common_struct
  type (converter_prob_pc_r_struct), pointer :: ppcr
  type (converter_distribution_struct), pointer :: dist
  type (spline_struct) dx_ds_spline(0:n_pt)
  real(rp) dx_ds_integ(0:n_pt)
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
real(rp) r_ran(5), pc_in, thickness, drd, drd2, azimuth_angle

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

call offset_particle (ele, param, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

! Find dists that straddle thickness

if (.not. allocated(ele%converter%dist(1)%sub_dist(1)%prob_pc_r%p_norm)) call probability_tables_setup(ele, err_flag)
if (err_flag) goto 9000

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
  out1%dx_ds  = drd2 * out1%dx_ds  + drd * out2%dx_ds
  out1%dy_ds  = drd2 * out1%dy_ds  + drd * out2%dy_ds
  out1%weight = drd2 * out1%weight + drd * out2%weight
endif

if (out1%lost) then
  orbit%state = lost$
  return
endif

azimuth_angle = twopi * r_ran(5)

orb0 = orbit
orbit%species = conv%species_out
orbit%p0c = ele%value(p0c$)
orbit%vec(6) = out1%pc_out / ele%value(p0c$) - 1
call convert_pc_to (out1%pc_out, orbit%species, beta = orbit%beta)

orbit%charge = out1%weight
orbit%vec(1) = orbit%vec(1) + out1%r * cos(azimuth_angle)
orbit%vec(2) = orbit%vec(2) + (out1%dx_ds * cos(azimuth_angle) - out1%dy_ds * sin(azimuth_angle)) * (1 + orbit%vec(6))
orbit%vec(3) = orbit%vec(3) + out1%r * sin(azimuth_angle)
orbit%vec(4) = orbit%vec(4) + (out1%dx_ds * sin(azimuth_angle) + out1%dy_ds * cos(azimuth_angle)) * (1 + orbit%vec(6))
orbit%vec(5) = orbit%vec(5) * orb0%beta / orbit%beta

call offset_particle (ele, param, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)
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
  call out_io (s_fatal$, r_name, 'INCOMING PARTICLE ENERGY OUT OF RANGE: \es12.4\ FOR CONVERTER: ' // ele%name, &
                                 'PARTICLE WILL BE MARKED AS LOST.')
  out%lost = .true.
  return
endif

ix_sd = bracket_index(pc_in, dist%sub_dist%pc_in, 1, rsd, restrict = .true.)
call calc_out_coords2 (ele, dist, dist%sub_dist(ix_sd), r_ran, out1, err_flag);  if (err_flag) return
call calc_out_coords2 (ele, dist, dist%sub_dist(ix_sd+1), r_ran, out2, err_flag);  if (err_flag) return

rsd2 = 1 - rsd
out%pc_out = rsd2 * out1%pc_out + rsd * out2%pc_out
out%r      = rsd2 * out1%r      + rsd * out2%r
out%dx_ds  = rsd2 * out1%dx_ds  + rsd * out2%dx_ds
out%dy_ds  = rsd2 * out1%dy_ds  + rsd * out2%dy_ds
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

! dx/ds calc

call calc_dir_out_params(dist, sub_dist, out, err_flag);  if (err_flag) return

com%r_ran = r_ran(3)
com%integ_prob_tot = out%integ_prob_tot
out%dx_ds = super_zbrent(dx_ds_func, -dist%dxy_ds_limit, dist%dxy_ds_limit, 0.0_rp, 1e-4_rp, status)

! dy/ds calc

b = 1 + out%alpha_x * (out%dx_ds - out%c_x)**2
k_const = sqrt(out%alpha_y * b) / (2 * atan(sqrt(out%alpha_y/b) * dist%dxy_ds_limit))
out%dy_ds = sqrt(b/out%alpha_y) * tan(sqrt(out%alpha_y * b) * (r_ran(4) - 0.5_rp) / k_const)

end subroutine calc_out_coords2

!------------------------------------------------
! contains

function dx_ds_func (x, status) result (value)

real(rp), intent(in) :: x
real(rp) value
integer status, ix

!

ix = bracket_index(x, com%dx_ds_spline%x0, 0, restrict = .true.)
value = (com%dx_ds_integ(ix) + spline1(com%dx_ds_spline(ix), x, -1)) - com%integ_prob_tot * com%r_ran

end function dx_ds_func

!------------------------------------------------
! contains

subroutine probability_tables_setup(ele, err_flag)

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
  dist%dxy_ds_limit = dist%dxy_ds_max
  if (ele%value(angle_out_max$) > 0) dist%dxy_ds_limit = min(atan(ele%value(angle_out_max$)), dist%dxy_ds_limit)
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
      err_flag = .true.
    endif

    do ipc = 1, npc
      do ir = 1, nr
        out%pc_out = ppcr%pc_out(ipc)
        out%r = ppcr%r(ir)
        dist%dxy_ds_limit = dist%dxy_ds_max
        call calc_dir_out_params (dist, sub_dist, out, err_flag);  if (err_flag) return
        unnorm_integ_prob_tot = out%integ_prob_tot   ! integ_prob_tot is not normalized by angle_max.
        if (ele%value(angle_out_max$) > 0) dist%dxy_ds_limit = min(atan(ele%value(angle_out_max$)), dist%dxy_ds_limit)
        call calc_dir_out_params (dist, sub_dist, out, err_flag);  if (err_flag) return
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

if (x < com%ppcr%pc_out_min .or. x > com%ppcr%pc_out_max .or. y > com%dist%dxy_ds_limit) then
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
  if (r(i) > com%dist%dxy_ds_limit) then
    value(i) = 0
    cycle
  endif
  ix = bracket_index(r(i), com%ppcr%r, 1, dr, restrict = .true.)
  value(i) = (1-dr)*com%ppcr%p_norm(com%ipc, ix) + dr*com%ppcr%p_norm(com%ipc, ix+1)
enddo


end function p1_norm_func

!------------------------------------------------
! contains

subroutine calc_dir_out_params (dist, sub_dist, out, err_flag)

type (converter_distribution_struct), target :: dist
type (converter_sub_distribution_struct), target :: sub_dist
type (converter_param_storage_struct) out
type (converter_beta_struct), pointer :: beta
type (converter_alpha_struct), pointer :: alpha
type (converter_c_x_struct), pointer :: c_x
type (spline_struct), pointer :: spn(:)

real(rp) dr, dx, x_min, x, rad, a_tan, b1, drad, arg
real(rp), pointer :: integ(:)

integer i, n, ix
logical err_flag

! beta calc

err_flag = .true.

beta => sub_dist%dir_out%beta
n = size(beta%fit_1d_r)
if (out%pc_out >= beta%fit_1d_r(n)%pc_out) then
  out%beta = poly_eval(beta%poly_pc, out%pc_out) * poly_eval(beta%poly_r, out%r)
else
  ix = bracket_index(out%pc_out, beta%fit_1d_r%pc_out, 1, dr, restrict = .true.)
  out%beta = (1 - dr) * poly_eval(beta%fit_1d_r(ix)%poly, out%r) + &
                  dr * poly_eval(beta%fit_1d_r(ix+1)%poly, out%r) 
endif

if (abs(out%beta) * dist%dxy_ds_limit > 1) then
  call out_io (s_error$, r_name, 'BETA VALUE TOO LARGE: \f12.4\ ', &
                                 '  AT (PC_OUT, R) = (\2es12.4\)', &
                                 '  PARTICLE WILL BE MARKED AS LOST.', r_array = [out%beta, out%pc_out, out%r])
  return
endif

! c_x calc

c_x => sub_dist%dir_out%c_x
out%c_x = poly_eval(c_x%poly_pc, out%pc_out) * poly_eval(c_x%poly_r, out%r)

! Alpha_x calc

alpha => sub_dist%dir_out%alpha_x
n = size(alpha%fit_1d_r)
if (out%pc_out >= alpha%fit_1d_r(n)%pc_out) then
  out%alpha_x = poly_eval(alpha%fit_2d_pc%poly, out%pc_out) * poly_eval(alpha%fit_2d_r%poly, out%r) * &
                exp(-(alpha%fit_2d_pc%k * out%pc_out + alpha%fit_2d_r%k * out%r))
else
  ix = bracket_index(out%pc_out, alpha%fit_1d_r%pc_out, 1, dr, restrict = .true.)
  out%alpha_x = (1 - dr) * poly_eval(alpha%fit_1d_r(ix)%poly, out%r) * exp(-alpha%fit_1d_r(ix)%k * out%r) + &
                   dr * poly_eval(alpha%fit_1d_r(ix+1)%poly, out%r) * exp(-alpha%fit_1d_r(ix+1)%k * out%r)
endif

if (out%alpha_x < 0) then
  call out_io (s_error$, r_name, 'ALPHA_X VALUE IS NEGATIVE: \f12.4\ ', &
                                 '  AT (PC_OUT, R) = (\2es12.4\)', &
                                 '  PARTICLE WILL BE MARKED AS LOST.', r_array = [out%alpha_x, out%pc_out, out%r])
  return
endif

! Alpha_y calc

alpha => sub_dist%dir_out%alpha_y
n = size(alpha%fit_1d_r)
if (out%pc_out >= alpha%fit_1d_r(n)%pc_out) then
  out%alpha_y = poly_eval(alpha%fit_2d_pc%poly, out%pc_out) * poly_eval(alpha%fit_2d_r%poly, out%r) * &
                exp(-(alpha%fit_2d_pc%k * out%pc_out + alpha%fit_2d_r%k * out%r))
else
  ix = bracket_index(out%pc_out, alpha%fit_1d_r%pc_out, 1, dr, restrict = .true.)
  out%alpha_y = (1 - dr) * poly_eval(alpha%fit_1d_r(ix)%poly, out%r) * exp(-alpha%fit_1d_r(ix)%k * out%r) + &
                   dr * poly_eval(alpha%fit_1d_r(ix+1)%poly, out%r) * exp(-alpha%fit_1d_r(ix+1)%k * out%r)
endif

if (out%alpha_y < 0) then
  call out_io (s_error$, r_name, 'ALPHA_Y VALUE IS NEGATIVE: \f12.4\ ', &
                                 '  AT (PC_OUT, R) = (\2es12.4\)', &
                                 '  PARTICLE WILL BE MARKED AS LOST.', r_array = [out%alpha_y, out%pc_out, out%r])
  return
endif

! For dx/ds use a spine fit.

spn => com%dx_ds_spline
integ => com%dx_ds_integ
integ(0) = 0
dx = 2 * dist%dxy_ds_limit / n_pt
x_min = -dist%dxy_ds_limit

do i = 0, n_pt
  x = x_min + i * dx
  rad = 1.0_rp / sqrt(out%alpha_y * (1 + out%alpha_x * (x - out%c_x)**2))
  drad = -out%alpha_x * out%alpha_y * (x - out%c_x)
  arg = out%alpha_y * rad * dist%dxy_ds_limit
  a_tan = atan(arg)
  b1 = 1 + out%beta * x

  spn(i)%x0 = x
  spn(i)%y0 = b1 * rad * a_tan
  spn(i)%coef(1) = out%beta * rad * a_tan + b1 * drad * &
                        (rad * out%alpha_y * dist%dxy_ds_limit /(1 + arg**2) + a_tan)

  if (i > 0) then
    spn(i-1) = create_a_spline([spn(i-1)%x0, spn(i-1)%y0], [spn(i)%x0, spn(i)%y0], spn(i-1)%coef(1), spn(i)%coef(1))
    integ(i) = integ(i-1) + spline1(spn(i-1), x, -1)
  endif
enddo

out%integ_prob_tot = 1 / integ(n_pt)
err_flag = .false.

end subroutine calc_dir_out_params

end subroutine track_a_converter
