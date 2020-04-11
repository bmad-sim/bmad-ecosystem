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
type this_common_struct
  type (ele_struct), pointer :: ele
  real(rp) E_out, r
  real(rp) A_dir
end type

type (this_common_struct) com
type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (converter_struct), pointer :: conv

real(rp), optional :: mat6(6,6)

integer i

logical, optional :: make_matrix
character(*), parameter :: r_name = 'track_a_converter'

!

conv => ele%converter
if (size(conv%dist) == 0) then
  call out_io (s_error$, r_name, 'CONVERTER ELEMENT DOES NOT HAVE ANY ASSOCIATED DISTRIBUTION(S) FOR: ' // ele%name, &
                                 'PARTICLE WILL BE MARKED AS LOST.')
  orbit%state = lost$
  return
endif

if (.not. allocated(ele%converter%dist(1)%sub_dist(1)%prob_E_r%p_norm)) call calc_integ_prob(ele)

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

! Calculate E_out




call offset_particle (ele, param, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

!------------------------------------------------
contains

subroutine calc_integ_prob(ele)
type (ele_struct), target :: ele
type (converter_distribution_struct), pointer :: dist
type (converter_direction_out_struct), pointer :: dir_out
type (converter_prob_E_r_struct), pointer :: per
real(rp) dxy_ds_max, angle_max, prob
real(rp) E_out, r, A_dir
integer id, isd, ie, ir, ne, nr
logical err_flag

!

do id = 1, size(ele%converter%dist)
  dist =>ele%converter%dist(id)
  dxy_ds_max = dist%dxy_ds_max
  if (ele%value(angle_out_max$) /= 0) angle_max = min(ele%value(angle_out_max$), angle_max)
  do isd = 1, size(dist%sub_dist)
    dir_out => dist%sub_dist(isd)%dir_out
    per => dist%sub_dist(isd)%prob_E_r
    ne = size(per%E_out)
    nr = size(per%r)
    if (.not. allocated(per%p_norm)) then
      allocate (per%p_norm(ne,nr), per%p_integ_E_out(ne), per%p_integ_r(ne,nr))
      !!! allocate
    endif

    do ie = 1, size(per%E_out)
      E_out = per%E_out(ie)
      do ir = 1, size(per%r)
        r = per%r(ir)
        com = this_common_struct(ele, E_out, r, 1.0_rp)
        prob = super_qromb_2D(p_angle, 0.0_rp, dxy_ds_max, 0.0_rp, dxy_ds_max, 1e-5_rp, 0.0_rp, 5, err_flag)
        A_dir = 1 / prob   ! A_dir is not normalized by angle_max.
        com = this_common_struct(ele, E_out, r, A_dir)
        prob = super_qromb_2D(p_angle, 0.0_rp, angle_max, 0.0_rp, angle_max, 1e-5_rp, 0.0_rp, 5, err_flag)
        per%p_norm(ie,ir) = per%prob(ie,ir) * prob  ! p_norm is now normalized by angle_max
      enddo
    enddo
  enddo
enddo

end subroutine calc_integ_prob

!------------------------------------------------
! contains

function p_angle (x, y) result (value)

real(rp) x, y, value
real(rp) beta, alpha_x, alpha_y, c_x

value = 1

end function p_angle

end subroutine
