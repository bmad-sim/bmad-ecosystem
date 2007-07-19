!+
! Subroutine slice_ele_calc (ele, param, i_slice, n_slice_tot, sliced_ele)
!
! Routine to create an element that represents a slice of another element.
! This routine can be used for detailed tracking through an element.
!
! i_slice should range from 1 to n_slice_tot.
! The longitudinal position from the beginning of the element of 
! the entrance face of the i^th slice is at:
!     s_start = length_ele * (i_slice-1) / n_slice_tot 
! The exit face of the i^th slice is at: 
!     s_end = length_ele * i_slice / n_slice_tot 
!
! It is assumed that this routine will be called in a loop so the set:
!     sliced_ele = ele
! will only be done when i_slice = 1.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele         -- Ele_struct: Element to slice and dice.
!   param       -- lat_param_struct: Lattice parameters
!   i_slice     -- Integer: Slice index
!   n_slice_tot -- Integer: Total number of slices.
!
! Output:
!   slice_ele -- Ele_struct: Sliced portion of ele.
!-

subroutine slice_ele_calc (ele, param, i_slice, n_slice_tot, sliced_ele)

use bookkeeper_mod, except_dummy => slice_ele_calc

implicit none

type (ele_struct) ele, sliced_ele
type (lat_param_struct) param

integer i_slice, n_slice_tot
real(rp) phase, s_start, E

logical bookit

!

bookit = .false.

if (i_slice == 1) then
  sliced_ele = ele
  sliced_ele%value(l$) = ele%value(l$) / n_slice_tot
  sliced_ele%value(hkick$) = ele%value(hkick$) / n_slice_tot
  sliced_ele%value(vkick$) = ele%value(vkick$) / n_slice_tot
  if (ele%key == hkicker$ .or. ele%key == vkicker$) &
                        sliced_ele%value(kick$) = ele%value(kick$) / n_slice_tot
  bookit = .true.
endif

if (ele%key == sbend$) then
  if (i_slice == 1) then
    sliced_ele%value(e1$) = ele%value(e1$)
    sliced_ele%value(e2$) = 0
  elseif (i_slice == n_slice_tot) then
    sliced_ele%value(e1$) = 0
    sliced_ele%value(e2$) = ele%value(e2$)
  else 
    sliced_ele%value(e1$) = 0
    sliced_ele%value(e2$) = 0
  endif
endif

if (ele%key == lcavity$) then
  sliced_ele%value(e_loss$) = sliced_ele%value(e_loss$) / n_slice_tot
  phase = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
  s_start = (i_slice - 1) * ele%value(l$) / n_slice_tot
  E = ele%value(E_TOT$) + ele%value(gradient$) * s_start * cos(phase)
  E = E - sliced_ele%value(e_loss$) * param%n_part * e_charge
  sliced_ele%value(E_TOT_START$) = E 
  call convert_total_energy_to (E, param%particle, pc = sliced_ele%value(p0c_start$))
  bookit = .true.
endif

if (bookit .and. .not. bmad_com%auto_bookkeeper) &
                              call attribute_bookkeeper (sliced_ele, param)

end subroutine
