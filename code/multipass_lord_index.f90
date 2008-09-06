!+
! Function multipass_lord_index (ix_ele, lat, ix_pass, ix_super_lord) result (ix_multi_lord)
!
! Routine to find the index of the multipass lord of a lattice element.
! ix_multi_lord will be positive for elements:
!   multipass_slaves
!   super_lords that are slaves of a multipass_lord
!   super_slaves whose super_lord is a slave of a multipass_lord
!
! Modules needed:
!   use bmad
!
! Input:
!   ix_ele   -- Integer: Index of the lattice element.
!   lat      -- Lat_struct: Lattice containing the element.
!
! Output:
!   ix_pass       -- Integer, optional: Multipass turn number.
!   ix_super_lord -- Integer, optional: Index of the super_lord of the element.
!                      Set to -1 if the element is not a super_slave.
!   ix_multi_lord -- Integer: Index of the multipass_lord if there is one.
!                      Set to -1 if there is no multipass_lord.
!-

function multipass_lord_index (ix_ele, lat, ix_pass, ix_super_lord) result (ix_multi_lord)

use bmad_struct
use bmad_interface

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

integer ix_ele, ix_multi_lord
integer, optional :: ix_super_lord, ix_pass
integer ix_sup, ix_mult, ic

!

ele => lat%ele(ix_ele)

ix_multi_lord = -1
if (present(ix_super_lord)) ix_super_lord = -1

if (ele%control_type == multipass_slave$) then
  ic = lat%ic(ele%ic1_lord)
  ix_multi_lord = lat%control(ic)%ix_lord
  if (present(ix_pass)) ix_pass = ic + 1 - lat%ele(ix_multi_lord)%ix1_slave
  return
endif

if (ele%control_type == super_slave$) then
  ic = lat%ic(ele%ic1_lord)
  ix_sup = lat%control(ic)%ix_lord
  if (lat%ele(ix_sup)%n_lord == 0) return
  ic = lat%ic(lat%ele(ix_sup)%ic1_lord)
  ix_mult = lat%control(ic)%ix_lord
  if (lat%ele(ix_mult)%control_type /= multipass_lord$) return
  ix_multi_lord = ix_mult
  if (present(ix_super_lord)) ix_super_lord = ix_sup
  if (present(ix_pass)) ix_pass = ic + 1 - lat%ele(ix_multi_lord)%ix1_slave
  return
endif

if (ele%control_type == super_lord$) then
  if (ele%n_lord == 0) return
  ic = lat%ic(ele%ic1_lord)
  ix_mult = lat%control(ic)%ix_lord
  if (lat%ele(ix_mult)%control_type /= multipass_lord$) return
  ix_multi_lord = ix_mult
  if (present(ix_pass)) ix_pass = ic + 1 - lat%ele(ix_multi_lord)%ix1_slave
endif

end function
