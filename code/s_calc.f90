!+
! Subroutine s_calc (lat)
!
! Subroutine to calculate the longitudinal distance S for the elements
! in a lattice.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat -- lat_struct:
!
! Output:
!   lat -- lat_struct:
!-

subroutine s_calc (lat)

use bmad_struct
use bmad_interface, except_dummy => s_calc
use lat_ele_loc_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, ele0, ele1
type (branch_struct), pointer :: branch

integer i, j, n, ic, icon, ix2, nt
real(8) ss, s_end

! Just go through all the elements and add up the lengths.
! The last super_slave (the super_slave at the exit end) of a super_lord
! has it's length adjusted to be compatable with the length of the super_lord

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (i > 0) branch%ele(0)%s = 0  ! Branches start from zero
  ss = branch%ele(0)%s
  nt = branch%n_ele_track
  do n = 1, nt
    ele => branch%ele(n)
    if (ele%slave_status == super_slave$) then
      do ic = ele%ic1_lord, ele%ic2_lord
        icon = lat%ic(ic)
        lord => lat%ele(lat%control(icon)%ix_lord)
        call find_element_ends (lat, lord, ele0, ele1)
        if (ele1%ix_ele == n) then ! Is last super_slave
          s_end = ele0%s + lord%value(l$)
          ! If the super_lord wraps around the lattice ends then must adjust s_end
          if (s_end > branch%ele(n)%s) s_end = s_end - (branch%ele(0)%s - branch%ele(n)%s) 
          ele%value(l$) = s_end - branch%ele(n-1)%s
        endif
      enddo
    endif
    ss = ss + ele%value(l$)
    lat%ele(n)%s = ss
  enddo
enddo

lat%param%total_length = ss - lat%ele(0)%s

! Now fill in the s positions of the super_lords and zero everyone else.
! Exception: A null_ele lord element is the result of a superposition on a multipass section.
! We need to preserve the s value of this element.

do n = lat%n_ele_track+1, lat%n_ele_max
  if (lat%ele(n)%key == null_ele$) cycle
  if (lat%ele(n)%lord_status == super_lord$) then
    ix2 = lat%control(lat%ele(n)%ix2_slave)%ix_slave
    lat%ele(n)%s = lat%ele(ix2)%s
  else
    lat%ele(n)%s = 0
  endif
enddo

end subroutine
