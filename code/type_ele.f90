!+
! Subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, 
!                              twiss_out, type_control, lattice, type_wake)
!
! Subroutine to type out information on an element. 
! See also the subroutine type2_ele.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele              -- Ele_struct: Element
!   type_zero_attrib -- Logical, optional: If False then surpress printing of
!                          attributes whose value is 0. Default is False.
!   type_mat6      -- Integer, optional:
!                         = 0   => Do not type ele%mat6
!                         = 4   => Type 4X4 xy submatrix
!                         = 6   => Type full 6x6 matrix (Default)
!   type_taylor    -- Logical, optional: Print out taylor map terms?
!                       If ele%taylor is not allocated then this is ignored.
!                       Default is False.
!   twiss_out      -- Integer, optional: Print the Twiss parameters at the 
!                         element end?
!                       = 0         => Do not print the Twiss parameters
!                       = radians$  => Print Twiss, phi in radians (Default).
!                       = degrees$  => Print Twiss, phi in degrees.
!                       = cycles$   => Print Twiss, phi in radians/2pi.
!   type_control   -- Logical, optional: If True then print control status.
!                       Default is False if lattice is not present. Otherwise True.
!   lattice        -- lat_struct, optional: Needed for control typeout.
!   type_wake      -- Logical, optional: If True then print the long-range and 
!                       short-range wakes information. If False then just print
!                       how many terms the wake has. Default is True.
!                       If ele%wake is not allocated then this is ignored.
!-

#include "CESR_platform.inc"

subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor,  &
                                  twiss_out, type_control, lattice, type_wake)

  use bmad_struct
  use bmad_interface, except => type_ele

  implicit none

  type (ele_struct)  ele
  type (lat_struct), optional :: lattice

  integer n_lines, i
  integer, optional :: type_mat6, twiss_out

  logical, optional :: type_control, type_zero_attrib, type_taylor, type_wake

  character(100), pointer :: lines(:) 

!

  nullify (lines)

  call type2_ele (ele, lines, n_lines, type_zero_attrib, type_mat6, &
                         type_taylor, twiss_out, type_control, lattice, type_wake)

  do i = 1, n_lines
    print '(1x, a)', trim(lines(i))
  enddo

  deallocate (lines)

end subroutine
