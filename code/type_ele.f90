!+
! Subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, twiss_out, 
!        type_control, type_wake, type_floor_coords, type_field, 
!        type_wall, nunit)
!
! Subroutine to type out information on an element. 
! See also the subroutine type2_ele.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele                -- Ele_struct: Element
!   type_zero_attrib   -- Logical, optional: If False then surpress printing of
!                            attributes whose value is 0. Default is False.
!   type_mat6          -- Integer, optional:
!                             = 0   => Do not type ele%mat6
!                             = 4   => Type 4X4 xy submatrix
!                             = 6   => Type full 6x6 matrix (Default)
!   type_taylor        -- Logical, optional: Print out taylor map terms?
!                           If ele%taylor is not allocated then this is ignored.
!                           Default is False.
!   twiss_out          -- Integer, optional: Print the Twiss parameters at the 
!                             element end?
!                           = 0         => Do not print the Twiss parameters
!                           = radians$  => Print Twiss, phi in radians (Default).
!                           = degrees$  => Print Twiss, phi in degrees.
!                           = cycles$   => Print Twiss, phi in radians/2pi.
!   type_control       -- Logical, optional: Print control status? Default is True.
!                           If ele%branch%lat is not associated cannot print status info.
!   lattice            -- lat_struct, optional: Needed for control typeout.
!   type_wake          -- Logical, optional: If True then print the long-range and 
!                           short-range wakes information. If False then just print
!                           how many terms the wake has. Default is True.
!                           If ele%wake is not allocated then this is ignored.
!   type_floor_coords  -- Logical, optional: If True then print the global ("floor")
!                           coordinates at the exit end of the element.
!                           Default is False.
!   type_field         -- Logical, optional: If True then print:
!                           Wiggler terms for a a map_type wiggler or
!                           RF field coefficients for a lcavity or rfcavity.
!                           Default is False.
!   type_wall          -- Logical, optional: If True then print wall info
!                           for a capillary. Default is False.
!   nunit              -- Integer, optional: Unit for writing:
!                             > 0 output to file only with unit = nunit
!                             = 0 output to terminal only (default)
!                             < 0 output to terminal and file with unit = abs(nunit).
!-

subroutine type_ele (ele, type_zero_attrib, type_mat6, type_taylor, twiss_out, &
      type_control, type_wake, type_floor_coords, type_field, &
      type_wall, nunit)

use bmad_interface, except_dummy => type_ele

implicit none

type (ele_struct)  ele

integer n_lines, i, iu
integer, optional :: type_mat6, twiss_out, nunit

logical, optional :: type_control, type_zero_attrib, type_taylor, type_wake
logical, optional :: type_floor_coords, type_field, type_wall

character(200), allocatable :: lines(:) 

!

call type2_ele (ele, lines, n_lines, type_zero_attrib, type_mat6, type_taylor, &
      twiss_out, type_control, type_wake, type_floor_coords, type_field, &
      type_wall)

iu = integer_option(0, nunit)

if (iu <= 0) then
  do i = 1, n_lines
    print '(1x, a)', trim(lines(i))
  enddo
endif

if (iu /= 0) then
  do i = 1, n_lines
    write (abs(iu), '(1x, a)') trim(lines(i))
  enddo
endif

end subroutine
