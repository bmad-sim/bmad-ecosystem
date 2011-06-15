!+
! Program: element_attributes
!
! Program to print which element attributes are defined for different types of elements.
! This is useful when new attributes need to be defined for a given element type.
! Generally only Bmad wizards are interested in this.
!-

program element_attributes

use bmad

implicit none

type (ele_struct) ele
integer i, j
character(40) a_name

!

do i = 1, n_key
  print *
  print *, '!---------------------------------'
  print *, key_name(i)
  ele%key = i
  do j = 1, n_attrib_special_maxx
    a_name = attribute_name (ele, j) 
    if (a_name(1:1) == '!') cycle
    print '(i10, 2x, a)', j, a_name
  enddo
enddo

end program
