!+
! Program: element_attributes
!
! Program to print which element attributes are defined for different types of elements.
! This is useful when new attributes need to be defined for a given element type.
! Generally only Bmad wizards are interested in this.
!-

program element_attributes

use basic_attribute_mod

implicit none

type (ele_struct) ele
integer i, j, n_used(n_attrib_special_maxx)
character(40) a_name

! indexes used for a given element key

open (1, file = 'element_attributes.dat')

n_used = 0

do i = 1, n_key
  if (i == overlay_lord$) cycle
  write (1, *) '!---------------------------------'
  write (1, *) key_name(i)
  ele%key = i
  do j = 1, n_attrib_special_maxx
    a_name = attribute_name (ele, j) 
    if (a_name == null_name$) cycle
    write (1, '(i10, 2x, a)') j, a_name
    n_used(j) = n_used(j) + 1
  enddo
  write (1, *)
enddo

! number of elements using an index

write (1, *) '!---------------------------------'
write (1, *) 'Index usage:'
write (1, *) '   Ix Count'
do i = 1, size(n_used)
  write (1, '(2i6)') i, n_used(i)
enddo

! List of elements using an index

do i = 1, n_attrib_special_maxx
  write (1, *)
  write (1, *) '!---------------------------------'
  write (1, '(a, i0)') 'Index: ', i
  do j = 1, n_key
    if (j == overlay_lord$) cycle
    ele%key = j
    a_name = attribute_name (ele, i) 
    if (a_name == null_name$) cycle
    write (1, *) '   ', key_name(ele%key)
  enddo
enddo

print *, 'Writen: element_attributes.dat'

end program
