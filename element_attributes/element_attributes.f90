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
integer i, j, n_used(n_attrib_special_maxx)
character(40) a_name

! indexes used for a given element key

n_used = 0

do i = 1, n_key
  if (i == overlay_lord$) cycle
  print *, '!---------------------------------'
  print *, key_name(i)
  ele%key = i
  do j = 1, n_attrib_special_maxx
    a_name = attribute_name (ele, j) 
    if (a_name == null_name$) cycle
    print '(i10, 2x, a)', j, a_name
    n_used(j) = n_used(j) + 1
  enddo
  print *
enddo

! number of elements using an index

print *, '!---------------------------------'
print *, 'Index usage:'
print *, '   Ix Count'
do i = 1, size(n_used)
  print '(2i6)', i, n_used(i)
enddo

! List of elements using an index

do i = 1, n_attrib_special_maxx
  print *
  print *, '!---------------------------------'
  print '(a, i0)', 'Index: ', i
  do j = 1, n_key
    if (j == overlay_lord$) cycle
    ele%key = j
    a_name = attribute_name (ele, i) 
    if (a_name == null_name$) cycle
    print *, '   ', key_name(ele%key)
  enddo
enddo

end program
