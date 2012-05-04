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
integer i, ik, j, n_used(n_attrib_special_maxx), key_indx(n_key)
character(40) a_name
character(20) date

! Element keys

open (1, file = 'element_attributes.dat')

call date_and_time_stamp (date, .true.)

write (1, *) date
write (1, *)

call indexx(key_name, key_indx)

do i = 1, n_key
  j = key_indx(i)
  write (1, '(i4, 3x, a, i4, 3x, a)') i, key_name(i), j, trim(key_name(j))
enddo

write (1, *) 

! indexes used for a given element key

n_used = 0

do ik = 1, n_key
  i = key_indx(ik)
  if (i == overlay_lord$) cycle
  write (1, *) '!---------------------------------'
  write (1, '(i3, 2x, a)') i, key_name(i)
  ele%key = i
  do j = 1, a0$
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

do i = 1, a0$
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

close (1)

!-----------------------------------------------------

open (1, file = 'element-attributes-list.tex')

do i = 1, n_key
  if (i == overlay_lord$) cycle
  write (1, *) '!---------------------------------'
  write (1, '(i3, 2x, a)') i, key_name(i)
  ele%key = i
  do j = 1, a0$
    a_name = attribute_name (ele, j) 
    if (a_name == null_name$) cycle
    write (1, '(i10, 2x, a)') j, a_name
    n_used(j) = n_used(j) + 1
  enddo
  write (1, *)
enddo



end program
