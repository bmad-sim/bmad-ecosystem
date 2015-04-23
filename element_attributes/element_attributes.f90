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

type entry_struct
  character(3) :: indx = ''
  character(40) :: name = ''
  character(3) :: attrib = ''
end type

type table_struct
  character(40) :: name = ''
  integer :: key = -1
  integer :: n_line = 0
  type (entry_struct) entry(num_ele_attrib_extended$)
end type

type (table_struct), target :: table(n_key$)
type (entry_struct), pointer :: e1, e2

type (ele_struct) ele
type (ele_attribute_struct) attrib
integer i, ik, j, n, n_used(num_ele_attrib_extended$), key_indx(n_key$)
integer ie, it, n_table, n_row, n_char
character(20) date
character(40) down_name, cap_name

! Element keys

open (1, file = 'element_attributes.dat')

call date_and_time_stamp (date, .true.)

write (1, *) date
write (1, *)

call indexx(key_name, key_indx)

do i = 1, n_key$
  j = key_indx(i)
  write (1, '(i4, 3x, a, i4, 3x, a)') i, key_name(i), j, trim(key_name(j))
enddo

write (1, *) 

! indexes used for a given element key

n_used = 0

do ik = 1, n_key$
  i = key_indx(ik)
  if (i == overlay$) cycle
  write (1, *) '!---------------------------------'
  write (1, '(i3, 2x, a)') i, key_name(i)
  ele%key = i
  do j = 1, a0$
    attrib = attribute_info(ele, j)
    if (attrib%name == null_name$) cycle
    if (attrib%type == private$) then
      write (1, '(i10, 2x, 2a)') j, attrib%name, '  [private]'
    else
      write (1, '(i10, 2x, a)') j, attrib%name
    endif
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
  do j = 1, n_key$
    if (j == overlay_lord$) cycle
    ele%key = j
    attrib = attribute_info(ele, i)
    if (attrib%name == null_name$) cycle
    write (1, *) '   ', key_name(ele%key)
  enddo
enddo

print *, 'Writen: element_attributes.dat'

close (1)

!-----------------------------------------------------
! Create LaTex file

! Setup table entries.

it = 0
n_char = 0

do n = 1, n_key$
  i = key_indx(n)
  if (i == overlay$) cycle

  it = it + 1
  table(it)%name = key_name(i)
  table(it)%key = i
  ele%key = i
  ie = 0
  do j = 1, a0$
    attrib = attribute_info(ele, j)

    select case (attrib%name)
    case (null_name$, 'ELE_BEGINNING', 'ELE_CENTER', 'ELE_END', &
          'REF_BEGINNING', 'REF_CENTER', 'REF_END')
      cycle
    end select

    if (attrib%name == 'CURVATURE_Z0_Y2') attrib%name = 'CURVATURE_Z0_Y2, ..., etc.'
    if (attrib%name == 'A0') attrib%name = 'A0 - A20, B0 - B20'
    if (attrib%name == 'K0L') attrib%name = 'K0L - K20L, T0 - T20'

    ie = ie + 1
    write (table(it)%entry(ie)%indx, '(i3)') j
    table(it)%entry(ie)%name = attrib%name
    n_char = max(n_char, len_trim(attrib%name))

    if (attrib%type == private$)        table(it)%entry(ie)%attrib = '[p]'
    if (attrib%type == dependent$)      table(it)%entry(ie)%attrib = '[d]'

    table(it)%n_line = table(it)%n_line + 1
  enddo
enddo

n_table = it

! Now write the file

open (1, file = 'list-element-attributes.tex')

do it = 1, n_table, 1
  down_name = downcase(table(it)%name)

  write (1, *) '%---------------------------------'
  write (1, *) '\section{', trim(down_name), '}'
  write (1, *) '\label{s:list.', trim(down_name), '}'
  write (1, *)
  write (1, *) '\begin{tabular}{lll||lll} \hline'
  ele%key = table(it)%key
  n_row = (table(it)%n_line+1)/2
  do ie = 1, n_row
    e1 => table(it)%entry(ie)
    e2 => table(it)%entry(ie+n_row)
    write (1, '(14a)') e1%indx, ' & ', e1%attrib, ' & ', e1%name(1:n_char), ' & ', &
                       e2%indx, ' & ', e2%attrib, ' & ', e2%name(1:n_char), ' \HH'
  enddo
  write (1, *) '\end{tabular}'
  write (1, *) '\vfill \break'
  write (1, *)
enddo



end program
