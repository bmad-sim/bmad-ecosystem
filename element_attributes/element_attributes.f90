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
  type (entry_struct) entry(num_ele_attrib_extended$+1)
end type

type (table_struct), target :: table(n_key$)
type (table_struct), pointer :: tab
type (entry_struct), pointer :: e1, e2, e3, e4
type (ele_struct) ele
type (ele_attribute_struct) attrib

integer i, ik, j, n, n_used(num_ele_attrib_extended$), key_indx(n_key$)
integer ie, it, n_table, n_row, n_char, ios, nl
integer indx(num_ele_attrib_extended$)

character(20) date, arg
character(40) cap_name, name
character(200) line

!

call cesr_getarg(1, arg)
call init_attribute_name_array()
call indexx(key_name, key_indx)

! Element keys

if (arg /= 'tex') then

  open (1, file = 'element_attributes.dat')

  call date_and_time_stamp (date, .true.)

  write (1, *) date
  write (1, *)

  do i = 1, n_key$
    j = key_indx(i)
    write (1, '(i4, 3x, a, i4, 3x, a)') i, key_name(i), j, trim(key_name(j))
  enddo

  write (1, *) 

  ! indexes used for a given element key

  n_used = 0

  do ik = 1, n_key$
    i = key_indx(ik)
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
      ele%key = j
      attrib = attribute_info(ele, i)
      if (attrib%name == null_name$) cycle
      write (1, *) '   ', key_name(ele%key)
    enddo
  enddo

  print *, 'Writen: element_attributes.dat'

  close (1)
endif

!-----------------------------------------------------
! Create LaTex file

! Setup table entries.

it = 0
n_char = 0

do n = 1, n_key$
  i = key_indx(n)

  it = it + 1
  tab => table(it)
  tab%name = key_name(i)
  tab%key = i
  ele%key = i
  ie = 0

  do j = 1, a0$
    attrib = attribute_info(ele, j)
    if (attrib%type == private$) cycle

    select case (attrib%name)
    case (null_name$, 'ELE_BEGINNING', 'ELE_CENTER', 'ELE_END', &
          'REF_BEGINNING', 'REF_CENTER', 'REF_END')
      cycle
    end select

    if (attrib%name == 'CURVATURE_Z0_Y2') attrib%name = 'CURVATURE_Z0_Y2, ..., etc.'
    if (attrib%name == 'A0') attrib%name = 'A0 - A20, B0 - B20'
    if (attrib%name == 'K0L') attrib%name = 'K0L - K20L, T0 - T20'

    ie = ie + 1
    tab%n_line = ie

    write (tab%entry(ie)%indx, '(i3)') j
    tab%entry(ie)%name = attrib%name
    n_char = max(n_char, len_trim(attrib%name))


!    if (attrib%type == private$)        tab%entry(ie)%attrib = '[p]'
!    if (attrib%type == dependent$)      tab%entry(ie)%attrib = '[d]'

  enddo
enddo

n_table = it

! Now write the file

open (1, file = 'list-element-attributes.tex')

write (1, '(a)') '\chapter{List of Element Attributes}'
write (1, '(a)') '\label{c:attrib.list}'
write (1, '(a)') ''
write (1, '(a)') 'Alphabetical list of element attributes for each type of element. '
write (1, '(a)') ''
write (1, '(a)') 'Note for programmers: The program that generates a file of attributes indexed by the'
write (1, '(a)') 'internal reference number is:'
write (1, '(a)') '\begin{example}'
write (1, '(a)') '  util_programs/element_attributes.f90 '
write (1, '(a)') '\end{example}'
write (1, '(a)') ''

!

do it = 1, n_table, 1
  tab => table(it)
  tab%name = key_name(tab%key)

  if (tab%key == def_bmad_com$) cycle
  if (tab%key == def_beam_start$) cycle
  if (tab%key == def_mad_beam$) cycle
  if (tab%key == def_parameter$) cycle
  if (tab%key == line_ele$) cycle
  if (tab%key == null_ele$) cycle
  if (tab%key == rbend$) cycle
  if (tab%key == rcollimator$) cycle
  if (tab%key == photon_fork$) cycle
  if (tab%key == monitor$) cycle
  if (tab%key == pipe$) cycle
  if (tab%key == vkicker$) cycle
  if (tab%key == undulator$) cycle
  if (tab%key == beginning_ele$) cycle

  nl = tab%n_line
  indx = nl+1  ! point to blank attribute name
  call indexx(tab%entry(1:nl)%name, indx(1:nl))


  write (1, *) '%---------------------------------'

  select case (tab%key)
  case (sbend$)
    write (1, *) '\section{Bends: Rbend and Sbend Element Attributes}'
    write (1, *) '\label{s:list.bend}'
  case (ecollimator$)
    write (1, *) '\section{Collimators: Ecollimator and Rcollimator Element Attributes}'
    write (1, *) '\label{s:list.collimator}'
  case (fork$)
    write (1, *) '\section{Fork and Photon_Fork Element Attributes}'
    write (1, *) '\label{s:list.fork}'
  case (instrument$)
    write (1, *) '\section{Instrument, Monitor, and Pipe Element Attributes}'
    write (1, *) '\label{s:list.instrument}'
  case (hkicker$)
    write (1, *) '\section{Kickers: Hkicker and Vkicker Element Attributes}'
    write (1, *) '\label{s:list.hvkicker}'
  case (wiggler$)
    write (1, *) '\section{:Wiggler and Undulator Element Attributes}'
    write (1, *) '\label{s:list.wiggler}'
  case default
    write (1, *) '\section{', trim(tab%name), ' Element Element Attributes}'
    name = downcase(key_name(tab%key))
    call str_substitute (name, '_', '.')
    write (1, *) '\label{s:list.', trim(name), '}'
  end select

  write (1, *)
  write (1, *) '\begin{tabular}{llll} \toprule'
  ele%key = tab%key
  n_row = (tab%n_line+3)/4
  do ie = 1, n_row
    e1 => tab%entry(indx(ie))
    e2 => tab%entry(indx(ie+n_row))
    e3 => tab%entry(indx(ie+2*n_row))
    e4 => tab%entry(indx(ie+3*n_row))
    write (1, '(14a)') downcase(e1%name(1:n_char)), ' & ', downcase(e2%name(1:n_char)), ' & ', &
                       downcase(e3%name(1:n_char)), ' & ', downcase(e4%name(1:n_char)), ' \\'
  enddo
  write (1, *) '\bottomrule'
  write (1, *) '\end{tabular}'
  write (1, *) '\vfill'
  write (1, *) 

enddo



end program
