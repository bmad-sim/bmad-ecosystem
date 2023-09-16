!+
! Program: element_attributes
!
! Program to print which element attributes are defined for different types of elements.
! This is useful when new attributes need to be defined for a given element type.
! Generally only Bmad wizards are interested in this.
!-

program element_attributes

use attribute_mod

implicit none

type entry_struct
  character(3) :: indx = ''
  character(40) :: name = ''
  character(20) :: units = ''
end type

type table_struct
  character(100) :: sort_name = ''
  character(40) :: label_ref = ''
  integer :: key = -1
  integer :: n_line = 0
  type (entry_struct) entry(num_ele_attrib_extended$+1)
end type

type (table_struct), target :: table(n_key$)
type (table_struct), pointer :: tab
type (entry_struct), pointer :: ee, e1, e2, e3, e4
type (ele_struct) ele
type (ele_attribute_struct) attrib

integer i, ik, j, n, n_used(num_ele_attrib_extended$), key_indx(n_key$)
integer ix, ie, it, n_row, n_char, ios, nl
integer indx(num_ele_attrib_extended$)

character(20) date, arg
character(40) cap_name, name
character(200) line

!

call get_command_argument(1, arg)
call init_attribute_name_array()
call indexer(key_name, key_indx)

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
      if (attrib%state == private$) then
        write (1, '(i10, 2x, 2a)') j, attrib%name, '  [private]'
      else
        write (1, '(i10, 2x, a)') j, attrib%name
        if (attribute_type(attrib%name, ele) == is_real$ .and. attribute_units(attrib%name, '???') == '???') then
          print *, 'No units for: ' // trim(attrib%name)
        endif
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

  print *, 'Written: element_attributes.dat'

  close (1)
endif

!-----------------------------------------------------
! Create LaTex file

! Setup table entries.

it = 0
n_char = 0

do n = 1, n_key$

  tab => table(n)
  tab%key = n
  ele%key = n

  select case (tab%key)
  case (sbend$)
    tab%sort_name = 'Bends: Rbend and Sbend Element Attributes'
    tab%label_ref = 'bend'
  case (ecollimator$)
    tab%sort_name = 'Collimators: Ecollimator and Rcollimator Element Attributes'
    tab%label_ref = 'collimator'
  case (fork$)
    tab%sort_name = 'Fork and Photon_Fork Element Attributes'
    tab%label_ref = 'fork'
  case (instrument$)
    tab%sort_name = 'Instrument, Monitor, and Pipe Element Attributes'
    tab%label_ref = 'instrument'
  case (hkicker$)
    tab%sort_name = 'Kickers: Hkicker and Vkicker Element Attributes'
    tab%label_ref = 'hvkicker'
  case (wiggler$)
    tab%sort_name = 'Wiggler and Undulator Element Attributes'
    tab%label_ref = 'wiggler'
  case (def_parameter$)
    tab%sort_name = 'Parameter Statement Attributes'
    tab%label_ref = 'parameter'
  case (def_line$)
    tab%sort_name = 'Line Statement Attributes'
    tab%label_ref = 'line'
  case (def_particle_start$)
    tab%sort_name = 'Particle_Start Statement Attributes'
    tab%label_ref = 'particle.start'
  case (def_bmad_com$)
    tab%sort_name = 'Bmad_Com Statement Attributes'
    tab%label_ref = 'bmad.com'
  case (beginning_ele$)
    tab%sort_name = 'Beginning Statement Attributes'
    tab%label_ref = 'beginning'
  case default
    tab%sort_name = trim(key_name(tab%key)) // ' Element Attributes'
    tab%label_ref = downcase(key_name(tab%key))
    call str_substitute (tab%label_ref, '_', '.')
  end select

  ie = 0

  do j = 1, a0$
    attrib = attribute_info(ele, j)
    if (attrib%state == private$) cycle

    select case (attrib%name)
    case (null_name$, 'ELE_BEGINNING', 'ELE_CENTER', 'ELE_END', &
          'REF_BEGINNING', 'REF_CENTER', 'REF_END')
      cycle
    end select

    ie = ie + 1
    tab%n_line = ie

    ee => tab%entry(ie)
    write (ee%indx, '(i3)') j

    select case (attrib%name)
    case ('CURVATURE_Z0_Y2'); ee%name = 'curvature_z0_y2, ..., etc.'
    case ('A0');              ee%name = 'a0 - a20, b0 - b20'
    case ('K0L');             ee%name = 'k0l - k20l, t0 - t20'
    case default
      ee%name = downcase(attrib%name)
      if (attribute_type(attrib%name) == is_real$) then
        ee%units = attribute_units(attrib%name, '???')
        if (ee%units /= '') then
          ee%name = trim(ee%name) // ' [' // trim(ee%units) // ']'
          ix = index(ee%name, '^')
          if (ix /= 0) then       ! Something like T/m^3 . The exponent is always a single digit.
            ee%name = ee%name(:ix-1) // '$' // ee%name(ix:ix+1) // '$' // ee%name(ix+2:)
          endif
        endif
      endif
    end select

    n_char = max(n_char, len_trim(ee%name))
  enddo
enddo

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

call indexer(table%sort_name, key_indx)

do it = 1, n_key$
  tab => table(key_indx(it))

  if (tab%key == def_mad_beam$) cycle
  if (tab%key == null_ele$) cycle
  if (tab%key == rbend$) cycle
  if (tab%key == rcollimator$) cycle
  if (tab%key == photon_fork$) cycle
  if (tab%key == monitor$) cycle
  if (tab%key == pipe$) cycle
  if (tab%key == vkicker$) cycle
  if (tab%key == undulator$) cycle

  nl = tab%n_line
  indx = nl+1  ! point to blank attribute name
  call indexer(tab%entry(1:nl)%name, indx(1:nl))


  write (1, *) '%---------------------------------'

  write (1, *) '\section{', trim(tab%sort_name), '}'
  write (1, *) '\label{s:list.', trim(tab%label_ref), '}'
  write (1, *)
  write (1, *) '\begin{tabular}{llll} \toprule'
  ele%key = tab%key
  n_row = (tab%n_line+3)/4
  do ie = 1, n_row
    e1 => tab%entry(indx(ie))
    e2 => tab%entry(indx(ie+n_row))
    e3 => tab%entry(indx(ie+2*n_row))
    e4 => tab%entry(indx(ie+3*n_row))
    write (1, '(14a)') e1%name(1:n_char), ' & ', e2%name(1:n_char), ' & ', e3%name(1:n_char), ' & ', e4%name(1:n_char), ' \\'
  enddo
  write (1, *) '\bottomrule'
  write (1, *) '\end{tabular}'
  write (1, *) '\vfill'
  write (1, *) 
enddo

end program
