module aml_parser_mod

use bmad_struct
use bmad_interface
use bmad_parser_mod
use uap_fortran
use random_mod

contains

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

subroutine transfer_attribute_info_to_ele (lat, node, ele)

implicit none

type (lat_struct) lat
type (uap_node_struct), target :: node
type (uap_node_struct), pointer :: attrib_node
type (uap_attribute_struct), pointer :: attrib
type (ele_struct) ele

integer i, ix, ios

character(40) value, str
character(32) :: r_name = 'transfer_attribute_info_to_ele'
logical found, err_flag

!

found = parser_set_attribute_value(node, 'name', value)
call str_upcase(ele%name, value)

node%ix = ele%ix_ele

call get_child_by_name (node, "bmad_attributes", attrib_node)
found = parser_set_attribute_value (attrib_node, 'class', value)
ele%key = key_name_to_key_index(value)
if (ele%key < 0) then
  call parser_error('CANNOT DETERMINE PROPER ELEMENT CLASS FOR ELEMENT: ' // ele%name)
  return
endif

do i = lbound(attrib_node%attributes, 1), ubound(attrib_node%attributes, 1)

  attrib => attrib_node%attributes(i)

  select case (attrib%name)
  case ('class', 'name', 'value', 'attribute')
    cycle
  case ('descrip')
    ele%descrip = attrib%value
    cycle
  case ('alias')
    ele%alias = attrib%value
    cycle
  case ('type')
    ele%type = attrib%value
    cycle
  end select

  if (attrib%name(1:3) == 'bl_' .or. attrib%name == 'b_' .or. &
      attrib%name(4:11) == 'gradient') ele%field_master = .true.

  call str_upcase (str, attrib%name)
  ix = attribute_index(ele, str)
  if (ix < 1) then
    call parser_error ('BAD ATTRIBUTE: ' // attrib%name, &
                                   'OF ELEMENT: ' // ele%name)
    cycle
  endif

  read (attrib%value, *, iostat = ios) ele%value(ix) 
  if (ios /= 0) then
    call parser_error ( &
                  'CANNOT CONVERT ATTRIBUTE: ' // attrib%name, &
                  'OF ELEMENT: ' // ele%name, &
                  'TO REAL NUMBER. ATTRIBUTE STRING IS: ' // attrib%value)
  endif

enddo

end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

subroutine read_parameter (param_node, value_str, value)

implicit none

type (uap_node_struct) param_node
character(*) value_str
real(rp) value
integer ios

!

read (value_str, *, iostat = ios) value
if (ios /= 0) then
  call parser_error ('BAD PARAMETER VALUE: ', node = param_node)
endif

end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

subroutine parser_error (what1, what2, what3, node)

implicit none

character(*) what1
character(*), optional :: what2, what3
type (uap_node_struct), optional :: node

! BP_COM%ERROR_FLAG is a common logical used so program will stop at end of parsing

print *, 'ERROR IN ', trim(bp_com%parser_name), ': ', trim(what1)

if (present(what2)) print '(22x, a)', trim(what2)
if (present(what3)) print '(22x, a)', trim(what3)
if (present(node)) call uap_print_node (node)

print *

bp_com%error_flag = .true.

end subroutine

end module

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine aml_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
!
! Subroutine to parse an AML input file and put the information in lat.
!
! Because of the time it takes to parse a file, AML_PARSER will save 
! LAT in a "digested" file with the name:
!               'digested_' // lat_file   ! For single precision AML version
! For subsequent calls to the same lat_file, AML_PARSER will just read in the
! digested file. AML_PARSER will always check to see that the digested file
! is up-to-date and if not the digested file will not be used.
!
! Input:
!   lat_file   -- Character(*): Name of the input file.
!   make_mats6 -- Logical, optional: Compute the 6x6 transport matrices for the
!                   Elements? Default is True.
!   use_line   -- Character(*), optional: If present then override the use 
!                   statement in the lattice file and use use_line instead.
!
! Output:
!   lat -- lat_struct: Lat structure. See AML_struct.f90 for more details.
!     %ele(:)%mat6  -- This is computed assuming an on-axis orbit 
!     %ele(:)%s     -- This is also computed.
!   digested_read_ok -- Logical, optional: Set True if the digested file was
!                        successfully read. False otherwise.
!
! Defaults:
!   lat%param%particle          = positron$
!   lat%param%geometry      = closed$
!-

subroutine aml_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)

use aml_parser_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, slave_ele
type (ele_struct) param_ele
type (uap_node_struct) expand_node
type (uap_node_struct), pointer :: control_node, machine_node, slave_node
type (uap_node_struct), pointer :: track_node, master_node, node, controller_node
type (uap_node_struct), pointer :: attrib_node, param_node, child_node
type (ele_struct), allocatable :: old_ele(:) 
type (all_pointer_struct), allocatable :: ptrs(:)

real(rp) dflt_coef, con_value, coef, value

integer digested_version, c_ok
integer i, j, k, n, ix, ios, iu, ix_lord
integer n_track, n_control, n_master, n_slave, n_con1, n_con2
integer ix_slave, ix_ele, ix0_att

character(*) lat_file
character(*), optional :: use_line

character(200) full_lat_file_name, digested_file
character(100) coef_str
character(40) name, value_str, slave_rank, lord_rank, dflt_attrib
character(40) con_value_str, attrib_str, param_str, class_str
character(16) :: r_name = 'aml_parser'
character(80) debug_line

logical, optional :: make_mats6, digested_read_ok
logical found, err_flag
logical :: debug = .false.

! see if digested file is open and current. If so read in and return.
! Note: The name of the digested file depends upon the real precision.

bp_com%write_digested = .true.
debug_line = ''
bp_com%parser_name = "AML_PARSER"

call form_digested_bmad_file_name (lat_file, digested_file, full_lat_file_name)
call read_digested_bmad_file (digested_file, lat, digested_version, err_flag, .true.)

! Must make sure that if use_line is present the digested file has used the 
! correct line

if (present(use_line)) then
  call str_upcase (name, use_line)
  if (name /= lat%use_name) bmad_status%ok = .false.
endif

if (bmad_status%ok) then
  if (lat%input_taylor_order /= 0) ptc_com%taylor_order_saved = lat%input_taylor_order
  call set_ptc (lat%ele(0)%value(e_tot$), lat%param%particle)
  if (lat%input_taylor_order == bmad_com%taylor_order) then
    if (present(digested_read_ok)) digested_read_ok = .true.
    return
  else
     call out_io (s_info$, r_name, 'Taylor_order has changed.', &
         'Taylor_order in digested file: \i4\ ', &
         'Taylor_order now:              \i4\ ', &
         i_array = [lat%input_taylor_order, bmad_com%taylor_order ])
     bmad_status%ok = .false.
    if (lat%input_taylor_order > bmad_com%taylor_order) bp_com%write_digested = .false.
  endif
endif

if (present(digested_read_ok)) digested_read_ok = .false.

! save all elements that have a taylor series

call save_taylor_elements (lat, old_ele)

!------------------------------------------------------------------------
! here if not OK bmad_status. So we have to do everything from scratch...

! Some defaults

lat%param%particle = positron$
lat%param%stable = .true.

! Put lattice information into root_node

call aml_parser_cpp (c_str(lat_file), expand_node, c_ok)

bmad_status%ok = f_logic(c_ok)

if (.not. bmad_status%ok) then
  if (global_com%exit_on_error) then
    call out_io (s_fatal$, r_name, 'AML_PARSER FINISHED. EXITING ON ERRORS')
    if (global_com%exit_on_error) call err_exit
  endif
  return
endif

if (debug) then
  iu = lunget()
  open (iu, file = 'expand_node.dat')
  call uap_print_tree (expand_node, iu)
  close (iu)
  call out_io (s_info, r_name, 'Written: expand_node.dat')
endif

! Transfer the information in the node tree to lat...

call get_child_by_name (expand_node, "control_list", control_node)
call get_child_by_name (expand_node, "machine", machine_node)
call get_child_by_name (machine_node, "tracking_lattice", track_node)
call get_child_by_name (machine_node, "master_list", master_node)
call get_child_by_name (expand_node, "param_list", param_node)

! The tracking lattice has a 0th node which is not counted.

n_track   = size(track_node%children) - 1
n_master  = size(master_node%children)
n_control = size(control_node%children) 

call init_lat (lat, n_track + n_master + n_control)
lat%n_ele_track = n_track
lat%n_ele_max = n_track 

! Init some stuff

lat%ele(0)%name = 'BEGINNING'
lat%ele(0)%key  = beginning_ele$
lat%param%geometry = closed$
lat%param%particle = positron$


!------------------------------------------------------------
! Transfer the tracking lattice attribute info

ix_ele = 0
do i = lbound(track_node%children, 1) + 1, ubound(track_node%children, 1)

  node => track_node%children(i)
  ix_ele = ix_ele + 1
  ele => lat%ele(ix_ele)

  ! Find the ele%slave_status. 

  found = parser_set_attribute_value (node, 'slave_rank', slave_rank)
  select case (slave_rank)

  case ('SUPER_SLAVE')
    ele%slave_status = super_slave$

  case ('MULTIPASS_SLAVE')
    ele%slave_status = multipass_slave$

  case ('')
    ele%slave_status = free$

  case default
    call parser_error('UNKNOWN SLAVE_RANK: ' // slave_rank)
  end select

  call transfer_attribute_info_to_ele (lat, node, ele)

  ! The length of a super slave needs to be saved

  if (slave_rank == 'SUPER_SLAVE') then
    call get_child_by_name (node, 'length', child_node)
    if (.not. associated(child_node)) cycle
    if (.not. parser_set_attribute_value (child_node, "actual", attrib_str)) cycle
    read (attrib_str, *) ele%value(l$)
  endif

enddo

!------------------------------------------------------------
! Transfer the master attribute info.
! super_lord drifts are excluded.


do i = lbound(master_node%children, 1), ubound(master_node%children, 1)
  node => master_node%children(i)

  lat%n_ele_max = lat%n_ele_max + 1
  ele => lat%ele(lat%n_ele_max)

  ! Find the ele%lord_type

  found = parser_set_attribute_value (node, 'lord_rank', lord_rank)
  select case (lord_rank)
  case ('SUPER_LORD')
    ele%lord_status = super_lord$
  case ('MULTIPASS_LORD')
    ele%lord_status = multipass_lord$
  case default
    call parser_error ('UNKNOWN LORD STATUS: ' // lord_rank, node = node)
  end select

  call transfer_attribute_info_to_ele (lat, node, ele)

  if (ele%key == drift$) then
    lat%n_ele_max = lat%n_ele_max - 1
    node%ix = -1
  endif

enddo

!------------------------------------------------------------
! Transfer the control attribute info

do i = lbound(control_node%children, 1), ubound(control_node%children, 1)
  node => control_node%children(i)

  lat%n_ele_max = lat%n_ele_max + 1
  ele => lat%ele(lat%n_ele_max)

  call transfer_attribute_info_to_ele (lat, node, ele)

  ! Find the ele%lord_status

  select case (ele%key)
  case (overlay$)
    ele%lord_status = overlay_lord$
  case (group$)
    ele%lord_status = group_lord$
  case default
    call parser_error ('UNKNOWN CONTROLLER LORD STATUS: ' // key_name(ele%key), node = node)
  end select

enddo

!------------------------------------------------------------
! Transfer the master slave info

do i = lbound(master_node%children, 1), ubound(master_node%children, 1)
  node => master_node%children(i)
  if (node%ix == -1) cycle   ! Ignore drifts

  ele => lat%ele(node%ix)
  ele%ixx = 0

  n_slave = size(node%slaves)
  if (n_slave == 0) return

  n_con1 = lat%n_control_max + 1
  n_con2 = n_con1 + n_slave - 1

  if (n_con2 >  size(lat%control)) call reallocate_control(lat, n_con2)
  lat%n_control_max = n_con2

  ele%n_slave = n_slave
  ele%ix1_slave = n_con1

  do j = lbound(node%slaves, 1), ubound(node%slaves, 1)
    n = j + n_con1 - lbound(node%slaves, 1)
    lat%control(n)%lord%ix_ele = ele%ix_ele
    ix_ele = node%slaves(j)%node%ix
    slave_ele => lat%ele(ix_ele)
    lat%control(n)slave%ix_ele = ix_ele
    lat%control(n)%coef = slave_ele%value(l$) / ele%value(l$)
  end do

enddo

!------------------------------------------------------------
! Transfer the control slave info

do i = lbound(control_node%children, 1), ubound(control_node%children, 1)
  node => control_node%children(i)
  ele => lat%ele(node%ix)

  call get_child_by_name (node, 'bmad_attributes', attrib_node)
  found = parser_set_attribute_value (attrib_node, 'value', con_value_str)
  found = parser_set_attribute_value (attrib_node, 'coef', coef_str)
  found = parser_set_attribute_value (attrib_node, 'attribute', dflt_attrib)

  if (ele%key == overlay$) then
    call str_upcase (dflt_attrib, dflt_attrib)
    ele%component_name = dflt_attrib
    ele%ix_value = attribute_index(ele, dflt_attrib)
    if (ele%ix_value < 1) then
      call parser_error ('BAD ATTRIBUTE: ' // dflt_attrib, 'FOR CONTROLLER: ' // ele%name)
    endif
  endif

  dflt_coef = 1
  if (coef_str /= '') then
    read (coef_str, *, iostat = ios) dflt_coef
    if (ios /= 0) then
      call parser_error ('BAD CONTROLLER COEF: ' // coef_str, &
                         'FOR CONTROLLER: ' // ele%name)
    endif
  endif

  con_value = 0
  if (con_value_str /= '') then
    read (con_value_str, *, iostat = ios) con_value
    if (ios /= 0) then
      call parser_error ('BAD CONTROLLER VALUE: ' // con_value_str, &
                         'FOR CONTROLLER: ' // ele%name)
    endif
  endif
  
  n_slave = 0
  do j = lbound(node%children, 1), ubound(node%children, 1)
    slave_node => node%children(j)
    n_slave = n_slave + size(slave_node%slaves)
  enddo
  if (n_slave == 0) return

  n_con1 = lat%n_control_max + 1
  n_con2 = n_con1 + n_slave - 1

  if (n_con2 >  size(lat%control)) call reallocate_control(lat, n_con2)
  lat%n_control_max = n_con2

  ele%n_slave = n_slave
  ele%ix1_slave = n_con1

  do j = lbound(node%children, 1), ubound(node%children, 1)
    slave_node => node%children(j)
    if (slave_node%name /= 'slave') cycle

    do k = lbound(slave_node%slaves, 1), ubound(slave_node%slaves, 1)
      n = k + n_con1 - lbound(slave_node%slaves, 1)
      ix_slave = slave_node%slaves(k)%node%ix
      lat%control(n)slave%ix_ele = ix_slave
      lat%control(n)%lord%ix_ele = ele%ix_ele

      found =  parser_set_attribute_value (attrib_node%children(j), 'coef', coef_str)
      coef = dflt_coef
      if (coef_str /= '') then
        read (coef_str, *, iostat = ios) coef_str
        if (ios /= 0) then
          call parser_error ( &
                          'BAD CONTROLLER COEF: ' // coef_str, &
                          'FOR CONTROLLER: ' // ele%name)
        endif
      endif
      lat%control(n)%coef = coef

      found = parser_set_attribute_value (attrib_node%children(j), 'attribute', attrib_str)
      if (attrib_str == '') attrib_str = dflt_attrib
      call str_upcase (attrib_str, attrib_str)
      ix = attribute_index (lat%ele(ix_slave), attrib_str)
      if (ix < 1) then
        call parser_error ('BAD CONTROLLER ATTRIBUTE: ' // attrib_str, &
                           'FOR CONTROLLER: ' // ele%name)
      endif
      lat%control(n)%ix_attrib = ix

    enddo  ! k

    n_con1 = n_con1 + size(slave_node%slaves)
  end do

enddo

!------------------------------------------------------------
! Add the lord info for controllers and masters

do i = 1, lat%n_control_max
  ix = lat%control(i)%lord%ix_ele
  if (lat%ele(ix)%key == group$) cycle
  ix = lat%control(i)slave%ix_ele
  lat%ele(ix)%n_lord = lat%ele(ix)%n_lord + 1
enddo

n = 0
do i = 1, lat%n_ele_max
  slave_ele => lat%ele(i)
  if (slave_ele%n_lord == 0) cycle
  slave_ele%ic1_lord = n + 1
  slave_ele%ixx = n + 1
  n = n + slave_ele%n_lord
enddo

do i = 1, lat%n_control_max
  ix = lat%control(i)%lord%ix_ele
  if (lat%ele(ix)%key == group$) cycle
  ix = lat%control(i)slave%ix_ele
  slave_ele => lat%ele(ix)
  if (slave_ele%n_lord == 0) cycle
  j = slave_ele%ixx
  lat%ic(j) = i
  slave_ele%ixx = j + 1
enddo

!------------------------------------------------------------
! Find super_slave key and correct the name.

do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (ele%slave_status /= super_slave$) cycle

  ele%name = ''
  do j = ele%ic1_lord, ele%ic1_lord+ele%n_lord-1
    ix = lat%ic(j)
    ix_lord = lat%control(ix)%lord%ix_ele
    ele%name = trim(ele%name) // '\' // trim(lat%ele(ix_lord)%name)  ! '
    call calc_superimpose_key (ele, lat%ele(ix_lord), ele)
    if (ele%key < 1) then
      call parser_error ('BAD SUPERPOSITION OF ELEMENTS! ' // ele%name)
    endif
  enddo

  if (ele%n_lord == 1) lat%ele(ix_lord)%ixx = lat%ele(ix_lord)%ixx + 1
  write (ele%name, '(2a, i0)') trim(ele%name), '\', lat%ele(ix_lord)%ixx   ! '

enddo

!------------------------------------------------------------
! Transfer the param info

do i = lbound(param_node%children, 1), ubound(param_node%children, 1)

  child_node => param_node%children(i)
  ix0_att = lbound(child_node%attributes, 1)

  call str_upcase (class_str, child_node%name)
  select case (class_str)
  case ('BEGINNING');  param_ele%key = beginning_ele$
  case ('PARTICLE_START'); param_ele%key = def_particle_start$
  case ('PARAMETER');  param_ele%key = def_parameter$
  case default
    call parser_error ('BAD PARAMETER GROUP: ', node = child_node)
    cycle
  end select

  call str_upcase (param_str, child_node%attributes(ix0_att)%name)
  value_str = child_node%attributes(ix0_att)%value

  select case (class_str)
  case ('BEGINNING', 'PARTICLE_START')
    call pointers_to_attribute (lat, class_str, param_str, .true., ptrs, err_flag)
    if (err_flag) cycle
    ptrs(1)%r = value
  case ('PARAMETER')
    select case (param_str)
    case ('N_PART')
      call read_parameter (child_node, value_str, value)
      lat%param%n_part = value
    case ("GEOMETRY")
      call match_word(value_str, geometry, ix)
      if (ix == 0) then
        call parser_error ('BAD GEOMETRY', node = child_node)
        cycle 
      endif
      lat%param%geometry = ix + (lbound(geometry, 1) - 1)
    case ("LATTICE")
      lat%lattice = value_str
    case ("E_TOT")
      call read_parameter (child_node, value_str, value)
      lat%ele(0)%value(e_tot$) = value
    case ("TAYLOR_ORDER")
      call read_parameter (child_node, value_str, value)
      lat%input_taylor_order = nint(value)
    case ("RAN_SEED")
      call read_parameter (child_node, value_str, value)
      call ran_seed_put (nint(value))
    case ("PARTICLE")
      call match_word(value_str, species_name, ix)
      if (ix < 1) then
        call parser_error ('BAD PARTICLE TYPE: ', value_str)
        cycle
      endif
      lat%param%particle = ix - lbound(species_name, 1) + 1

    case default
      call parser_error ('UNKNOWN PARAMETER: ' // param_str, node = child_node)
      
    end select
  end select

enddo

! End bookkeeping

call s_calc(lat)
lat%version = bmad_inc_version$

call lattice_bookkeeper (lat)
lat%input_taylor_order = bmad_com%taylor_order

!------------------------------------------------------------

if (debug_line /= '') call parser_debug_print_info (lat, debug_line)

! error check

if (bp_com%error_flag) then
  if (global_com%exit_on_error) then
    call out_io (s_fatal$, r_name, 'AML_PARSER FINISHED. EXITING ON ERRORS')
    stop
  else
    bmad_status%ok = .false.
    return
  endif
endif
              
call lat_sanity_check (lat, err_flag)

end subroutine

