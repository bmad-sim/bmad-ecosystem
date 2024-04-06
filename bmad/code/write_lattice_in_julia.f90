!+
! Subroutine write_lattice_in_julia(file_name, lat)
!
! Routine to create a Bmad-Julia lattice file.
!
! Input:
!   lat           -- lat_struct: Lattice
!   file_name     -- character(*): Output lattice file name.
!                     If the name does not have a .jl suffix then this suffix will be added.
!-

subroutine write_lattice_in_julia(file_name, lat)

use write_lat_file_mod, dummy => write_lattice_in_julia

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, lord, slave, slave2, multi_lord
type (multipass_region_lat_struct), target :: mult_lat
type (multipass_all_info_struct), target :: m_info
type (multipass_region_ele_struct), pointer :: mult_ele(:), m_ele
type (ele_pointer_struct), allocatable :: named_eles(:)  ! List of unique element names 

integer n, i, ix, ib, ie, iu, n_names, ix_match, ix_pass, ix_r
integer, allocatable :: an_indexx(:), index_list(:)

logical has_been_added, in_multi_region

character(*) file_name
character(40), allocatable :: names(:)
character(240) fname
character(1000) line
character(*), parameter :: r_name = 'write_lattice_in_julia'

! Open file

call fullfilename(file_name, fname)
call file_suffixer(fname, fname, '.jl', .true.)
iu = lunget()
open (iu, file = fname, status = 'unknown')

! Write element defs

write (iu, '(a)') '@eles begin'

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  ele_loop: do ie = 0, branch%n_ele_max
    ele => branch%ele(ie)
    if (ele%key == overlay$ .or. ele%key == group$ .or. ele%key == ramper$ .or. ele%key == girder$) cycle   ! Not currently handled
    if (ele%key == null_ele$) cycle

    multi_lord => pointer_to_multipass_lord (ele, ix_pass) 
    if (ele%lord_status == super_lord$ .and. ix_pass > 0) cycle
    if (ele%slave_status == super_slave$ .and. ix_pass > 1) cycle

    if (ele%slave_status == super_slave$) then
      lord => pointer_to_lord(ele, 1)
      slave => pointer_to_slave(lord, 1)
      slave2 => pointer_to_slave(lord, lord%n_slave)
      write (iu, '(2(a, i0), 2a)') 'slave_drift_', ib, '_', ele%ix_ele, ': drift, l = ', re_str(ele%value(l$))
      cycle
    endif

    if (ix_pass > 0) cycle

    ! Do not write anything for elements that have a duplicate name.

    call add_this_name_to_list (ele, names, an_indexx, n_names, ix_match, has_been_added, named_eles)
    if (.not. has_been_added) cycle

    ! Write element def

    line = trim(ele%name) // ' = ' // trim(key_name(ele%key)) // '('

    if (ele%type /= ' ') line = trim(line) // ', type = ' // quote(ele%type)
    if (ele%alias /= ' ') line = trim(line) // ', alias = ' // quote(ele%alias)
    if (associated(ele%descrip)) line = trim(line) // ', descrip = ' // quote(ele%descrip)

    !

    if (ele%key == fork$ .or. ele%key == photon_fork$) then
      n = nint(ele%value(ix_to_branch$))
      line = trim(line) // ', to_line = ' // trim(lat%branch(n)%name)
      if (ele%value(ix_to_element$) > 0) then
        i = nint(ele%value(ix_to_element$))
        line = trim(line) // ', to_element = ' // trim(lat%branch(n)%ele(i)%name)
      endif
    endif





    !

    ix = index(line, '(, ')
    line = line(1:ix) // line(ix+3:)

  enddo ele_loop
enddo

write (1, '(a)') 'end'

! Write branch lines
! First write multipass lines

call multipass_region_info(lat, mult_lat, m_info)

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  mult_ele => mult_lat%branch(ib)%ele
  in_multi_region = .false.

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    ix_pass = m_info%branch(ib)%ele(ie)%ix_pass
    if (ix_pass /= 1) cycle 

    if (mult_ele(ie)%region_start_pt) then
      if (in_multi_region) then
        call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #1! PLEASE REPORT THIS!')
      endif
      in_multi_region = .true.
      ix_r = mult_ele(ie)%ix_region
      write (iu, '(a)')
      write (line, '(a, i2.2, a)') 'multi_line_', ix_r, ' = beamline('
    endif

    if (mult_ele(ie)%ix_region /= ix_r) then
      call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #2! PLEASE REPORT THIS!')
    endif

    call write_line_element (line, iu, ele, lat, .true.)

    if (mult_ele(ie)%region_stop_pt) then
      line = line(:len_trim(line)-1) // ')'
      call write_lat_line (line, iu, .true.)
      in_multi_region = .false.
    endif
  enddo

  if (in_multi_region) then
    call out_io (s_error$, r_name, 'MULTIPASS BOOKKEEPING ERROR #3! PLEASE REPORT THIS!')
  endif

enddo  ! ib branch loop


!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  mult_ele => mult_lat%branch(ib)%ele

  in_multi_region = .false.

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    ix_pass = m_info%branch(ib)%ele(ie)%ix_pass
    if (ix_pass /= 1) cycle 
  enddo

  line = trim(branch%name) // ' = beamline(' // quote(branch%name) // ', [' 
  do ie = 0, branch%n_ele_max
    ele => branch%ele(ie)
    line = trim(line) // ' ' // trim(unique_name(ele, lat)) // ','
  enddo
  line = trim(line) // '], geometry = ' // trim(geometry_name(branch%param%geometry)) // ')'
  !!!call write_line(line)
enddo

!----------------------------------------------------------------------------------------------
contains

function unique_name(ele, lat) result (name)

type (lat_struct) lat
type (ele_struct) ele
integer n_match, ine, imax
character(50) name

!

imax = nametable_bracket_indexx(lat%nametable, ele%name, n_match)
if (n_match == 1) then
  name = ele%name
  return
endif

ine = ele_nametable_index(ele)
name = trim(ele%name) // '_' // int_str(ine - (imax - n_match))

end function unique_name

end subroutine
