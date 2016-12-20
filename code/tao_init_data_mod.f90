module tao_init_data_mod

use tao_mod

contains

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!+
! Subroutine tao_init_data (data_file)
!
! Subroutine to initialize the tao data structures.
! If data_file is not in the current directory then it 
! will be searched for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   data_file -- Character(*): Tao data initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_data (data_file)

use tao_data_and_eval_mod
use tao_lattice_calc_mod
use tao_input_struct
use random_mod
use nr, only: indexx
  
implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_input) d2_data
type (tao_d1_data_input) d1_data
type (tao_datum_input) datum(n_data_minn:n_data_maxx) 
type (ele_pointer_struct), allocatable :: eles(:)

real(rp) default_weight, def_weight        ! default merit function weight

integer ios, iu, i, j, j1, k, ix, n_uni, num
integer n, iostat, n_loc
integer n_d1_data, ix_ele, ix_min_data, ix_max_data, ix_d1_data

integer :: n_d2_data(lbound(s%u, 1) : ubound(s%u, 1))

character(*) data_file
character(40) :: r_name = 'tao_init_data'
character(200) file_name
character(40) name,  universe, d_typ, use_same_lat_eles_as
character(40) default_merit_type, default_data_source, def_merit_type, def_data_source
character(100) search_for_lat_eles
character(200) line, default_data_type, def_data_type

logical err, free, gang
logical :: good_unis(lbound(s%u, 1) : ubound(s%u, 1))
logical :: mask(lbound(s%u, 1) : ubound(s%u, 1))

namelist / tao_d2_data / d2_data, n_d1_data, universe, &
                default_merit_type, default_weight, default_data_type, default_data_source

namelist / tao_d1_data / d1_data, datum, ix_d1_data, &
               default_merit_type, default_weight, default_data_type, default_data_source, &
               use_same_lat_eles_as, search_for_lat_eles, ix_min_data, ix_max_data

!-----------------------------------------------------------------------
! Find out how many d2_data structures we need for each universe

call tao_hook_init_data() 
if (.not. s%com%init_data .or. data_file == '') then
  do i = lbound(s%u, 1), ubound(s%u, 1)
    call tao_init_data_in_universe (s%u(i), 0)
  enddo
  call tao_init_data_end_stuff ()
  return
endif

!---

call out_io (s_blank$, r_name, '*Init: Opening Data File: ' // file_name)
call tao_open_file (data_file, iu, file_name, s_fatal$)
if (iu == 0) then
  call out_io (s_fatal$, r_name, 'CANNOT OPEN DATA INIT FILE: ' // data_file)
  call err_exit
endif

n_d2_data = 0

do 
  universe = '*'
  d2_data%name = ''
  read (iu, nml = tao_d2_data, iostat = ios)
  if (ios > 0) then
    call out_io (s_error$, r_name, 'TAO_D2_DATA NAMELIST READ ERROR.')
    rewind (iu)
    do
      read (iu, nml = tao_d2_data)  ! force printing of error message
    enddo
  endif
  if (ios < 0 .and. d2_data%name == '') exit  ! Exit on end-of-file and no namelist read

  if (universe == '*') then
    good_unis = .true.
  else
    call location_decode (universe, good_unis, lbound(s%u, 1), num)
    if (num < 0) then
      call out_io (s_abort$, r_name, &
            'BAD UNIVERSE NUMBER IN TAO_D2_DATA NAMELIST: ' // d2_data%name)
      call err_exit
    endif
  endif

  where (good_unis) n_d2_data = n_d2_data + 1

enddo

rewind (iu)

! Allocate space for the data

do i = lbound(s%u, 1), ubound(s%u, 1)
  call tao_init_data_in_universe (s%u(i), n_d2_data(i))
enddo

! Init data

do 
  mask(:) = .true.      ! set defaults
  d2_data%name           = ''
  universe               = '*'
  default_merit_type     = ''
  default_weight         = 0      
  default_data_type      = ''
  default_data_source    = ''

  read (iu, nml = tao_d2_data, iostat = ios)
  if (ios < 0 .and. d2_data%name == '') exit    ! Exit on end-of-file and no namelist read
  call out_io (s_blank$, r_name, 'Init: Read tao_d2_data namelist: ' // d2_data%name)

  if (universe == '*') then
    good_unis = .true.
  else
    call location_decode (universe, good_unis, lbound(s%u, 1), num)
  endif

  uni_loop: do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. good_unis(i)) cycle  
    ! check if this data type has already been defined for this universe
    do k = 1, size(s%u(i)%d2_data)
      if (trim(s%u(i)%d2_data(k)%name) == trim(d2_data%name)) then
        mask(i) = .false.
        cycle uni_loop
      endif
    enddo
      
    call tao_d2_data_stuffit (s%u(i), d2_data%name, n_d1_data)
  enddo uni_loop

  def_merit_type  = default_merit_type   ! Save
  def_weight      = default_weight
  def_data_type   = default_data_type
  def_data_source = default_data_source

  do k = 1, n_d1_data
    use_same_lat_eles_as   = ''
    search_for_lat_eles    = ''
    d1_data%name           = ''
    default_merit_type     = def_merit_type
    default_weight         = def_weight
    default_data_type      = def_data_type
    default_data_source    = def_data_source
    ix_min_data            = int_garbage$
    ix_max_data            = int_garbage$

    datum(:) = tao_datum_input()

    ! Read datum/data

    read (iu, nml = tao_d1_data, iostat = ios)
    if (ios > 0) then
      call out_io (s_error$, r_name, 'TAO_D1_DATA NAMELIST READ ERROR.')
      rewind (iu)
      do
        read (iu, nml = tao_d1_data)  ! force printing of error message
      enddo
    endif

    ! Convert old format to new

    if (datum(0)%ele_name(1:7) == 'SEARCH:') then
      call string_trim(datum(0)%ele_name(8:), search_for_lat_eles, ix)
    elseif (datum(0)%ele_name(1:5) == 'SAME:') then
      call string_trim (datum(0)%ele_name(6:), use_same_lat_eles_as, ix)
    endif

    ! Check that we read the correct namelist

    if (ix_d1_data /= k) then
      write (line, '(a, 2i4)') ', k, ix_d1_data'
      call out_io (s_abort$, r_name, &
                'ERROR: IX_D1_DATA MISMATCH FOR D2_DATA: ' // d2_data%name, &
                '       THE D1_DATA HAD THE NAME: ' // d1_data%name, &
                '       I EXPECTED IX_D1_DATA TO BE: \i3\ ', &
                '       I READ IX_D1_DATA TO BE: \i3\ ', &
                i_array = (/ k, ix_d1_data /) )  
      call err_exit
    endif
    do i = lbound(datum, 1), ubound(datum, 1)
      if ((datum(i)%ele_ref_name /= '' .or. datum(i)%ele_start_name /= '') .and. datum(i)%ele_name == '') then
        write (line, '(4a, i0, a)') trim(d2_data%name), '.', trim(d1_data%name), '[', i, ']'
        call out_io (s_abort$, r_name, &
              'ERROR: ELE_NAME IS BLANK BUT ELE_REF_NAME OR ELE_START_NAME IS NOT FOR: ' // line)
        call err_exit
      endif

      if (index(datum(i)%data_type, 'dat::') /= 0) then
        call out_io (s_error$, r_name, &
                     'DATA_TYPE USES OLD "dat::" PREFIX. PLEASE CHANGE TO "data::": ' // datum(i)%data_type)
        call err_exit
      endif
    enddo
    call out_io (s_blank$, r_name, 'Init: Read tao_d1_data namelist: ' // d1_data%name)

    do i = lbound(s%u, 1), ubound(s%u, 1)
      if (.not. good_unis(i)) cycle
      if (.not. mask(i)) cycle
      call d1_data_stuffit (k, s%u(i), s%u(i)%n_d2_data_used)
    enddo

  enddo

enddo

close (iu)

! Init ix_data array

call tao_init_data_end_stuff ()

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
contains

subroutine d1_data_stuffit (i_d1, u, n_d2)

type (tao_universe_struct), target :: u
type (tao_d1_data_struct), pointer :: d1_this
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)
type (ele_pointer_struct), allocatable, save :: eles(:)
type (ele_struct), pointer :: ele
type (tao_data_struct), pointer :: dat

real(rp), allocatable :: s(:)

integer i, n1, n2, ix, k, ix1, ix2, j, jj, n_d2, i_d1
integer, allocatable :: ind(:)

character(20) fmt

!

d1_this => u%d2_data(n_d2)%d1(i_d1)  
if (d1_data%name == '') then
  write (d1_this%name, '(i0)') i_d1
else
  d1_this%name = d1_data%name    ! stuff in the data
endif

!-----------------------------------------
! Check if we are searching for elements or repeating elements
! and record the element names in the data structs.
    
if (search_for_lat_eles /= '') then
  call tao_init_find_elements (u, search_for_lat_eles, eles)
  if (size(eles) == 0) then
    call out_io (s_warn$, r_name, &
      'NO ELEMENTS FOUND IN SEARCH FOR: ' // search_for_lat_eles, &
      'WHILE SETTING UP DATA ARRAY: ' // tao_d2_d1_name(d1_this))
    return
  endif
  ! finish finding data array limits
  n1 = u%n_data_used + 1
  n2 = u%n_data_used + size(eles)
  call tao_allocate_data_array (u, n2)

  if (ix_min_data == int_garbage$) ix_min_data = 1
  ix1 = ix_min_data
  ix2 = ix1 + (n2 - n1)

  allocate (ind(size(eles)), s(size(eles)))
  do i = 1, size(eles)
    s(i) = eles(i)%ele%s
  enddo
  call indexx (s, ind)

  ! get element names

  jj = n1
  do k = lbound(eles, 1), ubound(eles, 1)
    if (jj .gt. n2) then
      call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
      call err_exit
    endif
    u%data(jj)%ele_name  = eles(ind(k))%ele%name
    u%data(jj)%ix_ele    = eles(ind(k))%ele%ix_ele
    u%data(jj)%ix_branch = eles(ind(k))%ele%ix_branch
    u%data(jj)%exists    = .true.
    jj = jj + 1
  enddo

  u%data(n1:n2)%good_user      = datum(ix1)%good_user
  u%data(n1:n2)%invalid_value  = datum(ix1)%invalid_value
  u%data(n1:n2)%ele_start_name = datum(ix1)%ele_start_name
  u%data(n1:n2)%ele_ref_name   = datum(ix1)%ele_ref_name
  u%data(n1:n2)%ix_bunch       = datum(ix1)%ix_bunch
  u%data(n1:n2)%invalid_value  = datum(ix1)%invalid_value
  u%data(n1:n2)%data_type      = datum(ix1)%data_type
  u%data(n1:n2)%merit_type     = datum(ix1)%merit_type
  u%data(n1:n2)%weight         = datum(ix1)%weight
  u%data(n1:n2)%s_offset       = datum(ix1)%s_offset
  u%data(n1:n2)%data_source    = datum(ix1)%data_source
  u%data(n1:n2)%meas_value     = 0  
  ! use default_data_type if given, if not, auto-generate the data_type
  if (default_data_type == '') default_data_type = trim(d2_data%name) // '.' // d1_data%name
  where (u%data(n1:n2)%data_type == '') u%data(n1:n2)%data_type = default_data_type
  ! eval_point
  call match_word (datum(ix1)%eval_point, anchor_pt_name(1:), ix, can_abbreviate = .false.)
  if (ix == 0) then
    call out_io (s_abort$, r_name, 'EVAL_POINT UNRECOGNIZED: ' // datum(ix1)%eval_point)
    stop
  endif
  u%data(n1:n2)%eval_point = ix

  deallocate (ind)

!-----------------------------------------
! use_same_lat_eles_as

elseif (use_same_lat_eles_as /= '') then
  call string_trim (use_same_lat_eles_as, name, ix)
  call tao_find_data (err, name, d1_array = d1_array, ix_uni = u%ix_uni)
  if (err .or. size(d1_array) /= 1) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
    call err_exit
  endif
  n1 = u%n_data_used + 1
  n2 = u%n_data_used + size(d1_array(1)%d1%d)
  call tao_allocate_data_array (u, n2)

  ix_min_data = lbound(d1_array(1)%d1%d, 1)
  ix1 = ix_min_data
  ix2 = ix1 + (n2 - n1)

  u%data(n1:n2)%ele_name        = d1_array(1)%d1%d%ele_name
  u%data(n1:n2)%ix_ele          = d1_array(1)%d1%d%ix_ele
  u%data(n1:n2)%ele_ref_name    = d1_array(1)%d1%d%ele_ref_name
  u%data(n1:n2)%ix_ele_ref      = d1_array(1)%d1%d%ix_ele_ref
  u%data(n1:n2)%ele_start_name  = d1_array(1)%d1%d%ele_start_name
  u%data(n1:n2)%ix_ele_start    = d1_array(1)%d1%d%ix_ele_start
  u%data(n1:n2)%exists          = d1_array(1)%d1%d%exists
  u%data(n1:n2)%s_offset        = d1_array(1)%d1%d%s_offset
  u%data(n1:n2)%eval_point      = d1_array(1)%d1%d%eval_point

  u%data(n1:n2)%invalid_value = datum(ix1)%invalid_value
  u%data(n1:n2)%meas_value    = 0  

  if (default_data_type /= '')    u%data(n1:n2)%data_type = default_data_type
  if (datum(ix1)%data_type /= '') u%data(n1:n2)%data_type = datum(ix1)%data_type

  if (default_data_source /= '')    u%data(n1:n2)%data_source = default_data_source
  if (datum(ix1)%data_source /= '') u%data(n1:n2)%data_source = datum(ix1)%data_source

  if (default_merit_type /= '')    u%data(n1:n2)%merit_type = default_merit_type
  if (datum(ix1)%merit_type /= '') u%data(n1:n2)%merit_type = datum(ix1)%merit_type

  if (default_weight /= 0)    u%data(n1:n2)%weight = default_weight
  if (datum(ix1)%weight /= 0) u%data(n1:n2)%weight = datum(ix1)%weight

  ! use default_data_type if given, if not, auto-generate the data_type
  if (default_data_type == '') default_data_type = trim(d2_data%name) // '.' // d1_data%name
  where (u%data(n1:n2)%data_type == '') u%data(n1:n2)%data_type = default_data_type

!-----------------------------------------
! Not SEARCH or SAME:

else

  if (ix_min_data == int_garbage$) ix_min_data = 1
  if (ix_max_data == int_garbage$) then
    do i = ubound(datum, 1), lbound(datum, 1), -1
      if (datum(i)%ele_name /= '' .or. datum(i)%data_type /= '') then
        ix_max_data = i
        exit
      endif
    enddo
  endif

  if (ix_max_data == int_garbage$) then
    call out_io (s_error$, r_name, 'NO DATA FOUND FOR: ' // tao_d2_d1_name(d1_this))
    return
  endif

  n1 = u%n_data_used + 1
  n2 = u%n_data_used + ix_max_data - ix_min_data + 1
  ix1 = ix_min_data
  ix2 = ix_max_data
  call tao_allocate_data_array (u, n2)

  ! Transfer info from the input structure

  u%data(n1:n2)%good_user      = datum(ix1:ix2)%good_user
  u%data(n1:n2)%weight         = datum(ix1:ix2)%weight
  u%data(n1:n2)%ele_name       = datum(ix1:ix2)%ele_name
  u%data(n1:n2)%invalid_value  = datum(ix1:ix2)%invalid_value
  u%data(n1:n2)%ele_ref_name   = datum(ix1:ix2)%ele_ref_name
  u%data(n1:n2)%ele_start_name = datum(ix1:ix2)%ele_start_name
  u%data(n1:n2)%ix_bunch       = datum(ix1:ix2)%ix_bunch
  u%data(n1:n2)%data_source    = datum(ix1:ix2)%data_source
  u%data(n1:n2)%data_type      = datum(ix1:ix2)%data_type
  u%data(n1:n2)%merit_type     = datum(ix1:ix2)%merit_type
  u%data(n1:n2)%invalid_value  = datum(ix1:ix2)%invalid_value
  u%data(n1:n2)%s_offset       = datum(ix1:ix2)%s_offset

  u%data(n1:n2)%meas_value = datum(ix1:ix2)%meas
  where (u%data(n1:n2)%meas_value == real_garbage$)  ! where %meas_value was set
    u%data(n1:n2)%meas_value = 0  
  elsewhere
    u%data(n1:n2)%good_meas = .true.
  end where

  ! use default_data_type if given, if not, auto-generate the data_type
  if (default_data_type == '') default_data_type = trim(d2_data%name) // '.' // d1_data%name
  where (u%data(n1:n2)%data_type == '') u%data(n1:n2)%data_type = default_data_type

  ! Find elements associated with the data

  do j = n1, n2

    call match_word (datum(j+ix1-n1)%eval_point, anchor_pt_name(1:), ix, can_abbreviate = .false.)
    if (ix == 0) then
      call out_io (s_abort$, r_name, 'EVAL_POINT UNRECOGNIZED: ' // datum(j+ix1-n1)%eval_point)
      stop
    endif
    u%data(j)%eval_point = ix

    if (data_type_needs_associated_element_to_exist(u%data(j)%data_type) .and. &
                                                         u%data(j)%ele_name == '') cycle

    u%data(j)%exists = .true.

    if (u%data(j)%ele_name /= '') then
      call str_upcase (u%data(j)%ele_name, u%data(j)%ele_name)
      call lat_ele_locator (u%data(j)%ele_name, u%design%lat, eles, n_loc)
      if (n_loc == 0) then
        call out_io (s_error$, r_name, 'ELEMENT NOT LOCATED: ' // u%data(j)%ele_name)
        u%data(j)%exists = .false.
        cycle
      endif
      u%data(j)%ix_ele    = eles(1)%ele%ix_ele
      u%data(j)%ix_branch = eles(1)%ele%ix_branch
    endif

    if (u%data(j)%ele_ref_name /= '') then
      call str_upcase (u%data(j)%ele_ref_name, u%data(j)%ele_ref_name)
      call lat_ele_locator (u%data(j)%ele_ref_name, u%design%lat, eles, n_loc)
      if (n_loc == 0) then
        call out_io (s_error$, r_name, 'ELE_REF NOT LOCATED: ' // u%data(j)%ele_ref_name)
        u%data(j)%exists = .false.
        cycle
      endif
      u%data(j)%ix_ele_ref = eles(1)%ele%ix_ele
    endif

    if (u%data(j)%ele_start_name /= '') then
      call str_upcase (u%data(j)%ele_start_name, u%data(j)%ele_start_name)
      call lat_ele_locator (u%data(j)%ele_start_name, u%design%lat, eles, n_loc)
      if (n_loc == 0) then
        call out_io (s_error$, r_name, 'ELE_START NOT LOCATED: ' // u%data(j)%ele_start_name)
        u%data(j)%exists = .false.
        cycle
      endif
      u%data(j)%ix_ele_start = eles(1)%ele%ix_ele
    endif

  enddo

endif

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! Point the %data back to the d1_data_struct

call tao_point_d1_to_data (d1_this, u%data(n1:n2), ix_min_data)

! In a d1_data array, not all the datums need to exist. 
! If a datum is not associated with an element, that generally means that
! it does not exist. There are, however, a few exceptions. EG: unstable.ring, etc.
! Here we mark data%exists for such datums.
! Also determine if we need to do the radiation integrals. This can save a lot of time.

do j = n1, n2
  dat => u%data(j)

  ! Use defaults if a component has not been set.

  if (dat%weight == 0)       dat%weight = default_weight
  if (dat%merit_type == '')  dat%merit_type = default_merit_type
  if (dat%merit_type == '')  dat%merit_type = 'target'
  if (dat%data_source == '') dat%data_source = default_data_source
  if (dat%data_source == '') dat%data_source = 'lat'

  ! Convert old style to new style

  ix = index(dat%data_type, 'emittance.')
  if (ix /= 0) dat%data_type = dat%data_type(1:ix-1) // 'emit.' // dat%data_type(ix+10:)
  if (dat%data_type(1:9) == 'unstable_') dat%data_type(9:9) = '.'

  !

  if (tao_rad_int_calc_needed(dat%data_type, dat%data_source)) then
    u%calc%rad_int_for_data = .true. 
    if (dat%ix_branch /= 0) then
      call out_io (s_fatal$, r_name, 'EVALUATING A DATUM OF TYPE: ' // dat%data_type, 'ON A BRANCH NOT YET IMPLEMENTED!')
      call err_exit
    endif
  endif

  if (tao_chrom_calc_needed(dat%data_type, dat%data_source)) then
    if (u%model%lat%branch(dat%ix_branch)%param%geometry == open$) then
      call out_io (s_warn$, r_name, 'CHROMATICITY DATUM NOT VALID FOR NON-CLOSED LATTICE!')
      dat%exists = .false.
    else
      u%calc%chrom_for_data = .true.
    endif  
  endif

  if (tao_beam_sigma_calc_needed(dat%data_type, dat%data_source)) then
    u%calc%beam_sigma_for_data = .true. 
  endif

  ! Some data types are global and are not associated with a particular element. Check for this.

  if (dat%data_type == 'unstable.orbit') then
    dat%exists = .true.
    if (dat%ele_name /= '') then
      call out_io (s_abort$, r_name, 'DATUM OF TYPE: ' // dat%data_type, &
                        'CANNOT HAVE AN ASSOCIATED ELEMENT: ' // dat%ele_name)
      call err_exit
    endif
  endif

  if (u%design%lat%param%geometry == closed$ .and. &
              (dat%data_type(1:12)  == 'chrom.dtune.' .or. dat%data_type(1:5)  == 'damp.' .or. &
               dat%data_type(1:17) == 'multi_turn_orbit.' .or. dat%data_type(1:5) == 'tune.' .or. &
               dat%data_type(1:13) == 'unstable.ring' .or. index(dat%data_type, 'emit.') /= 0)) then
    dat%exists = .true.
    if (dat%ele_name /= '') then
      call out_io (s_abort$, r_name, 'DATUM OF TYPE: ' // dat%data_type, &
                        'CANNOT HAVE AN ASSOCIATED ELEMENT IN A CIRCULAR LATTICE: ' // dat%ele_name)
      call err_exit
    endif
  endif

enddo

if (.not. any(u%data(n1:n2)%exists)) then
  call out_io (s_warn$, r_name, &
            'Note: All datums in: ' // tao_d2_d1_name(d1_this), &
            'are marked as non-existent')
endif

end subroutine d1_data_stuffit

end subroutine tao_init_data

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Defines what datums to evaluate at each element in specified universe

subroutine tao_init_data_end_stuff ()

implicit none

type (tao_universe_struct), pointer :: u
type (tao_data_struct), pointer :: data
type (tao_element_struct), pointer :: uni_ele(:)
integer i, ib, j, k, ix_ele, n_max

logical err

!

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  call tao_allocate_data_array (u, u%n_data_used, .true.) ! Trim u%data size
enddo

call tao_data_check (err)
if (err) stop

!----------------------------------------------------------------------
contains

function choose_ix_ele() result (ix_ele)
integer ix_ele

if (.not. data%exists) then
  ix_ele = int_garbage$
elseif (data%data_type(1:17) == 'multi_turn_orbit.') then
  ix_ele = int_garbage$ ! Does not get evaluated by tao_lattice_calc_mod
elseif (data%data_source /= 'beam') then
  ix_ele = -1
elseif (data%data_type(1:7) == 'rad_int') then
  ix_ele = -1
elseif (data%ix_ele > s%u(data%d1%d2%ix_uni)%model%lat%n_ele_track) then
  ix_ele = -1
elseif (data%ix_ele == -1) then
  ix_ele = -1
elseif (index(data%data_type, 'emit.') /= 0 .and. data%data_source == 'lat') then
  ix_ele = -1
elseif (data%data_type(1:6) == 'chrom.') then
  ix_ele = -1
elseif (data%ix_ele_ref > data%ix_ele) then
  ix_ele = u%model%lat%n_ele_track
else
  ix_ele = data%ix_ele
endif

end function choose_ix_ele

end subroutine tao_init_data_end_stuff

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_allocate_data_array (u, n_data, exact)

type (tao_universe_struct) :: u
type (tao_data_struct), allocatable :: data(:)
type (tao_d1_data_struct), pointer :: d1

integer i, j1, j2, n0, n_data
logical, optional :: exact  ! Default = False

! Exact means that size(u%data) must end up to be n_data.
! Not exact means that size(u%data) must be at least n_data.

u%n_data_used = n_data
  
if (n_data == size(u%data)) return 
if (.not. logic_option(.false., exact) .and. n_data < size(u%data)) return 

! Reallocate the data array. 
! If not exact then allocate more space than needed to reduce the number
! of times we need to reallocate stuff.

if (allocated(u%data)) then
  n0 = min(n_data, size(u%data))
  allocate (data(n0))
  data = u%data(1:n0)
  deallocate (u%data)
  if (logic_option(.false., exact)) then
    allocate (u%data(n_data))
  else
    allocate (u%data(2*n_data))
  endif
  u%data(1:n0) = data
  deallocate (data)
else
  allocate(u%data(n_data))
endif

! Since the data array gets reallocated the pointer from d1 to the datums must 
! be reestablished.

j2 = 0
do
  j1 = j2 + 1
  if (j1 > n0) exit
  d1 => u%data(j1)%d1
  if (.not. associated(d1)) exit
  do 
    if (j2 == n0) exit
    if (.not. associated(u%data(j2+1)%d1, d1)) exit
    j2 = j2 + 1
  enddo
  call tao_point_d1_to_data (d1, u%data(j1:j2), u%data(j1)%ix_d1)
enddo

! Set %ix_data. See the tao_data_struct for the defaults component values.

do i = n0+1, size(u%data)
  u%data(i)%ix_data = i
enddo

end subroutine tao_allocate_data_array

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_d2_data_stuffit (u, d2_name, n_d1_data)

type (tao_universe_struct), target :: u
type (tao_d2_data_struct), pointer :: d2

integer i, nn, n_d1_data
character(*) d2_name
character(40) :: r_name = 'tao_d2_data_stuffit'

! Setup another d2_data structure.

u%n_d2_data_used = u%n_d2_data_used + 1
nn = u%n_d2_data_used

if (size(u%d2_data) < nn) then
  call out_io (s_error$, r_name, 'D2_DATA ARRAY OVERFLOW!')
  call err_exit
endif

d2 => u%d2_data(nn)

d2%name = d2_name
d2%ix_uni = u%ix_uni

! allocate memory for the u%d1_data structures

if (allocated(d2%d1)) deallocate (d2%d1)
allocate(d2%d1(n_d1_data))

do i = 1, n_d1_data
  d2%d1(i)%d2 => d2
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine tao_init_data_in_universe (u, n_d2_data)

implicit none

type (tao_universe_struct) u
integer j, n_d2_data

!

allocate (u%data(0))
u%n_d2_data_used = 0      ! size of s%u(i)%d2_data(:) array
u%n_data_used = 0         ! size of s%u(i)%data(:) array
u%model%ix_rad_int_cache = 0
u%design%ix_rad_int_cache = 0
u%base%ix_rad_int_cache = 0

if (n_d2_data == 0) return
if (allocated(u%d2_data)) deallocate (u%d2_data)

allocate (u%d2_data(n_d2_data))

do j = 1, n_d2_data
  u%d2_data(j)%ix_d2_data = j
  u%d2_data(j)%descrip = ''
enddo

u%d2_data%name = ''  ! blank name means it doesn't exist

! This is needed to keep the totalview debugger happy.

if (allocated(u%dmodel_dvar)) deallocate (u%dmodel_dvar)
allocate (u%dmodel_dvar(1,1))

end subroutine tao_init_data_in_universe

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_point_d1_to_data (d1, data, n_min)
!
! Routine used for arbitrary data pointer indexing
!
! d1     -- tao_data_struct: the pointer
! data   -- tao_data_struct: the data
! n_min  -- integer: starting index for the pointer
!-

subroutine tao_point_d1_to_data (d1, data, n_min)

implicit none

integer n, n_min, i, n0, n1

type (tao_d1_data_struct), target :: d1
type (tao_data_struct), target :: data(n_min:)

d1%d => data

do n = lbound(data, 1), ubound(data, 1)
  data(n)%d1 => d1
  data(n)%ix_d1 = n
enddo

end subroutine tao_point_d1_to_data

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function data_type_needs_associated_element_to_exist (data_type) result (needs_element)
!
! Routne to determine whether a given data type needs to have an associated element to exist.
! For example, "phase.X" type data needs an associated element to exist.
!
! Input:
!   data_type -- Character(*): Type of data
!
! Output:
!   needs_element -- Logical: Set true if data type needs an associated element.
!-

function data_type_needs_associated_element_to_exist (data_type) result (needs_element)

character(*) data_type
character(40) head

logical needs_element
integer ix

!

head = data_type
ix = index(head, '.')
if (ix /= 0) head = head(1:ix)

needs_element = .true.

select case (head)
case ('alpha.')
case ('apparent_emit.', 'norm_apparent_emit.')
case ('beta.')
case ('bpm_cbar.')
case ('bpm_eta.')
case ('bpm_k.')
case ('bpm_orbit.')
case ('bpm_phase.')
case ('cbar.')
case ('chrom.')
case ('damp.')
case ('dpx_dx') 
case ('dpy_dy') 
case ('dpz_dz') 
case ('e_tot')
case ('element_attrib.')
case ('emit.', 'norm_emit.')
case ('eta.')
case ('etap.')
case ('expression:', 'expression.')
case ('floor.')
case ('gamma.')
case ('k.')
case ('momentum')
case ('momentum_compaction')
case ('n_particle_loss')
case ('normal.')
case ('orbit.')
case ('periodic.')
case ('phase.', 'phase_frac.')
case ('phase_frac_diff')
case ('photon.')
case ('r.')
case ('rad_int.')
case ('rad_int1.')
case ('ref_time')
case ('rel_floor.')
case ('s_position') 
case ('sigma.')
case ('spin.')
case ('t.', 'tt.')
case ('time')
case ('tune.')
case ('unstable.')
case ('wall.')
case ('wire.')  
case ('')
  needs_element = .false.
case default
  needs_element = .false.
end select

end function data_type_needs_associated_element_to_exist

end module
