module tao_init_data_mod

use tao_interface

integer, parameter, private :: n_datum_min     = -100   ! min index of datum per d1_data

contains

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!+
! Subroutine tao_init_data (data_file)
!
! Subroutine to initialize the tao data structures.
!
! Input:
!   data_file -- Character(*): Tao data initialization file.
!                  If blank, there is no file so just use the defaults.
!-

subroutine tao_init_data (data_file)

use tao_data_and_eval_mod
use tao_input_struct
use random_mod
  
implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_input) d2_data
type (tao_d1_data_input) d1_data
type (tao_datum_input), allocatable :: datum(:)
type (lat_struct), pointer :: lat

real(rp) default_weight, def_weight        ! default merit function weight

integer ios, iu, i, j, j1, k, ix, n_uni, num
integer n, iostat, n_loc, n_nml
integer n_d1_data, ix_ele, ix_min_data, ix_max_data, ix_d1_data

integer :: n_d2_data(lbound(s%u, 1):ubound(s%u, 1)), n_data(lbound(s%u, 1):ubound(s%u, 1))

character(*) data_file
character(40) :: r_name = 'tao_init_data'
character(200) file_name
character(40) name,  universe, d_typ, use_same_lat_eles_as
character(40) default_merit_type, default_data_source, def_merit_type, def_data_source
character(200) search_for_lat_eles
character(200) line, default_data_type, def_data_type

logical err, free, gang, found
logical :: good_uni(lbound(s%u, 1) : ubound(s%u, 1))

namelist / tao_d2_data / d2_data, n_d1_data, universe, &
                default_merit_type, default_weight, default_data_type, default_data_source

namelist / tao_d1_data / d1_data, datum, ix_d1_data, &
               default_merit_type, default_weight, default_data_type, default_data_source, &
               use_same_lat_eles_as, search_for_lat_eles, ix_min_data, ix_max_data

!-----------------------------------------------------------------------

n = 4000
do i = lbound(s%u, 1), ubound(s%u,1)
  lat => s%u(i)%model%lat
  do j = 0, ubound(lat%branch, 1)
    n = max(n, lat%branch(j)%n_ele_max+10)
    s%u(i)%model%tao_branch(j)%n_hterms = 0
  enddo
enddo

allocate (datum(n_datum_min:n))

! Find out how many d2_data structures we need for each universe

if (data_file == '') then
  do i = lbound(s%u, 1), ubound(s%u, 1)
    call tao_init_data_in_universe (s%u(i), 0)
  enddo
  call tao_hook_init_data()
  call tao_init_data_end_stuff()
  return
endif

!---

call out_io (s_blank$, r_name, '*Init: Opening Data File: ' // data_file)
call tao_open_file (data_file, iu, file_name, s_fatal$)
if (iu == 0) then
  call out_io (s_fatal$, r_name, 'CANNOT OPEN DATA INIT FILE: ' // data_file)
  return
endif

n_d2_data = 0
n_nml = 0

do 
  universe = '*'
  d2_data%name = ''
  n_nml = n_nml + 1
  read (iu, nml = tao_d2_data, iostat = ios)
  if (ios > 0) then
    call out_io (s_error$, r_name, 'TAO_D2_DATA NAMELIST READ ERROR IN FILE: ' // data_file, &
                                   'THIS IS THE ' // ordinal_str(n_nml) // ' TAO_D2_DATA NAMELIST IN THE FILE', &
                                   'WITH D2_DATA%NAME = ' // quote(d2_data%name))
    rewind (iu)
    do
      read (iu, nml = tao_d2_data)  ! force printing of error message
    enddo
  endif
  if (ios < 0 .and. d2_data%name == '') exit  ! Exit on end-of-file and no namelist read

  if (.not. tao_is_valid_name(d2_data%name, line)) then
    call out_io (s_error$, r_name, 'D2_DATA%NAME IN TAO_D2_DATA NAMELIST IS INVALID SINCE: ' // line, &
                                   'IN FILE: ' // data_file)
    cycle
  endif

  if (universe == '*') then
    good_uni = .true.
  else
    call location_decode (universe, good_uni, lbound(s%u, 1), num)
    if (num < 0) then
      call out_io (s_abort$, r_name, &
            'BAD UNIVERSE NUMBER IN TAO_D2_DATA NAMELIST: ' // quote(d2_data%name), &
            'IN FILE: ' // data_file)
      cycle
    endif
  endif

  where (good_uni) n_d2_data = n_d2_data + 1
enddo

rewind (iu)

! Allocate space for the data

do i = lbound(s%u, 1), ubound(s%u, 1)
  call tao_init_data_in_universe (s%u(i), n_d2_data(i))
enddo

!--------------------------------------------
! Now fill in the data.
! Loop over d2_data namelists.

do 
  d2_data%name           = ''
  universe               = '*'
  default_merit_type     = ''
  default_weight         = 0
  default_data_type      = ''
  default_data_source    = ''

  read (iu, nml = tao_d2_data, iostat = ios)
  if (ios < 0 .and. d2_data%name == '') exit    ! Exit on end-of-file and no namelist read
  call out_io (s_blank$, r_name, 'Init: Read tao_d2_data namelist: ' // quote(d2_data%name))

  if (universe == '*') then
    good_uni = .true.
  else
    call location_decode (universe, good_uni, lbound(s%u, 1), num)
  endif

  uni_loop: do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. good_uni(i)) cycle  
    ! check if this data type has already been defined for this universe
    do k = 1, size(s%u(i)%d2_data)
      if (trim(s%u(i)%d2_data(k)%name) /= trim(d2_data%name)) cycle
      good_uni(i) = .false.
      call out_io (s_error$, r_name, 'TWO D2 DATA STRUCTURES HAVE THE SAME NAME: ' // quote(d2_data%name), &
                                     'THE SECOND ONE WILL BE IGNORED!')
      cycle uni_loop
    enddo

    call tao_d2_data_stuffit (s%u(i), d2_data%name, n_d1_data)
    n_data(i) = s%u(i)%n_data_used
  enddo uni_loop

  def_merit_type  = default_merit_type   ! Save
  def_weight      = default_weight
  def_data_type   = default_data_type
  def_data_source = default_data_source

  !----------------------------
  ! Loop over d1_data namelists

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
    if (ios /= 0) then
      if (ios < 0) then
        call out_io (s_error$, r_name, 'TAO_D1_DATA NAMELIST END-OF-FILE READ ERROR (MISSING/BAD QUOTATION MARK PERHAPS?).', &
                                       'IN FILE: ' // data_file, &
                                       'THIS IS THE ' // ordinal_str(k) // ' TAO_D1_DATA NAMELIST AFTER THE ', &
                                       'TAO_D2_DATA NAMELIST WITH D2_DATA%NAME = ' // quote(d2_data%name))
      else
        call out_io (s_error$, r_name, 'TAO_D1_DATA NAMELIST READ ERROR.', &
                                       'IN FILE: ' // data_file, &
                                       'THIS IS THE ' // ordinal_str(k) // ' TAO_D1_DATA NAMELIST AFTER THE ', &
                                       'TAO_D2_DATA NAMELIST WITH D2_DATA%NAME = ' // quote(d2_data%name))
      endif
      rewind (iu)
      do
        read (iu, nml = tao_d1_data)  ! Force printing of error message
      enddo
    endif

    if (.not. tao_is_valid_name(d2_data%name, line)) then
      call out_io (s_error$, r_name, 'D1_DATA%NAME IN TAO_D1_DATA NAMELIST IS INVALID SINCE: ' // line, &
                                     'IN FILE: ' // data_file)
      cycle
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
      call out_io (s_error$, r_name, &
                'IX_D1_DATA MISMATCH FOR D2_DATA: ' // quote(d2_data%name), &
                'THE D1_DATA HAD THE NAME: ' // quote(d1_data%name), &
                'I EXPECTED IX_D1_DATA TO BE: \i3\ ', &
                'I READ IX_D1_DATA TO BE: \i3\ ', &
                'THIS D2_DATA WILL BE IGNORED!', i_array = [k, ix_d1_data])
      do i = lbound(s%u, 1), ubound(s%u, 1)
        if (.not. good_uni(i)) cycle
        s%u(i)%n_d2_data_used = s%u(i)%n_d2_data_used - 1
        s%u(i)%n_data_used = n_data(i)
      enddo
      exit
    endif

    do i = lbound(datum, 1), ubound(datum, 1)
      if (index(datum(i)%data_type, 'dat::') /= 0) then
        call out_io (s_error$, r_name, 'DATA_TYPE USES OLD "dat::" PREFIX. PLEASE CHANGE TO "data::": ' // quote(datum(i)%data_type))
        datum(i)%data_type = 'data::' // datum(i)%data_type(6:)
      endif
    enddo
    call out_io (s_blank$, r_name, 'Init: Read tao_d1_data namelist: ' // quote(d1_data%name))

    do i = lbound(s%u, 1), ubound(s%u, 1)
      if (.not. good_uni(i)) cycle
      call d1_data_stuffit (k, s%u(i), s%u(i)%n_d2_data_used, datum)
    enddo

  enddo  ! d1_data loop
enddo  ! d2_data loop

close (iu)

! Custom data setup?

call tao_hook_init_data()

! Init ix_data array

call tao_init_data_end_stuff ()

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
contains

subroutine d1_data_stuffit (i_d1, u, n_d2, datum)

use srdt_mod

type (tao_universe_struct), target :: u
type (tao_datum_input) :: datum(n_datum_min:)
type (tao_d1_data_struct), pointer :: d1_this
type (tao_d1_data_array_struct), allocatable :: d1_array(:)
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele
type (tao_data_struct), pointer :: dat
type (resonance_h_struct), allocatable :: h_temp(:)
type (bmad_normal_form_struct), pointer :: norm_form
real(rp), allocatable :: s(:)

integer i, n1, n2, ix, k, ix1, ix2, j, jj, n_d2, i_d1, has_associated_ele
integer, allocatable :: ind(:)
integer, pointer :: n_ht

character(20) fmt
character(6) h_str

!

d1_this => u%d2_data(n_d2)%d1(i_d1)
if (d1_data%name == '') then
  write (d1_this%name, '(i0)') i_d1
else
  d1_this%name = d1_data%name    ! stuff in the data
endif

! If given, use the default_data_type. If not, auto-generate the data_type.

if (default_data_type == '') default_data_type = trim(d2_data%name) // '.' // d1_data%name

!-----------------------------------------
! Check if we are searching for elements or repeating elements
! and record the element names in the data structs.
    
if (search_for_lat_eles /= '') then
  call tao_init_find_elements (u, search_for_lat_eles, eles)
  if (size(eles) == 0) then
    call out_io (s_warn$, r_name, &
      'NO ELEMENTS FOUND IN SEARCH FOR: ' // search_for_lat_eles, &
      'WHILE SETTING UP DATA ARRAY: ' // tao_d2_d1_name(d1_this))
      allocate (u%d2_data(n_d2)%d1(i_d1)%d(0))
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
  call indexer (s, ind)

  ! Get element names

  u%data(n1:n2)%good_user        = datum(ix1:ix2)%good_user
  u%data(n1:n2)%good_opt         = datum(ix1:ix2)%good_opt
  u%data(n1:n2)%invalid_value    = datum(ix1:ix2)%invalid_value
  u%data(n1:n2)%spin_map%axis_input   = datum(ix1:ix2)%spin_axis
  u%data(n1:n2)%ele_start_name   = datum(ix1:ix2)%ele_start_name
  u%data(n1:n2)%ele_ref_name     = datum(ix1:ix2)%ele_ref_name
  u%data(n1:n2)%ix_bunch         = datum(ix1:ix2)%ix_bunch
  u%data(n1:n2)%merit_type       = datum(ix1:ix2)%merit_type
  u%data(n1:n2)%weight           = datum(ix1:ix2)%weight
  u%data(n1:n2)%s_offset         = datum(ix1:ix2)%s_offset
  u%data(n1:n2)%data_source      = datum(ix1:ix2)%data_source
  u%data(n1:n2)%meas_value       = datum(ix1:ix2)%meas
  u%data(n1:n2)%error_rms        = datum(ix1:ix2)%error_rms

  do i = n1, n2
    j = i - n1 + ix1
    u%data(i)%data_type = trim(datum(j)%data_type)
    if (u%data(i)%data_type == '') u%data(i)%data_type = trim(default_data_type)
  enddo

  !

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

    call find_this_element (u%data(jj)%ele_ref_name,   'ELE_REF ELEMENT',   u, u%data(jj), u%data(jj)%ix_ele_ref)
    call find_this_element (u%data(jj)%ele_start_name, 'ELE_START ELEMENT', u, u%data(jj), u%data(jj)%ix_ele_start)

    ! eval_point
    call match_word (datum(ix1+jj-n1)%eval_point, anchor_pt_name(1:), ix, can_abbreviate = .false.)
    if (ix == 0) then
      call out_io (s_abort$, r_name, 'EVAL_POINT UNRECOGNIZED: ' // datum(ix1+jj-n1)%eval_point)
      stop
    endif
    u%data(jj)%eval_point = ix

    jj = jj + 1
  enddo

  deallocate (ind)

!-----------------------------------------
! use_same_lat_eles_as

elseif (use_same_lat_eles_as /= '') then
  call string_trim (use_same_lat_eles_as, name, ix)
  call tao_find_data (err, name, d1_array = d1_array, ix_uni = u%ix_uni)
  if (err .or. size(d1_array) /= 1) then
    call out_io (s_error$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
    return
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
  u%data(n1:n2)%meas_value      = 0

  u%data(n1:n2)%error_rms       = datum(ix1:ix2)%error_rms
  u%data(n1:n2)%invalid_value   = datum(ix1:ix2)%invalid_value
  u%data(n1:n2)%spin_map%axis_input  = datum(ix1:ix2)%spin_axis

  if (default_data_source /= '')  u%data(n1:n2)%data_source = default_data_source
  if (default_merit_type /= '')   u%data(n1:n2)%merit_type = default_merit_type
  if (default_weight /= 0)        u%data(n1:n2)%weight = default_weight

  do n = n1, n2
    ix = ix_min_data + n - n1

    if (default_data_type /= '')     u%data(n)%data_type   = trim(default_data_type)
    if (datum(ix)%data_type /= '')   u%data(n)%data_type   = trim(datum(ix)%data_type)
    if (u%data(n)%data_type == '')   u%data(n)%data_type   = trim(d2_data%name) // '.' // trim(d1_data%name)

    if (datum(ix)%weight /= 0)       u%data(n)%weight      = datum(ix)%weight
    if (datum(ix)%data_source /= '') u%data(n)%data_source = datum(ix)%data_source
    if (datum(ix)%merit_type /= '')  u%data(n)%merit_type  = datum(ix)%merit_type
  enddo

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

  ! If a global data type (no associated lattice element) then just construct a single datum.

  if (ix_max_data == int_garbage$) then
    has_associated_ele = tao_datum_has_associated_ele(tao_d2_d1_name(d1_this))
    if (has_associated_ele == yes$) then
      call out_io (s_error$, r_name, 'NO DATA FOUND FOR: ' // tao_d2_d1_name(d1_this))
      return
    else
      ix_min_data = 1
      ix_max_data = 1
    endif
  endif

  n1 = u%n_data_used + 1
  n2 = u%n_data_used + ix_max_data - ix_min_data + 1
  ix1 = ix_min_data
  ix2 = ix_max_data
  call tao_allocate_data_array (u, n2)

  ! Transfer info from the input structure

  u%data(n1:n2)%good_user        = datum(ix1:ix2)%good_user
  u%data(n1:n2)%good_opt         = datum(ix1:ix2)%good_opt
  u%data(n1:n2)%weight           = datum(ix1:ix2)%weight
  u%data(n1:n2)%ele_name         = datum(ix1:ix2)%ele_name
  u%data(n1:n2)%invalid_value    = datum(ix1:ix2)%invalid_value
  u%data(n1:n2)%error_rms        = datum(ix1:ix2)%error_rms
  u%data(n1:n2)%ele_ref_name     = datum(ix1:ix2)%ele_ref_name
  u%data(n1:n2)%ele_start_name   = datum(ix1:ix2)%ele_start_name
  u%data(n1:n2)%ix_bunch         = datum(ix1:ix2)%ix_bunch
  u%data(n1:n2)%data_source      = datum(ix1:ix2)%data_source
  u%data(n1:n2)%merit_type       = datum(ix1:ix2)%merit_type
  u%data(n1:n2)%spin_map%axis_input   = datum(ix1:ix2)%spin_axis
  u%data(n1:n2)%s_offset         = datum(ix1:ix2)%s_offset

  ! Find elements associated with the data

  do j = n1, n2
    ix = j - n1 + ix1
    ! use default_data_type if given, if not, auto-generate the data_type
    u%data(j)%data_type        = datum(ix)%data_type
    if (default_data_type == '') default_data_type = trim(d2_data%name) // '.' // trim(d1_data%name)
    if (u%data(j)%data_type == '') u%data(j)%data_type = trim(default_data_type)

    call match_word (datum(j+ix1-n1)%eval_point, anchor_pt_name(1:), ix, can_abbreviate = .false.)
    if (ix == 0) then
      call out_io (s_abort$, r_name, 'EVAL_POINT UNRECOGNIZED: ' // datum(j+ix1-n1)%eval_point)
      stop
    endif

    u%data(j)%eval_point = ix
    u%data(j)%exists = .true.

    call find_this_element (u%data(j)%ele_name,       'LATTICE ELEMENT',   u, u%data(j), u%data(j)%ix_ele, u%data(j)%ix_branch) 
    call find_this_element (u%data(j)%ele_ref_name,   'ELE_REF ELEMENT',   u, u%data(j), u%data(j)%ix_ele_ref)
    call find_this_element (u%data(j)%ele_start_name, 'ELE_START ELEMENT', u, u%data(j), u%data(j)%ix_ele_start)
  enddo
endif

u%data(n1:n2)%meas_value = datum(ix1:ix2)%meas
where (u%data(n1:n2)%meas_value == real_garbage$)  ! where %meas_value was not set
  u%data(n1:n2)%meas_value = 0  
elsewhere
  u%data(n1:n2)%good_meas = .true.
end where

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! Point the %data back to the d1_data_struct

call tao_point_d1_to_data (d1_this, u%data(n1:n2), ix_min_data)

! In a d1_data array, not all the datums need to exist. 
! If a datum is not associated with an element, that generally means that
! it does not exist. There are, however, a few exceptions. EG: unstable.lattice, etc.
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

  if (substr(dat%data_type, 1, 9) == 'unstable_') then
    call out_io (s_fatal$, r_name, &
         '"unstable_XXX" data source name (with an underscore) needs to be changed to "unstable.XXX" (with a dot).')
    stop
  endif

  ! Expressions may contain multiple norml.h terms

  u%calc%srdt_for_data = max(u%calc%srdt_for_data, tao_srdt_calc_needed(dat%data_type, dat%data_source))

  ix = 1
  do
    k = index(dat%data_type(ix:), 'normal.h.') 
    if (k == 0) exit
    
    if (dat%ix_branch /= 0) then
      call out_io (s_error$, r_name, 'EVALUATING A DATUM OF TYPE: ' // dat%data_type, 'ON A BRANCH NOT YET IMPLEMENTED!')        
    endif

    norm_form => u%model%tao_branch(dat%ix_branch)%bmad_normal_form
    n_ht => u%model%tao_branch(dat%ix_branch)%n_hterms

    h_str = substr(dat%data_type, ix+k+8, ix+k+13)
    found = .false.
    if (allocated(norm_form%h)) found = any(h_str == norm_form%h(1:n_ht)%id)

    if (.not. found) then
      n_ht = n_ht + 1
      call move_alloc (norm_form%h, h_temp)
      allocate (norm_form%h(n_ht))
      if (n_ht > 1) norm_form%h(1:n_ht-1) = h_temp
      norm_form%h(n_ht)%id = h_str
    endif

    if (len(dat%data_type) < ix + 20) exit
    ix = ix + k + 14
  enddo

  if (tao_rad_int_calc_needed(dat%data_type, dat%data_source)) then
    u%calc%rad_int_for_data = .true. 
    if (dat%ix_branch /= 0) then
      call out_io (s_error$, r_name, 'EVALUATING A DATUM OF TYPE: ' // dat%data_type, 'ON A BRANCH NOT YET IMPLEMENTED!')
    endif
  endif

  if (tao_lat_sigma_calc_needed(dat%data_type, dat%data_source)) u%calc%lat_sigma_for_data = .true. 
  if (tao_spin_matrices_calc_needed(dat%data_type, dat%data_source)) u%calc%spin_matrices = .true. 

  ! Some data types are global and are not associated with a particular element. Check for this.

  dat%exists = tao_data_sanity_check (dat, dat%exists, default_data_type)
  if (tao_chrom_calc_needed(dat%data_type, dat%data_source)) u%calc%chrom_for_data = .true.

enddo

if (.not. any(u%data(n1:n2)%exists) .and. n2 >= n1) then
  call out_io (s_warn$, r_name, &
            'Note: All datums in: ' // tao_d2_d1_name(d1_this), &
            'are marked as non-existent')
endif

end subroutine d1_data_stuffit

!------------------------------
! contains

subroutine find_this_element (ele_name, who, u, datum, ix_ele, ix_branch)

type (tao_universe_struct) u
type (tao_data_struct) datum
type (ele_pointer_struct), allocatable :: eles(:)
character(*) ele_name, who
integer ix_ele, n_loc
integer, optional :: ix_branch
logical exists

!

if (ele_name == '') return

call upcase_string (ele_name)
call lat_ele_locator (ele_name, u%design%lat, eles, n_loc)

if (n_loc == 0) then
  call out_io (s_warn$, r_name, who // ' NOT LOCATED: ' // ele_name, & 
                                 'FOR DATUM: ' // tao_datum_name(datum), &
                                 'WILL MARK THIS DATUM AS NOT EXISTING')
  datum%exists = .false.
  return
endif

if (n_loc > 1) then
  call out_io (s_error$, r_name, 'MULTIPLE LATTICE ELEMENTS OF THE SAME NAME: ' // ele_name, & 
                                 'FOR DATUM: ' // tao_datum_name(datum), &
                                 'WILL MARK THIS DATUM AS NOT EXISTING')
  datum%exists = .false.
  return
endif

ix_ele = eles(1)%ele%ix_ele
if (present(ix_branch)) ix_branch = eles(1)%ele%ix_branch

end subroutine find_this_element

end subroutine tao_init_data

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Defines what datums to evaluate at each element in specified universe

subroutine tao_init_data_end_stuff ()

implicit none

type (tao_universe_struct), pointer :: u
type (tao_data_struct), pointer :: data
type (tao_model_element_struct), pointer :: tao_model_ele(:)
integer i, ib, j, k, ix_ele, n_max

logical err

!

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  call tao_allocate_data_array (u, u%n_data_used, .true.) ! Trim u%data size
enddo

call tao_data_check (err)
if (err) stop

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

u%data%ix_uni = u%ix_uni

end subroutine tao_allocate_data_array

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine tao_d2_data_stuffit (u, d2_name, n_d1_data)

type (tao_universe_struct), target :: u
type (tao_d2_data_struct), pointer :: d2

integer i, nn, n_d1_data
character(*) d2_name
character(*), parameter :: r_name = 'tao_d2_data_stuffit'

! Setup another d2_data structure.

u%n_d2_data_used = u%n_d2_data_used + 1
nn = u%n_d2_data_used

if (size(u%d2_data) < nn) then
  call out_io (s_abort$, r_name, 'D2_DATA ARRAY OVERFLOW!')
  call err_exit
endif

d2 => u%d2_data(nn)

d2%name = d2_name
d2%ix_universe = u%ix_uni

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
subroutine tao_init_data_in_universe (u, n_d2_add, keep_existing_data)

implicit none

type (tao_universe_struct), target :: u
type (tao_d2_data_struct), allocatable :: d2_data(:)
type (tao_data_struct), allocatable :: data(:)
type (tao_d1_data_struct), pointer :: d1

integer n_d2_add, n_d2_data
integer i, j, k, n2
logical, optional :: keep_existing_data

!

if (logic_option(.false., keep_existing_data) .and. allocated(u%d2_data)) then
  n2 = size(u%d2_data)
  n_d2_data = n2 + n_d2_add
  call move_alloc(u%d2_data, d2_data)
  allocate (u%d2_data(n_d2_data))
  u%d2_data(1:n2) = d2_data
  do i = 1, n2
    do j = lbound(u%d2_data(i)%d1, 1), ubound(u%d2_data(i)%d1, 1)
      d1 => u%d2_data(i)%d1(j)
      d1%d2 => u%d2_data(i)
      do k = lbound(d1%d, 1), ubound(d1%d, 1)
        d1%d(k)%d1 => d1
      enddo
    enddo
  enddo

else
  n2 = 0
  n_d2_data = n2 + n_d2_add
  if (allocated(u%data)) deallocate(u%data)
  allocate (u%data(0))
  u%n_d2_data_used = 0      ! size of s%u(i)%d2_data(:) array
  u%n_data_used = 0         ! size of s%u(i)%data(:) array
  u%model%tao_branch%ix_rad_int_cache = 0
  u%design%tao_branch%ix_rad_int_cache = 0
  u%base%tao_branch%ix_rad_int_cache = 0

  if (n_d2_data == 0) return
  if (allocated(u%d2_data)) deallocate (u%d2_data)
  allocate (u%d2_data(n_d2_data))
endif

! 

u%d2_data(n2+1:n_d2_data)%name = ''  ! blank name means it doesn't exist

do j = 1, n_d2_data
  u%d2_data(j)%ix_d2_data = j
  if (j > n2) u%d2_data(j)%descrip = ''
enddo

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

end module
