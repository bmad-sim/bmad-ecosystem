!+
! Subroutine tao_init_global_and_universes (data_and_var_file)
!
! Subroutine to initialize the tao structures.
! If data_and_var_file is not in the current directory then it will be searched
! for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   data_and_var_file -- Character(*): Tao initialization file.

! Output:
!-

subroutine tao_init_global_and_universes (data_and_var_file)

  use tao_mod
  use tao_input_struct
  
  implicit none

  type (tao_d2_data_input) d2_data
  type (tao_d1_data_input) d1_data
  type (tao_data_input) data(n_data_minn:n_data_maxx) ! individual weight 
  type (tao_v1_var_input) v1_var
  type (tao_var_input) var(n_var_minn:n_var_maxx)
  type (tao_global_struct) global
  type (tao_d1_data_struct), pointer :: d1_ptr

  real(rp) :: default_weight        ! default merit function weight
  real(rp) :: default_step          ! default "small" step size
  real(rp) default_low_lim, default_high_lim

  integer ios, iu, i, j, k, ix, n_uni
  integer n_data_max, n_var_max, n_d2_data_max, n_v1_var_max
  integer n, n_universes, iostat, universe
  integer ix_min_var, ix_max_var, n_d1_data
  integer ix_min_data, ix_max_data, ix_d1_data
  integer, parameter :: ele_name$ = 1, ele_key$ = 2

  character(*) data_and_var_file
  character(40) :: r_name = 'tao_init_global_and_universes'
  character(200) file_name
  character(16) name,  default_universe, default_data_type
  character(16) default_merit_type, default_attribute
  character(100) line

  logical err
  logical counting, searching

  namelist / tao_params / global, n_data_max, n_var_max, n_d2_data_max, n_v1_var_max
         
  namelist / tao_d2_data / d2_data, n_d1_data, default_merit_type, universe
  namelist / tao_d1_data / d1_data, data, ix_d1_data, ix_min_data, &
                           ix_max_data, default_weight, default_data_type
  namelist / tao_var / v1_var, var, default_weight, default_step, &
                      ix_min_var, ix_max_var, default_universe, default_attribute, &
                      default_low_lim, default_high_lim

!-----------------------------------------------------------------------
! Init lattaces
! read global structure from tao_params namelist

  global%valid_plot_who(1:5) = &
                (/ 'model ', 'base  ', 'ref   ', 'design', 'meas  ' /)
  call tao_open_file ('TAO_INIT_DIR', data_and_var_file, iu, file_name)
  call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)
  read (iu, nml = tao_params)
  call out_io (s_blank$, r_name, 'Init: Read tao_params namelist')
  close (iu)

  s%global = global  ! transfer global to s
  
  n = size(s%u)
  do i = 1, size(s%u)
    call init_universe (s%u(i))
  enddo

  allocate (s%var(n_var_max))
  allocate (s%v1_var(n_v1_var_max))

  s%v1_var%name = ' '  ! blank name means it doesn't (yet) exist
  s%var(:)%good_opt  = .true.
  s%var(:)%exists    = .false.
  s%var(:)%good_var  = .true.
  s%var(:)%good_user = .false.

  s%n_var_used = 0
  s%n_v1_var_used = 0       ! size of s%v1_var(:) array

!-----------------------------------------------------------------------
! Init data

  call tao_open_file ('TAO_INIT_DIR', data_and_var_file, iu, file_name)

  do 
    d2_data%name = ' '      ! set default
    universe = 0
    default_merit_type = 'target'
    read (iu, nml = tao_d2_data, iostat = ios, err = 9100)
    if (ios < 0) exit         ! exit on end-of-file
    call out_io (s_blank$, r_name, &
                      'Init: Read tao_d2_data namelist: ' // d2_data%name)

    n_uni = universe      ! universe to use 
    if (n_uni == 0) then
      do i = 1, size(s%u)
        call d2_data_stuffit (s%u(i))
      enddo
    else
      call d2_data_stuffit (s%u(n_uni))
    endif

    do k = 1, n_d1_data
      default_weight = 0      ! set default
      default_data_type  = ' '
      data(:)%name       = ' '
      data(:)%data_type  = ' '
      data(:)%merit_type = ' '
      data(:)%ele_name   = ' '
      data(:)%ele2_name  = ' '
      data(:)%meas_value = 0
      data(:)%weight     = 0
      data(:)%good_data  = .false.
      read (iu, nml = tao_d1_data, err = 9150)
      if (ix_d1_data /= k) then
        write (line, '(a, 2i4)') 'IX_D1_DATA MISMATCH:', k, ix_d1_data
        call out_io (s_abort$, r_name, line, 'FOR: ' // d2_data%name)
        call err_exit
      endif
      call out_io (s_blank$, r_name, &
                      'Init: Read tao_d1_data namelist: ' // d1_data%name)
      if (n_uni == 0) then          ! 0 => use all universes
        do i = 1, size(s%u)
          call d1_data_stuffit (k, s%u(i), s%u(i)%n_d2_data_used)
        enddo
      else
        call d1_data_stuffit (k, s%u(n_uni), s%u(n_uni)%n_d2_data_used)
      endif
    enddo

  enddo

!-----------------------------------------------------------------------
! Init vars

  rewind (iu)

  do
    v1_var%name = " "         ! set default
    default_merit_type = 'limit'
    default_weight = 0     ! set default
    default_step = 0       ! set default
    default_attribute = ' '
    default_universe = ' '
    default_low_lim = -1e30
    default_high_lim = 1e30
    var%name = ' '
    var%ele_name = ' '
    var%merit_type = ' '
    var%weight = 0         ! set default
    var%step = 0           ! set default
    var%attribute = ' '
    var%universe = ' '
    var%low_lim = default_low_lim
    var%high_lim = default_high_lim

    read (iu, nml = tao_var, iostat = ios, err = 9200)
    if (ios < 0) exit         ! exit on end-of-file
    call out_io (s_blank$, r_name, &
                        'Init: Read tao_var namelist: ' // v1_var%name)

    if (v1_var%name == ' ') cycle

    if (default_universe == ' ' .and. all(var%universe == ' ')) &
                                                  default_universe = 'gang'

    if (default_universe == 'clone') then
      do i = 1, size(s%u)
        call var_stuffit_common
        write (s%v1_var(s%n_v1_var_used)%name, '(2a, i1)') &
                                s%v1_var(s%n_v1_var_used)%name, ';', i
        call var_stuffit (i)
      enddo

    elseif (default_universe == 'gang') then
      call var_stuffit_common
      call var_stuffit_all_uni

    else
      if (default_universe == ' ') then
        n = -1
      else
        read (default_universe, *, iostat = ios) n
        if (ios /= 0) then
          call out_io (s_abort$, r_name, &
              'CANNOT READ DEFAULT_UNIVERSE INDEX: ' // default_universe, &
              'FOR VARIABLE: ' // v1_var%name)
          call err_exit
        endif
      endif
      call var_stuffit_common
      call var_stuffit (n)
    endif

  enddo

  close (iu)
  return

! namelist read error.

9100 continue
  call out_io (s_error$, r_name, 'TAO_D2_DATA NAMELIST READ ERROR.')
  rewind (iu)
  do
    read (iu, nml = tao_d2_data)  ! force printing of error message
  enddo

9150 continue
  call out_io (s_error$, r_name, 'TAO_D1_DATA NAMELIST READ ERROR.')
  rewind (iu)
  do
    read (iu, nml = tao_d1_data)  ! force printing of error message
  enddo

9200 continue
  call out_io (s_error$, r_name, 'TAO_VAR NAMELIST READ ERROR.')
  rewind (iu)
  do
    read (iu, nml = tao_var)  ! force printing of error message
  enddo

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
contains

subroutine init_universe (u)

  type (tao_universe_struct) :: u

!

  n = u%design%n_ele_max
  allocate (u%model_orb(0:n), u%design_orb(0:n), u%base_orb(0:n))

  u%n_d2_data_used = 0      ! size of s%u(i)%d2_data(:) array
  u%n_data_used = 0         ! size of s%u(i)%data(:) array

! For linacs, specify initial conditions

  u%design_orb(0)%vec = 0.0
  
  call twiss_and_track (u%design, u%design_orb)
  u%model  = u%design; u%model_orb  = u%design_orb
  u%base = u%design; u%base_orb = u%design_orb

! allocate and set defaults

  if (n_d2_data_max /= 0) then
    allocate (u%d2_data(n_d2_data_max))
    u%d2_data%name = ' '  ! blank name means it doesn't exist
  endif

  if (n_data_max /= 0) then
    allocate (u%data(n_data_max))
    u%data(:)%exists = .false.       ! set default
    u%data(:)%good_data  = .false.   ! set default
    u%data(:)%good_ref   = .false.   ! set default
    u%data(:)%good_user  = .true.    ! set default
    u%data(:)%good_opt   = .true.
    u%data(:)%merit_type = 'target'  ! set default
    u%data(:)%ele_name   = ' '
    u%data(:)%ix_ele     = -1
    u%data(:)%ele2_name  = ' '
    u%data(:)%ix_ele2    = -1
  endif

! This is needed to keep the totalview debugger happy.

  allocate (u%dmodel_dvar(1,1))
  
end subroutine init_universe

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine d2_data_stuffit (u)

type (tao_universe_struct), target :: u

integer nn

! Setup another d2_data structure.

  u%n_d2_data_used = u%n_d2_data_used + 1
  nn = u%n_d2_data_used

  if (size(u%d2_data) < nn) then
    call out_io (s_error$, r_name, &
              'N_D2_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif

  u%d2_data(nn)%name = d2_data%name 

! allocate memory for the u%d1_data structures

  allocate(u%d2_data(nn)%d1(n_d1_data))

end subroutine

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine d1_data_stuffit (i_d1, u, n_d2)

type (tao_universe_struct), target :: u
integer i, n1, n2, ix, k, ix1, ix2, j, jj, n_d2

integer i_d1, num_hashes

character(20) count_name1, count_name2, ix_char
character(20) search_string
character(20) fmt

logical, allocatable :: found_one(:)

!

u%d2_data(n_d2)%d1(i_d1)%d2 => u%d2_data(n_d2)  ! point back to the parent

! are we counting elements and forming data names?

if (index(data(0)%name, 'COUNT:') /= 0) then
  counting = .true.
  call form_count_name (data(0)%name(7:), num_hashes, count_name1, count_name2)
! if using SAME: then use the specified d1_data to count datums below...
elseif (index(data(0)%ele_name, 'SAME:') .eq. 0) then
  counting = .false.
  n1 = u%n_data_used + 1
  n2 = u%n_data_used + ix_max_data - ix_min_data + 1
  ix1 = ix_min_data
  ix2 = ix_max_data
  u%n_data_used = n2
  if (n2 > size(u%data)) then
    call out_io (s_abort$, r_name, &
                'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif
endif

u%d2_data(n_d2)%d1(i_d1)%name = d1_data%name  ! stuff in the data

! now check if we are searching for elements or repeating elements
! and record the element names in the data structs
    
if (index(data(0)%ele_name, 'SEARCH') .ne. 0) then
  allocate (found_one(u%design%n_ele_max))
  if (index(data(0)%ele_name, 'SEARCH_KEY:') .ne. 0) then
    call string_trim(data(0)%ele_name(12:), search_string, ix)
    call find_elements (u, search_string, ele_key$, found_one)
  elseif  (index(data(0)%ele_name, 'SEARCH:') .ne. 0) then
    call string_trim(data(0)%ele_name(8:), search_string, ix)
    call find_elements (u, search_string, ele_name$, found_one)
  else
    call out_io (s_abort$, r_name, 'Syntax Error in data(0)%ele_name SEARCH string')
    call err_exit
  endif
  ! finish finding data array limits
  if (counting) then
    n1 = u%n_data_used + 1
    n2 = u%n_data_used + count(found_one)
    ix1 = ix_min_data
    ix2 = (count(found_one) - (1-ix_min_data))
    u%n_data_used = n2
    if (n2 > size(u%data)) then
      call out_io (s_abort$, r_name, &
                  'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
      call err_exit
    endif
  endif
  ! get element names
  jj = n1
  do j = 1, size(found_one)
    if (found_one(j)) then
      if (jj .gt. n2) then
        call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
        call err_exit
      endif
      u%data(jj)%ele_name = u%design%ele_(j)%name
      u%data(jj)%ix_ele   = j
      u%data(jj)%exists   = .true.
      jj = jj + 1
    endif
  enddo
  u%data(n1:n2)%meas_value = 0 
  u%data(n1:n2)%data_type  = default_data_type
  u%data(n1:n2)%merit_type = default_merit_type 
  u%data(n1:n2)%good_data  = .false.

elseif (index(data(0)%ele_name, 'SAME:') /= 0) then
  call string_trim (data(0)%ele_name(6:), name, ix)
  call tao_find_data (err, u, name, d1_ptr = d1_ptr)
  if (err) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
    call err_exit
  endif
  n1 = u%n_data_used + 1
  n2 = n1 + size(d1_ptr%d) - 1
  u%n_data_used = n2
  if (n2 > size(u%data)) then
    call out_io (s_abort$, r_name, &
                'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif
  u%data(n1:n2)%data_type  = default_data_type
  u%data(n1:n2)%ele_name   = d1_ptr%d%ele_name
  u%data(n1:n2)%ix_ele     = d1_ptr%d%ix_ele
  u%data(n1:n2)%ele2_name  = d1_ptr%d%ele2_name
  u%data(n1:n2)%ix_ele2    = d1_ptr%d%ix_ele2
  u%data(n1:n2)%exists     = d1_ptr%d%exists
 
  u%data(n1:n2)%meas_value = d1_ptr%d%meas_value
  u%data(n1:n2)%merit_type = d1_ptr%d%merit_type
  u%data(n1:n2)%good_data  = d1_ptr%d%good_data
  u%data(n1:n2)%weight     = d1_ptr%d%weight
else
  u%data(n1:n2)%ele_name  = data(ix1:ix2)%ele_name
  u%data(n1:n2)%ele2_name = data(ix1:ix2)%ele2_name
  do j = n1, n2
    if (u%data(j)%ele_name == ' ') cycle
    call str_upcase (u%data(j)%ele_name, u%data(j)%ele_name)
    call element_locator (u%data(j)%ele_name, u%design, ix)
    if (ix < 0) then
      call out_io (s_abort$, r_name, 'ELEMENT NOT LOCATED: ' // &
                                                       u%data(j)%ele_name)
      call err_exit
    endif
    u%data(j)%ix_ele = ix
    u%data(j)%exists = .true.

    if (u%data(j)%ele2_name == ' ') cycle
    call str_upcase (u%data(j)%ele2_name, u%data(j)%ele2_name)
    call element_locator (u%data(j)%ele2_name, u%design, ix)
    if (ix < 0) then
      call out_io (s_abort$, r_name, 'ELEMENT2 NOT LOCATED: ' // &
                                                       u%data(j)%ele2_name)
      call err_exit
    endif
    u%data(j)%ix_ele2 = ix
  enddo
  u%data(n1:n2)%meas_value = data(ix1:ix2)%meas_value
  u%data(n1:n2)%data_type  = data(ix1:ix2)%data_type
  u%data(n1:n2)%merit_type = data(ix1:ix2)%merit_type
  u%data(n1:n2)%good_data  = data(ix1:ix2)%good_data
  u%data(n1:n2)%weight     = data(ix1:ix2)%weight
endif

! use default_data_type if given, if not, auto-generate the data_type
if (default_data_type .eq. ' ') then
  where (u%data(n1:n2)%data_type == ' ') u%data(n1:n2)%data_type = &
                            trim(d2_data%name) // ':' // d1_data%name
else
  where (u%data(n1:n2)%data_type == ' ') u%data(n1:n2)%data_type = &
                                                    default_data_type
endif


			    

! Create data names

if (index(data(0)%name, 'COUNT:') /= 0) then
  jj = ix1
  do j = n1, n2
    if (jj .gt. ix2) then
      call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
      call err_exit
    endif
    write(fmt, '(a,i1.1,a,i1.1,a)') '(a, I', num_hashes, '.', num_hashes, ', a)'
    write(u%data(j)%name, fmt) trim(count_name1), jj, trim(count_name2)
    jj = jj + 1
  enddo

elseif (index(data(0)%name, 'SAME:') /= 0) then
  call string_trim (data(0)%name(6:), name, ix)
  call tao_find_data (err, u, name, d1_ptr = d1_ptr)
  if (err) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
    call err_exit
  endif
  n2 = n1 + size(d1_ptr%d) - 1
  u%data(n1:n2)%name = d1_ptr%d%name
else
  u%data(n1:n2)%name = data(ix1:ix2)%name
endif


! now for some family guidance...
! point the children to the grandchildren in the big data array

call tao_point_d1_to_data (u%d2_data(n_d2)%d1(i_d1)%d, &
                                      u%data(n1:n2), ix_min_data, n1)

! point the %data back to the d1_data_struct

do j = n1, n2
  u%data(j)%d1 => u%d2_data(n_d2)%d1(i_d1)
  if (u%data(j)%weight == 0) u%data(j)%weight = default_weight
  if (u%data(j)%merit_type == ' ') u%data(j)%merit_type =  default_merit_type
enddo

! point the children back to the mother    

u%d2_data(n_d2)%d1(i_d1)%d2 => u%d2_data(n_d2)
if (allocated(found_one)) deallocate (found_one)  
  
end subroutine d1_data_stuffit

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine var_stuffit (ix_u_in)

  type (tao_var_struct), pointer :: s_var
  integer i, ix_u, n, ios, ix_u_in

! point the children back to the mother

  n = s%n_v1_var_used
  do i = lbound(s%v1_var(n)%v, 1), ubound(s%v1_var(n)%v, 1)
    s_var => s%v1_var(n)%v(i)

    ! universe to use
    ix_u = ix_u_in
    if (.not. (counting .and. searching)) then
      if (var(i)%universe /= ' ') then
        read (var(i)%universe, *, iostat = ios) ix_u
        if (ios /= 0) then
          call out_io (s_abort$, r_name, &
              'CANNOT READ DEFAULT_UNIVERSE INDEX: ' // default_universe, &
              'FOR VARIABLE: ' // v1_var%name)
          call err_exit
        endif
      endif
    endif
     
    allocate (s_var%this(1))
    if (s_var%ele_name == ' ') cycle
    call tao_pointer_to_var_in_lattice (s_var, s_var%this(1), ix_u)
    s_var%model_value = s_var%this(1)%model_ptr
    s_var%design_value = s_var%model_value
    s_var%base_value = s_var%this(1)%base_ptr
    s_var%exists = .true.
  enddo

end subroutine var_stuffit


!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine var_stuffit_all_uni 

  type (tao_universe_struct), pointer :: u
  type (tao_var_struct), pointer :: s_var

  integer i, j, n1, n2, ix1, ix2
  logical err

! point the children back to the mother


  n = s%n_v1_var_used
  
  do i = lbound(s%v1_var(n)%v, 1), ubound(s%v1_var(n)%v, 1)
    s_var => s%v1_var(n)%v(i)
    allocate (s_var%this(size(s%u)))
    if (s_var%ele_name == ' ') cycle
    do j = 1, size(s%u)
      u => s%u(j)
      call tao_pointer_to_var_in_lattice (s_var, s_var%this(j), j)
    enddo
    s_var%model_value = s_var%this(1)%model_ptr
    s_var%design_value = s_var%this(1)%model_ptr
    s_var%base_value = s_var%this(1)%base_ptr
    s_var%exists = .true.
  enddo

end subroutine

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! stuff common to all universes

subroutine var_stuffit_common

character(20) count_name1, count_name2, ix_char
character(20) fmt, search_string

integer i, j, jj, nn, n1, n2, ix1, ix2, num_hashes, ix

logical, allocatable :: found_one(:)

! count number of v1 entries

  s%n_v1_var_used = s%n_v1_var_used + 1
  nn = s%n_v1_var_used

  ! are we searching for and couting elements?
  if (index(var(0)%name, 'COUNT:') /= 0) then
    counting = .true.
    call form_count_name (var(0)%name(7:), num_hashes, count_name1, count_name2)
    if (index(var(0)%ele_name, 'SEARCH') /= 0) then
      searching = .true.
      allocate (found_one(s%u(1)%design%n_ele_max))
      if (index(var(0)%ele_name, 'SEARCH_KEY:') .ne. 0) then
        call string_trim(var(0)%ele_name(12:), search_string, ix)
        call find_elements (s%u(1), search_string, ele_key$, found_one)
      elseif (index(var(0)%ele_name, 'SEARCH:') .ne. 0) then 
        call string_trim(var(0)%ele_name(8:), search_string, ix)
        call find_elements (s%u(1), search_string, ele_name$, found_one)
      else
	call out_io (s_abort$, r_name, 'Syntax Error in var(0)%ele_name SEARCH string')
	call err_exit
      endif
    else
      call out_io (s_abort$, r_name, 'If you are counting elements you should &
                                      also be searching for them')
      call err_exit
    endif
    n1 = s%n_var_used + 1
    n2 = s%n_var_used + count(found_one)
    ix1 = ix_min_var
    ix2 = (count(found_one) - (1-ix_min_var))
    s%n_var_used = n2
    !get element names
    jj = n1
    do j = 1, size(found_one)
      if (found_one(j)) then
        if (jj .gt. n2) then
          call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
          call err_exit
        endif
        s%var(jj)%ele_name = s%u(1)%design%ele_(j)%name
        jj = jj + 1
      endif
    enddo
    ! Create var names
    jj = ix1
    do j = n1, n2
      if (jj .gt. ix2) then
        call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
        call err_exit
      endif
      write(fmt, '(a,i1.1,a,i1.1,a)') '(a, I', num_hashes, '.', num_hashes, ', a)'
      write(s%var(j)%name, fmt) trim(count_name1), jj, trim(count_name2)
      jj = jj + 1
    enddo
    s%var(n1:n2)%attrib_name = default_attribute
    s%var(n1:n2)%weight = default_weight
    s%var(n1:n2)%step = default_step
    s%var(n1:n2)%merit_type = default_merit_type
    s%var(n1:n2)%low_lim = default_low_lim
    s%var(n1:n2)%high_lim = default_high_lim
  else
    counting = .false.
    n1 = s%n_var_used + 1
    n2 = s%n_var_used + ix_max_var - ix_min_var + 1
    ix1 = ix_min_var
    ix2 = ix_max_var
 
    s%n_var_used = n2
 
    s%var(n1:n2)%ele_name = var(ix1:ix2)%ele_name
    s%var(n1:n2)%name = var(ix1:ix2)%name

    s%var(n1:n2)%attrib_name = var(ix1:ix2)%attribute
    where (s%var(n1:n2)%attrib_name == ' ') s%var(n1:n2)%attrib_name = default_attribute
 
    s%var(n1:n2)%weight = var(ix1:ix2)%weight
    where (s%var(n1:n2)%weight == 0) s%var(n1:n2)%weight = default_weight
 
    s%var(n1:n2)%step = var(ix1:ix2)%step
    where (s%var(n1:n2)%step == 0) s%var(n1:n2)%step = default_step
 
    s%var(n1:n2)%merit_type = var(ix1:ix2)%merit_type
    where (s%var(n1:n2)%merit_type == ' ') s%var(n1:n2)%merit_type = default_merit_type
 
    s%var(n1:n2)%low_lim = var(ix1:ix2)%low_lim
    where (s%var(n1:n2)%low_lim == -1e30) s%var(n1:n2)%low_lim = default_low_lim
 
    s%var(n1:n2)%high_lim = var(ix1:ix2)%high_lim
    where (s%var(n1:n2)%high_lim == 1e30) s%var(n1:n2)%high_lim = default_high_lim
  endif
 
  s%v1_var(nn)%name = v1_var%name

! now for some family guidance...
! point the v1_var mother to the appropriate children in the big data array

  call tao_point_v1_to_var (s%v1_var(nn), s%var(n1:n2), ix_min_var, n1)

  if (abs(lbound(s%v1_var(s%n_v1_var_used)%v, 1) - &
                          ubound(s%v1_var(s%n_v1_var_used)%v, 1)) .gt. 1000) then
    call out_io (s_blank$, r_name, "Initilizing a large number of variables.")
    call out_io (s_blank$, r_name, "This may take a while...")
    call out_io (s_blank$, r_name, " ")
  endif

end subroutine

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! this forms the name used in the variable or data where the number of hashes is
! replaced by the element index

subroutine form_count_name (count_name, num_hashes, count_name1, count_name2)

implicit none

character(*) count_name, count_name1, count_name2
integer num_hashes, ix


  ! 'COUNT:' is 6 characters long
  call string_trim (count_name, count_name1, ix)
  ix = index (count_name1, '#')
  if (ix .eq. 0) then
    call out_io (s_abort$, r_name, "WHEN USING 'COUNT:' MUST HAVE '#' &
                    WILDCARD IN NAME")
    call err_exit
  endif
  call tao_count_strings (count_name1, '#', num_hashes)
  count_name2 = count_name1(ix+num_hashes:)
  count_name1 = count_name1(:ix-1) 

end subroutine form_count_name

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! This searches the lattice for the specified element and flags found_one(:)
!
! Attribute can be either ele_name$ or ele_key$

subroutine find_elements (u, search_string, attribute, found_one)

type (tao_universe_struct) :: u
character(*) search_string
integer attribute, key, found_key
logical found_one(:)

integer j

  found_one = .false.
  if (attribute .eq. ele_name$) then
    do j = 1, u%design%n_ele_max
      if (match_wild(u%design%ele_(j)%name, search_string)) &
      found_one(j) = .true.
    enddo
  elseif (attribute .eq. ele_key$) then
    found_key = 0
    call upcase_string(search_string)
    do j = 1, size(key_name)
      if (key_name(j)(1:len(trim(search_string))) .eq. search_string) then
!      if (index(key_name(j), trim(search_string)) .ne. 0) then
	found_key = found_key + 1
       	key = j
      endif
    enddo
    if (found_key .ne. 1) then
      call out_io (s_abort$, r_name, "Ambiguous or non-existant key name")
      call err_exit
    endif
    do j = 1, u%design%n_ele_max
      if (u%design%ele_(j)%key .eq. key) &
      found_one(j) = .true.
    enddo
  else 
    !bug in call to subroutine
    call out_io (s_abort$, r_name, "Internal Error in find_elements!")
    call err_exit
  endif

end subroutine find_elements

end subroutine
