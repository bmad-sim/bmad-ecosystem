!+
! Subroutine tao_init_global_and_universes (s, data_and_var_file)
!
! Subroutine to initialize the tao structures.
! If data_and_var_file is not in the current directory then it will be searched
! for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   data_and_var_file -- Character(*): Tao initialization file.

! Output:
!   s -- Tao_super_universe_struct:
!-

subroutine tao_init_global_and_universes (s, data_and_var_file)

  use tao_mod
  use tao_input_struct
  
  implicit none

  integer ios, iu, i, j

  type (tao_super_universe_struct) s
  type (tao_d2_data_input) d2_data
  type (tao_d1_data_input) d1_data(n_d1_data_maxx)
  type (tao_data_input) data(n_d1_data_maxx)
  type (tao_v1_var_input) v1_var
  type (tao_var_input) var
  type (tao_global_struct) global
  type (tao_d1_data_struct), pointer :: d1_ptr

  type n_uni_struct
    integer d2_data
    integer data
  end type
  type (n_uni_struct), pointer :: nu(:)

  integer n_data_max, n_var_max, n_d2_data_max, n_v1_var_max
  integer n, n_universes, iostat

  character(*) data_and_var_file
  character(40) :: r_name = 'TAO_INIT_GLOBAL_AND_UNIVERSES'
  character(200) file_name
  character(16) name

  logical err

  namelist / tao_params / global, n_data_max, n_var_max, n_d2_data_max, n_v1_var_max
         
  namelist / tao_data / d2_data, d1_data, data
  namelist / tao_var / v1_var, var

! Init lattaces
! read global structure from tao_params namelist
! nu(i) keeps track of the sizes of allocatable pointers in universe s%u(i).

  global%valid_plot_who(1:5) =  (/ 'model ', 'base  ', 'ref   ', 'design', 'data  ' /)
  call tao_open_file ('TAO_INIT_DIR', data_and_var_file, iu, file_name)
  read (iu, nml = tao_params)
  call out_io (s_blank$, r_name, 'Init: Read tao_params namelist')
  close (iu)

  s%global = global  ! transfer global to s
  
  n = size(s%u)
  allocate (nu(n))
  nu%d2_data = 0      ! size of s%u(i)%d2_data(:) array
  nu%data = 0         ! size of s%u(i)%data(:) array

  do i = 1, size(s%u)
    call init_universe (s%u(i))
  enddo

  allocate (s%var(n_var_max))
  allocate (s%v1_var(n_v1_var_max))

  s%v1_var(:)%good_opt = .true.
  s%var(:)%exists    = .false.
  s%var(:)%good_var  = .true.
  s%var(:)%good_user = .false.
  s%var(:)%merit_type = 'target'

  s%n_var_used = 0
  s%n_v1_var_used = 0       ! size of s%v1_var(:) array


! Init data

  call tao_open_file ('TAO_INIT_DIR', data_and_var_file, iu, file_name)

  do
    d1_data(:)%name = ' '      ! set default
    data(:)%default_weight = 0.0      ! set default
    do i = lbound(d1_data, 1), ubound(d1_data, 1)
      data(i)%weight = 0.0        ! set default
    enddo
    read (iu, nml = tao_data, iostat = ios, err = 9100)
    if (ios < 0) exit         ! exit on end-of-file
    call out_io (s_blank$, r_name, &
                      'Init: Read tao_data namelist: ' // d2_data%name)
    n = d2_data%universe      ! universe to use 
    if (n == 0) then          ! 0 => use all universes
      do i = 1, size(s%u)
        call data_stuffit (s%u(i), nu(i), d2_data, d1_data, data)
      enddo
    else
      call data_stuffit (s%u(n), nu(n), d2_data, d1_data, data)
    endif
  enddo

! Init vars

  rewind (iu)

  do
    v1_var%name = " "         ! set default
    var%default_weight = 0     ! set default
    var%default_step = 0       ! set default
    var%weight = 0         ! set default
    var%step = 0           ! set default
    read (iu, nml = tao_var, iostat = ios, err = 9200)
    if (ios < 0) exit         ! exit on end-of-file
    call out_io (s_blank$, r_name, &
                        'Init: Read tao_var namelist: ' // v1_var%name)

    if (v1_var%universe == 'all') then
      do i = 1, size(s%u)
        call var_stuffit (i, nu(i), v1_var, var)
      enddo
    elseif (v1_var%universe == 'common') then
      call var_stuffit_common (v1_var, var)
    else
      read (v1_var%universe, *, iostat = ios) n
      if (ios /= 0) then
        call out_io (s_fatal$, r_name, &
              'CANNOT READ UNIVERSE INDEX: ' // v1_var%universe, &
              'FOR VARIABLE: ' // v1_var%name)
        call err_exit
      endif
      call var_stuffit (n, nu(n), v1_var, var)
    endif

  enddo

  close (iu)
  return

! namelist read error.

9100 continue
  call out_io (s_error$, r_name, 'TAO_DATA NAMELIST READ ERROR.')
  rewind (iu)
  do
    read (iu, nml = tao_data)  ! force printing of error message
  enddo

9200 continue
  call out_io (s_error$, r_name, 'TAO_VAR NAMELIST READ ERROR.')
  rewind (iu)
  do
    read (iu, nml = tao_var)  ! force printing of error message
  enddo

!--------------------------------------------------------------------------
contains

subroutine init_universe (u)

  type (tao_universe_struct) :: u

!

  n = u%design%n_ele_max
  allocate (u%model_orb(0:n), u%design_orb(0:n), u%base_orb(0:n))

!For linacs, specify initial conditions

  u%design_orb(0)%vec = 0.0
  
  call twiss_and_track (u%design, u%design_orb)
  u%model  = u%design; u%model_orb  = u%design_orb
  u%base = u%design; u%base_orb = u%design_orb

! allocate and set defaults

  if (n_d2_data_max /= 0) then
    allocate (u%d2_data(n_d2_data_max))
    u%d2_data(:)%good_opt = .true.
  endif

  if (n_data_max /= 0) then
    allocate (u%data(n_data_max))
    u%data(:)%exists = .false.       ! set default
    u%data(:)%good_data  = .false.   ! set default
    u%data(:)%good_ref   = .false.   ! set default
    u%data(:)%good_user  = .true.    ! set default
    u%data(:)%merit_type = 'target'  ! set default
  endif
  
end subroutine init_universe

!----------------------------------------------------------------
! contains

subroutine data_stuffit (u, nu, d2_data, d1_data, data)

type (tao_universe_struct), target :: u
type (n_uni_struct) nu
integer i, nn, n1, n2, ix, k, ix1, ix2, j, jj

type (tao_d2_data_input) d2_data
type (tao_d1_data_input) d1_data(:)
type (tao_data_input) data(:)

integer i_d2, i_d1, nd1

character(20) count_name1, count_name2, ix_char
character(20) search_string

logical counting, searching
logical, allocatable :: found_one(:)

! Count number of d2 data entries and
! transfer info from input structures to the universe structure.

  nu%d2_data = nu%d2_data + 1
  nn = nu%d2_data

  if (size(u%d2_data) < nn) then
    call out_io (s_error$, r_name, &
              'N_D2_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif

  u%d2_data(nn)%name = d2_data%name 

! allocate memory for the u%d1_data structures

  nd1 = count(d1_data(:)%name /= ' ') 
  allocate(u%d2_data(nn)%d1(nd1))

  do i = 1, nd1

    if (d1_data(i)%name == ' ') exit

    u%d2_data(nn)%d1(i)%d2 => u%d2_data(nn)  ! point back to the parent

! are we counting elements and forming data names?
    if (index(data(i)%name(0), 'COUNT:') /= 0) then
      counting = .true.
      count_name1 = data(i)%name(0)
      call string_trim (count_name1(7:), count_name1, ix)
      call string_trim (count_name1, count_name1, ix)
      ix = index (count_name1, '#')
      if (ix .eq. 0) then
      call out_io (s_abort$, r_name, "WHEN USING 'SAME:' MUST HAVE '#' &
                  WILDCARD IN NAME")
      call err_exit
      endif
      count_name2 = count_name1(ix+1:)
      count_name1 = count_name1(:ix-1) 
    else
      counting = .false.
      n1 = nu%data + 1
      n2 = nu%data + data(i)%ix_max - data(i)%ix_min + 1
      ix1 = data(i)%ix_min
      ix2 = data(i)%ix_max
      nu%data = n2
      if (n2 > size(u%data)) then
        call out_io (s_abort$, r_name, &
                'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
        call err_exit
      endif
    endif

    u%d2_data(nn)%d1(i)%name = d1_data(i)%name  ! stuff in the data

! now check if we are searching for elements or repeating elements
! and record the element names in the data structs
    
    if (index(data(i)%ele_name(0), 'SEARCH:') /= 0) then
      search_string = data(i)%ele_name(0)
      call string_trim (search_string(8:), search_string, ix)
      allocate (found_one(u%design%n_ele_max))
      found_one = .false.
      do j = 1, u%design%n_ele_max
        if (match_wild(u%design%ele_(j)%name, search_string)) &
        found_one(j) = .true.
        !keep track of element index in lattice
        u%design%ele_(j)%ix_pointer = j
      enddo
      !finish finding data array limits
      if (counting) then
        n1 = nu%data + 1
        n2 = nu%data + count(found_one)
        ix1 = data(i)%ix_min
        ix2 = (count(found_one) - (1-data(i)%ix_min))
        nu%data = n2
        if (n2 > size(u%data)) then
          call out_io (s_abort$, r_name, &
                  'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
          call err_exit
        endif
      endif
      !get element names
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
    elseif (index(data(i)%ele_name(0), 'SAME:') /= 0) then
      call string_trim (data(i)%ele_name(0)(6:), name, ix)
      call tao_find_data (err, u, name, d1_ptr = d1_ptr)
      if (err) then
        call out_io (s_abort$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
        call err_exit
      endif
      u%data(n1:n2)%ele_name = d1_ptr%d%ele_name
      u%data(n1:n2)%ix_ele   = d1_ptr%d%ix_ele
      u%data(n1:n2)%exists   = d1_ptr%d%exists
    else
      u%data(n1:n2)%ele_name = data(i)%ele_name(ix1:ix2)
      do j = n1, n2
        if (u%data(j)%ele_name == ' ') cycle
        call element_locator (u%data(j)%ele_name, u%design, ix)
        if (ix < 0) then
          call out_io (s_abort$, r_name, 'ELEMENT NOT LOCATED: ' // &
                                                       u%data(j)%ele_name)
          call err_exit
        endif
        u%data(j)%ix_ele = ix
        u%data(j)%exists = .true.
      enddo
    endif

    u%data(n1:n2)%class = trim(d2_data%name) // ':' // d1_data(i)%name

! Create data names
    if (index(data(i)%name(0), 'COUNT:') /= 0) then
      jj = ix1
      do j = n1, n2
        if (jj .gt. ix2) then
          call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
          call err_exit
        endif
        write (ix_char, '(I)') jj
        u%data(j)%name = trim(count_name1) // trim(ix_char) // count_name2
      jj = jj + 1
      enddo

    elseif (index(data(i)%name(0), 'SAME:') /= 0) then
      call string_trim (data(i)%name(0)(6:), name, ix)
      call tao_find_data (err, u, name, d1_ptr = d1_ptr)
      if (err) then
        call out_io (s_abort$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
        call err_exit
      endif
      u%data(n1:n2)%name = d1_ptr%d%name
    else
      u%data(n1:n2)%name = data(i)%name(ix1:ix2)      
    endif

    if (ix2-ix1 .le. size(data(i)%weight) .and. .not. counting) &
          u%data(n1:n2)%weight   = data(i)%weight(ix1:ix2)

! now for some family guidance...
! point the children to the grandchildren in the big data array

    call tao_point_d1_to_data (u%d2_data(nn)%d1(i)%d, &
                                      u%data(n1:n2), data(i)%ix_min, n1)

! point the %data back to the d1_data_struct

    do j = n1, n2
      u%data(j)%d1 => u%d2_data(nn)%d1(i)  
      if (u%data(j)%weight == 0) u%data(j)%weight = data(i)%default_weight
    enddo

! point the children back to the mother    

    u%d2_data(nn)%d1(i)%d2 => u%d2_data(nn)
    if (allocated(found_one)) deallocate (found_one)

  enddo
  
  
end subroutine data_stuffit

!----------------------------------------------------------------
! contains

subroutine var_stuffit (ix_u, nu, v1_var, var)

  type (n_uni_struct) nu
  type (tao_v1_var_input) v1_var
  type (tao_var_input) var

  integer i, j, nn, n1, n2, ix1, ix2, ix_u
  logical err

! count number of v1 entries

  s%n_v1_var_used = s%n_v1_var_used + 1
  nn = s%n_v1_var_used
  
! stuff in the data

  s%v1_var(nn)%name = v1_var%name

  if (v1_var%name == ' ') return
  n1 = s%n_var_used + 1
  n2 = s%n_var_used + var%ix_max - var%ix_min + 1
  ix1 = var%ix_min
  ix2 = var%ix_max

  s%n_var_used = n2

  s%var(n1:n2)%ele_name = var%ele_name(ix1:ix2)
  s%var(n1:n2)%name = var%name(ix1:ix2)

! now for some family guidance...
! point the v1_var mother to the appropriate children in the big data array

  call tao_point_v1_to_var (s%v1_var(nn)%v, s%var(n1:n2), var%ix_min, n1)

  s%var(n1:n2)%attrib_name = v1_var%attribute
  s%var(n1:n2)%weight = var%weight
  where (s%var(n1:n2)%weight == 0) s%var(n1:n2)%weight = var%default_weight
  s%var(n1:n2)%step = var%step
  where (s%var(n1:n2)%step == 0) s%var(n1:n2)%step = var%default_step

! point the children back to the mother

  do i = n1, n2
    s%var(i)%v1 => s%v1_var(nn)
    allocate (s%var(i)%this(1))
    if (s%var(i)%ele_name == ' ') cycle
    call tao_pointer_to_var_in_lattice (s, s%var(i), s%var(i)%this(1), ix_u)
    s%var(i)%model_value = s%var(i)%this(1)%model_ptr
    s%var(i)%design_value = s%var(i)%model_value
    s%var(i)%base_value = s%var(i)%this(1)%base_ptr
    s%var(i)%exists = .true.
  enddo

end subroutine


!----------------------------------------------------------------
! contains

subroutine var_stuffit_common (v1_var, var)

  type (tao_universe_struct), pointer :: u
  type (tao_v1_var_input) v1_var
  type (tao_var_input) var

  integer i, j, nn, n1, n2, ix1, ix2
  logical err

! stuff in the data

  if (v1_var%name == ' ') return

  s%n_v1_var_used = s%n_v1_var_used + 1
  nn = s%n_v1_var_used

  n1 = s%n_var_used + 1
  n2 = s%n_var_used + var%ix_max - var%ix_min + 1
  ix1 = var%ix_min
  ix2 = var%ix_max

  s%n_var_used = n2

  s%var(n1:n2)%ele_name = var%ele_name(ix1:ix2)
  s%var(n1:n2)%name = var%name(ix1:ix2)

! now for some family guidance...
! point the v1_var mother to the appropriate children in the big data array

  call tao_point_v1_to_var (s%v1_var(nn)%v, s%var(n1:n2), var%ix_min, n1)

  s%var(n1:n2)%attrib_name = v1_var%attribute
  s%var(n1:n2)%weight = var%weight
  where (s%var(n1:n2)%weight == 0) s%var(n1:n2)%weight = var%default_weight
  s%var(n1:n2)%step = var%step
  where (s%var(n1:n2)%step == 0) s%var(n1:n2)%step = var%default_step

! point the children back to the mother

  do i = n1, n2

    allocate (s%var(i)%this(size(s%u)))
    if (s%var(i)%ele_name == ' ') cycle
    s%var(i)%v1 => s%v1_var(nn)

    do j = 1, size(s%u)
      u => s%u(j)
      call tao_pointer_to_var_in_lattice (s, s%var(i), s%var(i)%this(j), j)
      s%var(i)%exists = .true.
    enddo

    s%var(i)%model_value = s%var(i)%this(1)%model_ptr
    s%var(i)%design_value = s%var(i)%this(1)%model_ptr
    s%var(i)%base_value = s%var(i)%this(1)%base_ptr

  enddo

end subroutine

end subroutine
