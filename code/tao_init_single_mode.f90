!+
! Subroutine tao_init_single_mode (single_mode_file)
!
! Subroutine to initialize the tao single mode stuff.
! If the single_mode_file is not in the current directory then it will be searched
! for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   single_mode_file -- Character(*): Tao initialization file.

! Output:
!-

subroutine tao_init_single_mode (single_mode_file)

  use tao_mod
  use tao_input_struct
  
  implicit none

  type (tao_key_input) key(500)

  integer ios, iu, i, j, n, n1, n2, nn, i_max, ix_u, ix

  character(*) single_mode_file
  character(40) :: r_name = 'tao_init_single_mode'
  character(200) file_name
  character(16) name

  logical err

  namelist / key_bindings / key

! Init lattaces
! read global structure from tao_params namelist
! nu(i) keeps track of the sizes of allocatable pointers in universe s%u(i).

  key%ele_name = ' '
  key%universe = '0'
  key%merit_type = s%global%default_key_merit_type

  call tao_open_file ('TAO_INIT_DIR', single_mode_file, iu, file_name)
  read (iu, nml = key_bindings, iostat = ios)
  if (ios < 0) then
    call out_io (s_blank$, r_name, 'Init: No key_bindings namelist found')
    return 
  elseif (ios > 0) then
    call out_io (s_error$, r_name, 'KEY_BINDINGS NAMELIST READ ERROR.')
    return 
  endif
  call out_io (s_blank$, r_name, 'Init: Read key_bindings namelist')
  close (iu)

! associate keys with elements and attributes.

  i_max = 0
  key_loop: do i = 1, size(key)
    if (key(i)%ele_name == ' ') cycle
    i_max = i
    call str_upcase (key(i)%ele_name,    key(i)%ele_name)
    call str_upcase (key(i)%attrib_name, key(i)%attrib_name)
    call str_upcase (key(i)%universe,     key(i)%universe)
    if (key(i)%universe(1:2) == 'U:') key(i)%universe = key(i)%universe(3:)
  enddo key_loop

! allocate the key table. Make it at least s%global%n_key_table_max big.

  i_max = max(s%global%n_key_table_max, i_max)
  allocate (s%key(i_max))
  n1 = s%n_var_used + 1
  s%n_var_used = s%n_var_used + i_max
  n2 = s%n_var_used
  s%var(n1:n2)%exists = .false.
  s%key(:)%ix_var = 0

  do i = 1, i_max

    if (key(i)%ele_name == ' ') cycle
    write (name, *) i
    call string_trim (name, name, ix)
    n = i + n1 - 1
    s%var(n)%name        = 'key-' // name
    s%var(n)%alias       = name
    s%var(n)%ele_name    = key(i)%ele_name
    s%var(n)%attrib_name = key(i)%attrib_name
    s%var(n)%weight      = key(i)%weight
    s%var(n)%step        = key(i)%small_step  
    s%var(n)%high_lim    = key(i)%high_lim
    s%var(n)%low_lim     = key(i)%low_lim
    s%var(n)%good_user   = key(i)%good_opt
    s%var(n)%merit_type  = key(i)%merit_type
    s%var(n)%exists      = .true.

    read (key(i)%universe, *, iostat = ios) ix_u
    if (ios /= 0) then
      call out_io (s_abort$, r_name, 'BAD UNIVERSE INDEX FOR KEY: ' // key(i)%ele_name)
      call err_exit
    endif

    if (ix_u == 0) then
      allocate (s%var(n)%this(size(s%u)))
      do j = 1, size(s%u)
        call tao_pointer_to_var_in_lattice (s%var(n), s%var(n)%this(j), j)
      enddo
    else
      allocate (s%var(n)%this(1))
      call tao_pointer_to_var_in_lattice (s%var(n), s%var(n)%this(1), ix_u)
    endif

    s%var(n)%model_value = s%var(n)%this(1)%model_ptr
    s%var(n)%design_value = s%var(n)%model_value
    s%var(n)%base_value = s%var(n)%this(1)%base_ptr
    s%var(n)%exists = .true.


    s%key(i)%val0 = s%var(n)%model_value
    s%key(i)%delta = key(i)%delta
    s%key(i)%ix_var = n

  enddo

! setup s%v1_var

  nn = s%n_v1_var_used + 1
  s%n_v1_var_used = nn
  call tao_point_v1_to_var (s%v1_var(nn), s%var(n1:n2), 1, n1)
  s%v1_var(nn)%name = 'key'
  s%var(n1:n2)%good_opt = .true.

end subroutine
