!+
! Subroutine tao_init_global_and_universes (init_file, data_file, var_file)
!
! Subroutine to initialize the tao structures.
! If init_file, data_file or var_file is not in the current directory then it 
! will be searched for in the directory:
!   TAO_INIT_DIR
!
! Input:
!   init_file      -- Character(*): Tao initialization file.
!   data_file      -- Character(*): Tao data initialization file.
!   var_file       -- Character(*): Tao variable initialization file.
!
! Output:
!-

subroutine tao_init_global_and_universes (init_file, data_file, var_file)

use tao_mod
use tao_data_mod
use tao_lattice_calc_mod
use tao_input_struct
use macroparticle_mod
use bmad_parser_mod
use random_mod
use csr_mod, only: csr_com
use spin_mod
  
implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_input) d2_data
type (tao_d1_data_input) d1_data
type (tao_data_input) data(n_data_minn:n_data_maxx) ! individual weight 
type (tao_v1_var_input) v1_var
type (tao_var_struct), pointer :: var_ptr
type (tao_this_var_struct), pointer :: this
type (tao_var_input) var(n_var_minn:n_var_maxx)
type (tao_global_struct), save :: global, default_global
type (tao_d1_data_struct), pointer :: d1_ptr
type (beam_init_struct) beam_init
type (macro_init_struct) macro_init
type (tao_coupled_uni_input) coupled
type (spin_polar_struct) spin

real(rp) :: default_weight        ! default merit function weight
real(rp) :: default_step          ! default "small" step size
real(rp) default_low_lim, default_high_lim
real(rp) :: default_data_noise         ! default noise for data type

integer ios, iu, i, j, i2, j2, k, ix, n_uni
integer n_data_max, n_var_max, n_d2_data_max, n_v1_var_max
integer n, n_universes, iostat, ix_universe
integer ix_min_var, ix_max_var, n_d1_data
integer ix_min_data, ix_max_data, ix_d1_data
integer, parameter :: ele_name$ = 1, ele_key$ = 2

character(*) init_file, data_file, var_file
character(40) :: r_name = 'tao_init_global_and_universes'
character(200) file_name
character(40) name,  universe, default_universe, default_data_type
character(40) default_merit_type, default_attribute
character(100) line

logical err, free
logical searching
logical calc_emittance
logical, allocatable :: found_one(:), mask(:)


namelist / tao_params / global, bmad_com, csr_com, &
          n_data_max, n_var_max, n_d2_data_max, n_v1_var_max, spin
  
namelist / tao_coupled_uni_init / ix_universe, coupled
  
namelist / tao_beam_init / ix_universe, calc_emittance, beam_init
         
namelist / tao_macro_init / ix_universe, calc_emittance, macro_init
         
namelist / tao_d2_data / d2_data, n_d1_data, default_merit_type, universe, &
                           default_data_noise
  
namelist / tao_d1_data / d1_data, data, ix_d1_data, ix_min_data, &
                           ix_max_data, default_weight, default_data_type
                     
namelist / tao_var / v1_var, var, default_weight, default_step, &
                    ix_min_var, ix_max_var, default_universe, default_attribute, &
                    default_low_lim, default_high_lim

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Init lattaces
! read global structure from tao_params namelist

global = default_global         ! establish defaults
global%valid_plot_who(1:5) = (/ 'model ', 'base  ', 'ref   ', 'design', 'meas  ' /)
global%default_key_merit_type = 'limit'

call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)
if (iu == 0) then
  call out_io (s_abort$, r_name, "Error opening init file")
  call err_exit
endif

read (iu, nml = tao_params)
call out_io (s_blank$, r_name, 'Init: Read tao_params namelist')

close (iu)

s%global = global  ! transfer global to s%global
  
n = size(s%u)
do i = 1, size(s%u)
  call init_universe (s%u(i))
  s%u(i)%ix_uni = i
  s%u(i)%do_synch_rad_int_calc = .false.
  s%u(i)%do_chrom_calc         = .false.
enddo

if (associated(s%var)) deallocate (s%var)
if (associated(s%v1_var)) deallocate (s%v1_var)
allocate (s%var(n_var_max))
allocate (s%v1_var(n_v1_var_max))

s%v1_var%name = ' '  ! blank name means it doesn't (yet) exist
s%var(:)%good_opt  = .true.
s%var(:)%exists    = .false.
s%var(:)%good_var  = .true.
s%var(:)%good_user = .true.

s%n_var_used = 0
s%n_v1_var_used = 0       ! size of s%v1_var(:) array

!-----------------------------------------------------------------------
! Seed random number generator

call ran_seed (s%global%random_seed)

!-----------------------------------------------------------------------
! allocate lattice coord_structs and equate model and base to design

call init_lattices ()
  
!-----------------------------------------------------------------------
! Init coupled universes

call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)

! defaults
do i = 1, size(s%u)
  s%u(i)%coupling%coupled = .false.
  s%u(i)%coupling%match_to_design = .false.
  s%u(i)%coupling%use_coupling_ele = .false.
  s%u(i)%coupling%from_uni = -1
  s%u(i)%coupling%from_uni_s = -1
  s%u(i)%coupling%from_uni_ix_ele = -1
enddo

do
  ix_universe = -1
  coupled%from_universe = -1
  coupled%at_element = ' '
  coupled%at_ele_index = -1
  coupled%at_s = -1
  coupled%match_to_design = .false.

  read (iu, nml = tao_coupled_uni_init, iostat = ios)

  if (ios == 0) then
    if (ix_universe == -1) then
      call out_io (s_abort$, r_name, &
            'INIT: READ TAO_COUPLED_UNI_INIT NAMELIST HAS NOT SET IX_UNIVERSE!')
      call err_exit
    endif
    call out_io (s_blank$, r_name, &
        'Init: Read tao_coupled_uni_init namelist for universe \i3\ ', ix_universe)
    i = ix_universe
    call init_coupled_uni (s%u(i), coupled, i)
    cycle
  elseif (ios > 0) then
    call out_io (s_abort$, r_name, 'INIT: TAO_COUPLED_UNI_INIT NAMELIST READ ERROR!')
    rewind (iu)
    do
      read (iu, nml = tao_coupled_uni_init)  ! generate an error message
    enddo
  endif

  close (iu)
  exit

enddo


!-----------------------------------------------------------------------
! Init Beam

! Do not initialize both beam and macro
if (s%global%track_type == 'beam') then

  call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
  call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)

  ! defaults
  do 
    ix_universe = -1
    beam_init%a_norm_emitt  = -1
    beam_init%b_norm_emitt  = 0.0
    beam_init%dPz_dz = 0.0
    beam_init%center(:) = 0.0
    beam_init%bunch_charge = 0.0
    beam_init%ds_bunch = 1
    beam_init%sig_z   = 0.0
    beam_init%sig_e   = 0.0
    beam_init%renorm_center = .true.
    beam_init%renorm_sigma = .true.
    beam_init%n_bunch = 1
    beam_init%n_particle  = 1
    calc_emittance = .false.
    read (iu, nml = tao_beam_init, iostat = ios)

    if (ios == 0) then
      if (beam_init%a_norm_emitt == -1) then
        call out_io (s_abort$, r_name, &
              'TAO_BEAM_INIT NAMELIST: BEAM_INIT%A_NORM_EMITT NOT SET !')
        call err_exit
      endif
      call out_io (s_blank$, r_name, &
              'Init: Read tao_beam_init namelist for universe \i3\ ', ix_universe)
      if (ix_universe == -1) then
        do i = 1, size(s%u)
          call init_beam(s%u(i), beam_init, calc_emittance)
        enddo
      else
        i = ix_universe
        call init_beam(s%u(i), beam_init, calc_emittance)
      endif
      cycle
    elseif (ios > 0) then
      call out_io (s_abort$, r_name, 'INIT: TAO_BEAM_INIT NAMELIST READ ERROR!')
      rewind (iu)
      do
        read (iu, nml = tao_beam_init)  ! generate an error message
      enddo
    endif

    close (iu)
    exit

  enddo


!-----------------------------------------------------------------------
! Init macroparticles
 
elseif(s%global%track_type == 'macro') then
  call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
  call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)

  ! defaults
  do
    ix_universe = -1
    macro_init%x%norm_emit  = 0.0
    macro_init%y%norm_emit  = 0.0
    macro_init%dPz_dz = 0.0
    macro_init%center(:) = 0.0
    macro_init%ds_bunch = 1
    macro_init%sig_z   = 10e-6
    macro_init%sig_e   = 10e-3
    macro_init%sig_e_cut = 3
    macro_init%sig_z_cut = 3
    macro_init%n_bunch = 1
    macro_init%n_slice = 1
    macro_init%n_macro = 1
    macro_init%n_part  = 1e10
    calc_emittance = .false.
    read (iu, nml = tao_macro_init, iostat = ios)
    if (ios == 0) then
      if (ix_universe == -1) then
        call out_io (s_abort$, r_name, &
                'INIT: READ TAO_MACRO_INIT NAMELIST HAS NOT SET IX_UNIVERSE!')
        call err_exit
      endif
      call out_io (s_blank$, r_name, &
              'Init: Read tao_macro_init namelist for universe \i3\ ', ix_universe)
      i = ix_universe
      call init_macro(s%u(i), macro_init, calc_emittance)  ! generate an error message
      cycle
    elseif (ios > 0) then
      call out_io (s_abort$, r_name, 'INIT: TAO_MACRO_INIT NAMELIST READ ERROR!')
      rewind (iu)
      do
        read (iu, nml = tao_macro_init)  ! generate an error message
      enddo
    endif

    close (iu)
    exit

  enddo

endif

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Init data

call tao_open_file ('TAO_INIT_DIR', data_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)

allocate (mask(size(s%u)))
      
do 
  mask(:) = .true.      ! set defaults
  d2_data%name = ' '
  universe = '*'
  default_merit_type = 'target'
  default_data_noise = 0.0
  read (iu, nml = tao_d2_data, iostat = ios, err = 9100)
  if (ios < 0) exit         ! exit on end-of-file
  call out_io (s_blank$, r_name, &
                      'Init: Read tao_d2_data namelist: ' // d2_data%name)
    
  if (universe == '*') then
    uni_loop1: do i = 1, size(s%u)

    ! check if this data type has already been defined for this universe
    do k = 1, size(s%u(i)%d2_data)
      if (trim(s%u(i)%d2_data(k)%name) == trim(d2_data%name)) then
        mask(i) = .false.
        cycle uni_loop1
      endif
    enddo
      
    call d2_data_stuffit (s%u(i), i)
    enddo uni_loop1
  else
    read (universe, *, iostat = ios) n_uni
    if (ios /= 0 .or. n_uni > size(s%u)) then
      call out_io (s_abort$, r_name, &
            'BAD UNIVERSE NUMBER IN TAO_D2_DATA NAMELIST: ' // d2_data%name)
      call err_exit
    endif
    if (n_uni == 0) then
      call out_io (s_error$, r_name, &
              '"UNIVERSE == 0" MUST BE REPLACED BY "UNIVERSE = *"', &
              ' IN TAO_D2_DATA NAMELIST: ' // d2_data%name)
      call err_exit
    endif
    call d2_data_stuffit (s%u(n_uni), n_uni)
  endif

  do k = 1, n_d1_data
    default_weight = 0      ! set default
    default_data_type  = ' '
    data(:)%data_type  = default_data_type
    data(:)%merit_type = default_merit_type 
    data(:)%name       = ' '
    data(:)%merit_type = ' '
    data(:)%ele_name   = ' '
    data(:)%ele0_name  = ' '
    data(:)%meas_value = real_garbage$  ! used to tag when %meas_value is set in file
    data(:)%weight     = 0.0
    data(:)%data_noise  = real_garbage$
    data(:)%good_user  = .true.
    read (iu, nml = tao_d1_data, err = 9150)
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
    do i = lbound(data, 1), ubound(data, 1)
      if (data(i)%ele0_name /= ' ' .and. data(i)%ele_name == ' ') then
        write (line, '(4a, i0, a)') trim(d2_data%name), '.', trim(d1_data%name), '[', i, ']'
        call out_io (s_abort$, r_name, &
              'ERROR: ELE_NAME IS BLANK BUT ELE0_NAME IS NOT FOR: ' // line)
        call err_exit
      endif
    enddo
    call out_io (s_blank$, r_name, &
                      'Init: Read tao_d1_data namelist: ' // d1_data%name)
    if (universe == '*') then          ! * => use all universes
      uni_loop2: do i = 1, size(s%u)

      ! check if this data type has already been defined for this universe
      if (.not. mask(i)) cycle uni_loop2
      
        call d1_data_stuffit (k, s%u(i), s%u(i)%n_d2_data_used)
      enddo uni_loop2
    else
      call d1_data_stuffit (k, s%u(n_uni), s%u(n_uni)%n_d2_data_used)
    endif
  enddo

enddo

if (allocated(mask)) deallocate(mask)
close (iu)

!-----------------------------------------------------------------------
! Init vars

call tao_open_file ('TAO_INIT_DIR', var_file, iu, file_name)
call out_io (s_blank$, r_name, '*Init: Opening File: ' // file_name)

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
  var%good_user = .true.

  read (iu, nml = tao_var, iostat = ios, err = 9200)
  if (ios < 0) exit         ! exit on end-of-file
  call out_io (s_blank$, r_name, &
                        'Init: Read tao_var namelist: ' // v1_var%name)
  call str_upcase (default_attribute, default_attribute)
  do i = lbound(var, 1), ubound(var, 1)
    call str_upcase (var(i)%attribute, var(i)%attribute)
    call str_upcase (var(i)%ele_name, var(i)%ele_name)
  enddo

  if (v1_var%name == ' ') cycle

  if (default_universe == ' ' .and. all(var%universe == ' ')) &
                                                  default_universe = 'gang'

  if (default_universe == 'clone') then
    do i = 1, size(s%u)
      call var_stuffit_common
      write (s%v1_var(s%n_v1_var_used)%name, '(2a, i0)') &
                                trim(s%v1_var(s%n_v1_var_used)%name), '_u', i
      call var_stuffit_1_uni (i)
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
    call var_stuffit_1_uni (n)
  endif

enddo

close (iu)

!-----------------------------------------------------------------------
! Init ix_data array

do i = 1, size(s%u)
  call init_ix_data (s%u(i))
enddo

return

!-----------------------------------------------------------------------
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
integer i

!

u%is_on = .true.          ! turn universe on
u%n_d2_data_used = 0      ! size of s%u(i)%d2_data(:) array
u%n_data_used = 0         ! size of s%u(i)%data(:) array
u%ix_rad_int_cache = 0

! allocate and set defaults

if (n_d2_data_max /= 0) then
  if (associated(u%d2_data)) deallocate (u%d2_data)
  allocate (u%d2_data(n_d2_data_max))
  do i = 1, n_d2_data_max
    u%d2_data(i)%descrip = ' '
  enddo
  u%d2_data%name = ' '  ! blank name means it doesn't exist
endif

if (n_data_max /= 0) then
  if (associated(u%data)) deallocate (u%data)
  allocate (u%data(n_data_max))
  u%data(:)%exists = .false.       ! set default
  u%data(:)%good_meas  = .false.   ! set default
  u%data(:)%good_ref   = .false.   ! set default
  u%data(:)%good_user  = .true.    ! set default
  u%data(:)%good_opt   = .true.
  u%data(:)%merit_type = 'target'  ! set default
  u%data(:)%ele_name   = ' '
  u%data(:)%ix_ele     = -1
  u%data(:)%ele0_name  = ' '
  u%data(:)%ix_ele0    = 0 ! by default, data relative to beginning of lattice
endif

! This is needed to keep the totalview debugger happy.

if (associated(u%dmodel_dvar)) deallocate (u%dmodel_dvar)
allocate (u%dmodel_dvar(1,1))
  
end subroutine init_universe

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine init_lattices ()

implicit none

integer i

do i = 1, size(s%u)
  n = s%u(i)%design%lat%n_ele_max
  if (allocated(s%u(i)%model%orb)) then
    deallocate (s%u(i)%model%orb, s%u(i)%model%bunch_params)
    deallocate (s%u(i)%design%orb, s%u(i)%design%bunch_params)
    deallocate (s%u(i)%base%orb, s%u(i)%base%bunch_params)
  endif
  allocate (s%u(i)%model%orb(0:n), s%u(i)%model%bunch_params(0:n))
  allocate (s%u(i)%design%orb(0:n), s%u(i)%design%bunch_params(0:n))
  allocate (s%u(i)%base%orb(0:n), s%u(i)%base%bunch_params(0:n))
  ! Specify initial conditions
  s%u(i)%design%orb(0)%vec = 0.0
  call polar_to_spinor (spin, s%u(i)%design%orb(0))
  call init_ring (s%u(i)%model%lat, s%u(i)%design%lat%n_ele_max)
  call init_ring (s%u(i)%base%lat, s%u(i)%design%lat%n_ele_max)
  s%u(i)%model%lat = s%u(i)%design%lat
  s%u(i)%base%lat  = s%u(i)%design%lat
enddo
  
end subroutine init_lattices


!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine d2_data_stuffit (u, ix_uni)

type (tao_universe_struct), target :: u

integer nn, ix_uni

! Setup another d2_data structure.

u%n_d2_data_used = u%n_d2_data_used + 1
nn = u%n_d2_data_used

if (size(u%d2_data) < nn) then
  call out_io (s_error$, r_name, &
              'N_D2_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
  call err_exit
endif

u%d2_data(nn)%name = d2_data%name 
u%d2_data(nn)%ix_uni = ix_uni

! allocate memory for the u%d1_data structures

if (associated(u%d2_data(nn)%d1)) deallocate (u%d2_data(nn)%d1)
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
character(40) search_string
character(20) fmt

!

u%d2_data(n_d2)%d1(i_d1)%d2 => u%d2_data(n_d2)  ! point back to the parent
u%d2_data(n_d2)%d1(i_d1)%name = d1_data%name    ! stuff in the data

! Check if we are searching for elements or repeating elements
! and record the element names in the data structs.
    
if (data(0)%ele_name(1:6) == 'SEARCH') then
  allocate (found_one(u%design%lat%n_ele_max))
  if (data(0)%ele_name(1:11) == 'SEARCH_KEY:') then
    call string_trim(data(0)%ele_name(12:), search_string, ix)
    call find_elements (u, search_string, ele_key$, data(0)%ele0_name, found_one)
  elseif  (data(0)%ele_name(1:7) == 'SEARCH:') then
    call string_trim(data(0)%ele_name(8:), search_string, ix)
    call find_elements (u, search_string, ele_name$, data(0)%ele0_name, found_one)
  else
    call out_io (s_abort$, r_name, 'Syntax Error in data(0)%ele_name SEARCH string')
    call err_exit
  endif
  ! finish finding data array limits
  n1 = u%n_data_used + 1
  n2 = u%n_data_used + count(found_one)
  u%n_data_used = n2
  ix1 = ix_min_data
  ix2 = ix_min_data + (n2 - n1)
  if (n2 > size(u%data)) then
    call out_io (s_abort$, r_name, &
                  'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif
  ! get element names
  jj = n1
  do j = 1, size(found_one)
    if (found_one(j)) then
      if (jj .gt. n2) then
        call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
        call err_exit
      endif
      u%data(jj)%ele_name = u%design%lat%ele_(j)%name
      u%data(jj)%ix_ele   = j
      u%data(jj)%exists   = .true.
      jj = jj + 1
    endif
  enddo
  u%data(n1:n2)%meas_value = 0 
  u%data(n1:n2)%data_type  = default_data_type
  u%data(n1:n2)%merit_type = default_merit_type 
  u%data(n1:n2)%good_meas  = .false.

! SAME:

elseif (data(0)%ele_name(1:5) == 'SAME:') then
  call string_trim (data(0)%ele_name(6:), name, ix)
  call tao_find_data (err, name, d1_ptr = d1_ptr, ix_uni = u%ix_uni)
  if (err .or. .not. associated(d1_ptr)) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
    call err_exit
  endif
  n1 = u%n_data_used + 1
  n2 = u%n_data_used + size(d1_ptr%d)
  u%n_data_used = n2
  ix1 = ix_min_data
  ix2 = ix_min_data + (n2 - n1)
  if (n2 > size(u%data)) then
    call out_io (s_abort$, r_name, &
                'N_DATA_MAX NOT LARGE ENOUGH IN INPUT FILE: ' // file_name)
    call err_exit
  endif

  u%data(n1:n2)%ele_name   = d1_ptr%d%ele_name
  u%data(n1:n2)%ix_ele     = d1_ptr%d%ix_ele
  u%data(n1:n2)%ele0_name  = d1_ptr%d%ele0_name
  u%data(n1:n2)%ix_ele0    = d1_ptr%d%ix_ele0
  u%data(n1:n2)%exists     = d1_ptr%d%exists

  u%data(n1:n2)%merit_type = default_merit_type
  u%data(n1:n2)%weight     = default_weight
  u%data(n1:n2)%data_type  = default_data_type
  u%data(n1:n2)%meas_value = 0 

! Not SEARCH or SAME:

else

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

  ! Transfer info from the input structure

  u%data(n1:n2)%data_type  = data(ix1:ix2)%data_type
  u%data(n1:n2)%merit_type = data(ix1:ix2)%merit_type
  u%data(n1:n2)%good_user  = data(ix1:ix2)%good_user
  u%data(n1:n2)%weight     = data(ix1:ix2)%weight
  u%data(n1:n2)%ele_name   = data(ix1:ix2)%ele_name
  u%data(n1:n2)%ele0_name  = data(ix1:ix2)%ele0_name

  ! Find elements associated with the data

  do j = n1, n2

    call tao_hook_does_data_exist (u%data(j))
    if (u%data(j)%exists) cycle

    if (u%data(j)%ele_name == ' ') cycle
    call str_upcase (u%data(j)%ele_name, u%data(j)%ele_name)
    call element_locator (u%data(j)%ele_name, u%design%lat, ix)
    if (ix < 0) then
      call out_io (s_abort$, r_name, 'ELEMENT NOT LOCATED: ' // &
                                                       u%data(j)%ele_name)
      call err_exit
    endif
    u%data(j)%ix_ele = ix
    u%data(j)%exists = .true.

    if (u%data(j)%ele0_name == ' ') cycle
    call str_upcase (u%data(j)%ele0_name, u%data(j)%ele0_name)
    call element_locator (u%data(j)%ele0_name, u%design%lat, ix)
    if (ix < 0) then
      call out_io (s_abort$, r_name, 'ELEMENT2 NOT LOCATED: ' // &
                                                       u%data(j)%ele0_name)
      call err_exit
    endif
    u%data(j)%ix_ele0 = ix
  enddo

  ! If %meas_value was set then %good_meas is set to True
  u%data(n1:n2)%meas_value = data(ix1:ix2)%meas_value
  where (u%data(n1:n2)%meas_value == real_garbage$)  ! where %meas_value was set
    u%data(n1:n2)%meas_value = 0  
  elsewhere
    u%data(n1:n2)%good_meas = .true.
  end where

endif

!--------------------------------------------
! use default_data_type if given, if not, auto-generate the data_type
if (default_data_type == ' ') then
  where (u%data(n1:n2)%data_type == ' ') u%data(n1:n2)%data_type = &
                            trim(d2_data%name) // '.' // d1_data%name
else
  where (u%data(n1:n2)%data_type == ' ') u%data(n1:n2)%data_type = &
                                                    default_data_type
endif


! set data noise (only applicable to all data types)
if (d2_data%name .eq. "bpm") then
  do j = n1, n2
    u%design%lat%ele_(u%data(j)%ix_ele)%r(1,1) = default_data_noise
  enddo
  do j = lbound(data,1), ubound(data,1)
    if (data(j)%data_noise /= real_garbage$) &
      u%design%lat%ele_(u%data(n1+j-ix1)%ix_ele)%r(1,1) = data(j)%data_noise
  enddo
endif                   

! Create data names
! are we counting elements and forming data names?

if (data(0)%name(1:6) == 'COUNT:') then
  call form_count_name (data(0)%name(7:), num_hashes, count_name1, count_name2)
  jj = ix1
  do j = n1, n2
    if (jj .gt. ix2) then
      call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
      call err_exit
    endif
    write(fmt, '(a,i0,a,i0,a)') '(a, I', num_hashes, '.', num_hashes, ', a)'
    write(u%data(j)%name, fmt) trim(count_name1), jj, trim(count_name2)
    jj = jj + 1
  enddo

elseif (data(0)%name(1:5) == 'SAME:') then
  call string_trim (data(0)%name(6:), name, ix)
  call tao_find_data (err, name, d1_ptr = d1_ptr, ix_uni = u%ix_uni)
  if (err) then
    call out_io (s_abort$, r_name, 'CANNOT MATCH "SAME:" NAME: ' // name)
    call err_exit
  endif
  n2 = n1 + size(d1_ptr%d) - 1
  u%data(n1:n2)%name = d1_ptr%d%name
elseif (index(data(0)%name, '*') /= 0) then
  ix = index(data(0)%name, '*')
  do j = n1, n2
    u%data(j)%name = data(0)%name(1:ix-1) // &
                  trim(u%data(j)%ele_name) // trim(data(0)%name(ix+1:))
  enddo
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

! do we need to do the radiation integrals?

do j = n1, n2
  if (u%data(j)%data_type(1:10) == 'emittance.') u%do_synch_rad_int_calc = .true. 
  if (u%data(j)%data_type(1:2) == 'i5') u%do_synch_rad_int_calc = .true. 
  if (u%data(j)%data_type(1:6) == 'chrom.') u%do_chrom_calc = .true.
  if (u%data(j)%data_type(1:10)     == 'emittance.' .or. &
          u%data(j)%data_type(1:6)  == 'chrom.' .or. &
          u%data(j)%data_type(1:13) == 'unstable_ring' .or. &
          u%data(j)%data_type(1:2)  == 'i5' .or. &
          u%data(j)%data_type(1:17) == 'multi_turn_orbit.') then
    u%data(j)%exists = .true.
    if (u%data(j)%ele_name /= ' ') then
      call out_io (s_abort$, r_name, 'DATUM OF TYPE: ' // u%data(j)%data_type, &
                        'CANNOT HAVE AN ASSOCIATED ELEMENT: ' // u%data(j)%ele_name)
      call err_exit
    endif
    cycle
  endif
enddo

end subroutine d1_data_stuffit

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine var_stuffit_1_uni (ix_u_in)

type (tao_var_struct), pointer :: s_var
integer i, j, ix_u, n, ios, ix_u_in, j_save, ie, n_ele

! point the children back to the mother

n = s%n_v1_var_used
j_save = 1
do i = lbound(s%v1_var(n)%v, 1), ubound(s%v1_var(n)%v, 1)
  s_var => s%v1_var(n)%v(i)

  if (allocated(s_var%this)) deallocate (s_var%this)
  if (s_var%ele_name == ' ') then
    allocate (s_var%this(0))
    s_var%exists = .false.
    cycle
  endif

  ! universe to use
  ix_u = ix_u_in
  if (.not. searching) then
    if (var(i)%universe /= ' ') then
      read (var(i)%universe, *, iostat = ios) ix_u
      if (ios /= 0) then
        call out_io (s_abort$, r_name, &
              'CANNOT READ DEFAULT_UNIVERSE INDEX: ' // default_universe, &
              'FOR VARIABLE: ' // v1_var%name)
        call err_exit
      endif
    endif
    call tao_locate_elements (s_var, ix_u, n_ele)
    if (n_ele == 0) then
      call out_io (s_error$, r_name, 'ELEMENT DOES NOT EXIST: ' // s_var%ele_name)
      s_var%exists = .false.
      cycle
    endif
    allocate (s_var%this(n_ele))
    do ie = 1, n_ele
      call tao_pointer_to_var_in_lattice (s_var, s_var%this(ie), ix_u, &
                                      ix_ele = s%u(ix_u)%model%lat%ele_(ie)%ixx)
    enddo
  else ! use found_one array
    found_one_loop: do j = j_save, size(found_one)
      if (found_one(j)) then
        allocate (s_var%this(1))
        call tao_pointer_to_var_in_lattice (s_var, s_var%this(1), ix_u, ix_ele = j)
        j_save = j+1
        exit found_one_loop
      endif
      if (j == size(found_one)) then
        call out_io (s_abort$, r_name, &
                     "Internal error in counting variables")
        call err_exit
      endif
    enddo found_one_loop
  endif
     
  s_var%model_value = s_var%this(1)%model_ptr
  s_var%design_value = s_var%model_value
  s_var%base_value = s_var%this(1)%base_ptr
  s_var%exists = .true.
enddo

end subroutine var_stuffit_1_uni

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains

subroutine var_stuffit_all_uni 

type (tao_var_struct), pointer :: s_var

integer i, j, n1, n2, iu, j_save, n_tot, n_ele, ie
logical err

! point the children back to the mother


n = s%n_v1_var_used
j_save = 1
  
s_loop: do i = lbound(s%v1_var(n)%v, 1), ubound(s%v1_var(n)%v, 1)
  s_var => s%v1_var(n)%v(i)
  if (allocated(s_var%this)) deallocate (s_var%this)
  if (s_var%ele_name == ' ') then
    allocate (s_var%this(0))
    s_var%exists = .false.
    cycle
  endif
  if (.not. searching) then
    n_tot = 0
    do iu = 1, size(s%u)
      call tao_locate_elements (s_var, iu, n_ele)
      s%u(iu)%ixx = n_ele
      n_tot = n_tot + n_ele
    enddo
    if (n_tot == 0) then
      call out_io (s_error$, r_name, 'ELEMENT DOES NOT EXIST: ' // s_var%ele_name)
      s_var%exists = .false.
      cycle
    endif
    allocate (s_var%this(n_tot))
    n_tot = 0
    do iu = 1, size(s%u)
      do ie = 1, s%u(iu)%ixx
        call tao_pointer_to_var_in_lattice (s_var, s_var%this(ie+n_tot), iu, &
                                          ix_ele = s%u(iu)%model%lat%ele_(ie)%ixx)
      enddo
      n_tot = n_tot + s%u(iu)%ixx
    enddo
  else
    found_one_loop: do j = j_save, size(found_one(1:s%u(1)%design%lat%n_ele_max))
      if (found_one(j)) then
        allocate (s_var%this(size(s%u)))
        do iu = 1, size(s%u)
          call tao_pointer_to_var_in_lattice (s_var, s_var%this(iu), iu, ix_ele = j)
        enddo
        j_save = j+1
        exit found_one_loop
      endif
    enddo found_one_loop
  endif
  s_var%model_value = s_var%this(1)%model_ptr
  s_var%design_value = s_var%this(1)%model_ptr
  s_var%base_value = s_var%this(1)%base_ptr
  s_var%exists = .true.
enddo s_loop

end subroutine

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! stuff common to all universes

subroutine var_stuffit_common

character(20) count_name1, count_name2, ix_char
character(20) fmt

integer i, iu, j, jj, nn, n1, n2, ix1, ix2, num_hashes, ix
integer num_ele, ios, ixx1, ixx2


! count number of v1 entries

s%n_v1_var_used = s%n_v1_var_used + 1
nn = s%n_v1_var_used
s%v1_var(nn)%name = v1_var%name

! are we searching for and counting elements?

if (var(0)%ele_name(1:6) == 'SEARCH') then
  searching = .true.
  ! search through all universes specified
  num_ele = 0
  if (default_universe == 'gang' .or. default_universe == 'clone') then
    do iu = 1, size(s%u)
      num_ele = num_ele + s%u(iu)%design%lat%n_ele_max 
    enddo
    if (allocated(found_one)) deallocate (found_one)  
    allocate (found_one(num_ele))
    ! search in all universes
    call search_for_vars (0, found_one) 
  elseif (any(var%universe /= ' ')) then
    call out_io (s_abort$, r_name, &
           "Cannot specify individual universes when searching for variables")
    call err_exit
  else
    read (default_universe, '(i)', iostat = ios) iu
    if (ios /= 0) then
      call out_io (s_abort$, r_name, &
            "DEFAULT_UNIVERSE MUST BE 'gang', 'clone' OR A NUMBER FOR: " // &
            v1_var%name)
      call err_exit
    endif
    if (iu < 1 .or. iu > size(s%u)) then
      call out_io (s_abort$, r_name, &
              'DEFAULT_UNIVERSE NUMBER OUT OF RANGE FOR: ' // v1_var%name)
    endif
    if (allocated(found_one)) deallocate (found_one)  
    allocate (found_one(s%u(iu)%design%lat%n_ele_max))
    call search_for_vars (iu, found_one)
  endif
  n1 = s%n_var_used + 1
  n2 = s%n_var_used + count(found_one)
  ix1 = ix_min_var
  ix2 = (count(found_one) - (1-ix_min_var))
  s%n_var_used = n2
  ! get element names
  ! if num_ele /= 0 then searching through multiple universes
  jj = n1
  if (num_ele /= 0) then
    ixx2 = 1
    do iu = 1, size(s%u)
      ixx1 = ixx2
      ixx2 = ixx1 + s%u(iu)%design%lat%n_ele_max - 1
      do j = 1, size(found_one(ixx1:ixx2))
        if (found_one(ixx1+j-1)) then
          if (jj .gt. n2) then
            call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT SEARCHING")
            call err_exit
          endif
          s%var(jj)%ele_name = s%u(iu)%design%lat%ele_(j)%name
          s%var(jj)%s = s%u(iu)%design%lat%ele_(j)%s
          jj = jj + 1
        endif
      enddo
    enddo
  else  
    do j = 1, size(found_one)
      if (found_one(j)) then
        if (jj .gt. n2) then
          call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT SEARCHING")
          call err_exit
        endif
        s%var(jj)%ele_name = s%u(iu)%design%lat%ele_(j)%name
        s%var(jj)%s = s%u(iu)%design%lat%ele_(j)%s
        jj = jj + 1
      endif
    enddo
  endif

  s%var(n1:n2)%attrib_name = default_attribute
  s%var(n1:n2)%weight = default_weight
  s%var(n1:n2)%step = default_step
  s%var(n1:n2)%merit_type = default_merit_type
  s%var(n1:n2)%low_lim = default_low_lim
  s%var(n1:n2)%high_lim = default_high_lim

! If not searching...

else  
  searching = .false.
  n1 = s%n_var_used + 1
  n2 = s%n_var_used + ix_max_var - ix_min_var + 1
  ix1 = ix_min_var
  ix2 = ix_max_var
 
  s%n_var_used = n2
 
  s%var(n1:n2)%ele_name    = var(ix1:ix2)%ele_name
  s%var(n1:n2)%name        = var(ix1:ix2)%name
  s%var(n1:n2)%good_user   = var(ix1:ix2)%good_user
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
 
! If counting...

if (var(0)%name(1:6) == 'COUNT:') then
  call form_count_name (var(0)%name(7:), num_hashes, count_name1, count_name2)
  ! Create var names
  jj = ix1
  do j = n1, n2
    if (jj .gt. ix2) then
      call out_io (s_abort$, r_name, "INTERNAL ERROR DURING ELEMENT COUNTING")
      call err_exit
    endif
    write(fmt, '(a,i0,a,i0,a)') '(a, I', num_hashes, '.', num_hashes, ', a)'
    write(s%var(j)%name, fmt) trim(count_name1), jj, trim(count_name2)
    jj = jj + 1
  enddo
elseif (index(var(0)%name, '*') /= 0) then
  ix  = index(var(0)%name, '*')
  do j = n1, n2
    s%var(j)%name = var(0)%name(1:ix-1) // &
                              trim(s%var(j)%ele_name) // trim(var(0)%name(ix+1:))
  enddo
endif

! now for some family guidance...
! point the v1_var mother to the appropriate children in the big data array

call tao_point_v1_to_var (s%v1_var(nn), s%var(n1:n2), ix_min_var, n1)

end subroutine

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! This will search each universe for variables matching name or key
!
! The found_one array contains logicals for all universes to be searched
! through

subroutine search_for_vars (uni, found_one)

implicit none

integer uni, jj, ix1, ix2
logical :: found_one(:)
character(40) search_string

ix2 = 0

if (uni == 0) then
  do jj = 1, size(s%u)
    ix1 = ix2 + 1
    ix2 = s%u(jj)%design%lat%n_ele_max
    if (index(var(0)%ele_name, 'SEARCH_KEY:') /= 0) then
      call string_trim(var(0)%ele_name(12:), search_string, ix)
      call find_elements (s%u(jj), search_string, ele_key$, 'no_slaves', found_one(ix1:ix2))
    elseif (index(var(0)%ele_name, 'SEARCH:') /= 0) then 
      call string_trim(var(0)%ele_name(8:), search_string, ix)
      call find_elements (s%u(jj), search_string, ele_name$, 'no_slaves', found_one(ix1:ix2))
    else
      call out_io (s_abort$, r_name, 'Syntax Error in var(0)%ele_name SEARCH string')
      call err_exit
    endif
  enddo
else
  if (index(var(0)%ele_name, 'SEARCH_KEY:') /= 0) then
    call string_trim(var(0)%ele_name(12:), search_string, ix)
    call find_elements (s%u(uni), search_string, ele_key$, 'no_slaves', found_one)
  elseif (index(var(0)%ele_name, 'SEARCH:') /= 0) then 
    call string_trim(var(0)%ele_name(8:), search_string, ix)
    call find_elements (s%u(uni), search_string, ele_name$, 'no_slaves', found_one)
  else
    call out_io (s_abort$, r_name, 'Syntax Error in var(0)%ele_name SEARCH string')
    call err_exit
  endif

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
if (ix == 0) then
  call out_io (s_abort$, r_name, &
        "WHEN USING 'COUNT:' MUST HAVE '#' WILDCARD IN NAME")
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

subroutine find_elements (u, search_string, attribute, restriction, found_one)

type (tao_universe_struct), target :: u
type (ele_struct), pointer :: ele

character(*) search_string, restriction
integer attribute, key, found_key
logical found_one(:)

integer j, n_max
logical no_slaves

!

if (restriction == 'no_lords') then
  n_max = u%design%lat%n_ele_use
  no_slaves = .false.
elseif (restriction == 'no_slaves') then
  n_max = u%design%lat%n_ele_max
  no_slaves = .true.
elseif (restriction == ' ') then
  n_max = u%design%lat%n_ele_max
  no_slaves = .false.  
else
  call out_io (s_abort$, r_name, "BAD SEARCH RESTRICTION: " // restriction)
  call err_exit
endif

!

found_one = .false.

if (attribute == ele_name$) then
  do j = 1, n_max
    ele => u%design%lat%ele_(j)
    if (no_slaves .and. &
            (ele%control_type == super_slave$ .or. ele%control_type == multipass_slave$)) cycle
    if (.not. match_wild(ele%name, search_string)) cycle
    found_one(j) = .true.
  enddo

elseif (attribute == ele_key$) then
  call match_word (search_string, key_name, key)
  if (key < 1) then
    call out_io (s_abort$, r_name, "Ambiguous or non-existant key name: " // search_string)
    call err_exit
  endif

  do j = 1, n_max
    ele => u%design%lat%ele_(j)
    if (ele%key /= key) cycle
    if (no_slaves .and. &
            (ele%control_type == super_slave$ .or. ele%control_type == multipass_slave$)) cycle
    found_one(j) = .true.
  enddo

else 
  !bug in call to subroutine
  call out_io (s_abort$, r_name, "Internal Error in find_elements!")
  call err_exit
endif

end subroutine find_elements

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! Initialize universe coupling
!

subroutine init_coupled_uni (u, coupled, this_uni_index)

implicit none

type (tao_universe_struct) u
type (tao_universe_struct), pointer ::  from_uni
type (tao_coupled_uni_input) coupled 
integer this_uni_index

character(40) ele_name

integer j, ix

!
if (coupled%from_universe .eq. 0 .or. coupled%at_element .eq. "none") then
  u%coupling%coupled = .false.
  return
endif

if (coupled%from_universe .ge. this_uni_index) then
  call out_io (s_abort$, r_name, &
        "A universe can only inject into a universe with a greater universe index")
  call err_exit
endif
  
u%coupling%coupled = .true.
u%coupling%from_uni = coupled%from_universe
from_uni => s%u(coupled%from_universe)

call init_ele (u%coupling%coupling_ele)
u%coupling%coupling_ele%key = match$
u%coupling%match_to_design = coupled%match_to_design
if (u%coupling%match_to_design) u%coupling%use_coupling_ele = .true.
          
! find extraction element
call string_trim (coupled%at_element, ele_name, ix)
if (ix /= 0) then
  if (coupled%at_s /= -1 .or. coupled%at_ele_index /= -1) then
    call out_io (s_error$, r_name, &
        "INIT Coupling: cannot specify an element, it's index or position at same time!")
    call out_io (s_blank$, r_name, "Will use element name.")
  endif
  if (ele_name == "end") then
    u%coupling%from_uni_s  = from_uni%design%lat%ele_(from_uni%design%lat%n_ele_use)%s
    u%coupling%from_uni_ix_ele = from_uni%design%lat%n_ele_use
  else
    ! using element name 
    ! find last element with name
    do j = from_uni%design%lat%n_ele_use, 0, -1
      if (ele_name(1:ix) == trim(from_uni%design%lat%ele_(j)%name)) then
        u%coupling%from_uni_s = from_uni%design%lat%ele_(j)%s
        u%coupling%from_uni_ix_ele = j
        return
      endif
      if (j == 0) then
        call out_io (s_abort$, r_name, &
                    "Couldn't find coupling element in universe \I\ ", &
                     coupled%from_universe)
        call err_exit
      endif
    enddo
  endif
elseif (coupled%at_ele_index /= -1) then
  if (coupled%at_s /= -1) then
    call out_io (s_error$, r_name, &
        "INIT Coupling: cannot specify an element, it's index or position at same time!")
    call out_io (s_blank$, r_name, "Will use element index.")
  endif
    u%coupling%from_uni_s = from_uni%design%lat%ele_(coupled%at_ele_index)%s
    u%coupling%from_uni_ix_ele = coupled%at_ele_index
else
  ! using s position
  if (s%global%track_type /= 'single' ) then
    call out_io (s_abort$, r_name, &
     "Cannot specify arbitrary s position for coupling if not tracking a single particle")
    call err_exit
  endif
  !FIX_ME: get ix_ele for element right before this s position
  u%coupling%from_uni_s = coupled%at_s
endif

end subroutine init_coupled_uni

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! Initialize the beams. Determine which element to track beam to
!

subroutine init_beam (u, beam_init, calc_emittance)

implicit none

type (tao_universe_struct) u
type (beam_init_struct) beam_init
logical calc_emittance

!
  
if (u%design%lat%param%lattice_type == circular_lattice$) then
  call out_io (s_blank$, r_name, "***")
  call out_io (s_blank$, r_name, &
                 "Beam tracking through a circular lattice.")
  call out_io (s_blank$, r_name, &
         "Twiss parameters and initial orbit will be found from the closed orbit.")
  if (calc_emittance) then
    call out_io (s_blank$, r_name, &
                  "Emittance will be found using the radiation integrals.")
  else 
    call out_io (s_blank$, r_name, &
                  "Emittance will be as set in tao_beam_init.")
  endif
  u%macro_beam%calc_emittance = calc_emittance
  call out_io (s_blank$, r_name, "***")
elseif (calc_emittance) then
  call out_io (s_blank$, r_name, "***")
  call out_io (s_warn$, r_name, &
                "Calc_emittance is only applicable to circular lattices!")
  call out_io (s_blank$, r_name, "***")
endif
  
u%beam%beam_init = beam_init
u%design%orb(0)%vec = beam_init%center

! No initialization for a circular lattice
if (u%design%lat%param%lattice_type == circular_lattice$) return
  
! This is just to get things allocated
call init_beam_distribution (u%design%lat%ele_(0), beam_init, u%beam%beam)
if (u%coupling%coupled) &
    call init_beam_distribution (u%design%lat%ele_(0), beam_init, u%coupling%injecting_beam)

end subroutine init_beam

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! Initialize the macroparticles. Determine which element to track beam to
!

subroutine init_macro(u, macro_init, calc_emittance)

implicit none

type (tao_universe_struct) u
type (macro_init_struct) macro_init
logical calc_emittance

!

if (u%design%lat%param%lattice_type == circular_lattice$) then
  call out_io (s_blank$, r_name, "***")
  call out_io (s_blank$, r_name, &
                 "Macroparticle tracking through a circular lattice.")
  call out_io (s_blank$, r_name, &
         "Twiss parameters and initial orbit will be found from the closed orbit.")
  if (calc_emittance) then
    call out_io (s_blank$, r_name, &
                  "Emittance will be found using the radiation integrals.")
  else 
    call out_io (s_blank$, r_name, &
                  "Emittance will be as set in tao_macro_init.")
  endif
  u%macro_beam%calc_emittance = calc_emittance
  call out_io (s_blank$, r_name, "***")
elseif (calc_emittance) then
  call out_io (s_blank$, r_name, "***")
  call out_io (s_warn$, r_name, &
                "Calc_emittance is only applicable to circular lattices!")
  call out_io (s_blank$, r_name, "***")
endif

u%macro_beam%macro_init = macro_init
u%design%orb(0)%vec = macro_init%center

! Don't initialize beams in circular lattice
if (u%design%lat%param%lattice_type == circular_lattice$) return
    
! This is just to get things allocated
call init_macro_distribution (u%macro_beam%beam, macro_init, u%design%lat%ele_(0), .true.)
if (u%coupling%coupled) &
  call init_macro_distribution (u%coupling%injecting_macro_beam, &
                                             macro_init, u%design%lat%ele_(0), .true.)

! keep track of where macros are lost
if (associated (u%macro_beam%ix_lost)) deallocate (u%macro_beam%ix_lost)
allocate (u%macro_beam%ix_lost(macro_init%n_bunch, macro_init%n_slice, macro_init%n_macro))
u%macro_beam%ix_lost(:,:,:) = -1

end subroutine init_macro

!----------------------------------------------------------------
!----------------------------------------------------------------
! contains
!
! Defines what datums to evaluate at each element in specified universe

subroutine init_ix_data (u)

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct), pointer :: data

integer, automatic :: n_data(-2:u%design%lat%n_ele_max)
integer, automatic :: ix_next(-2:u%design%lat%n_ele_max)

integer j, k, ix_ele

n_data(:) = 0

! allocate the ix_data array
if (associated(u%ix_data)) deallocate(u%ix_data)
allocate(u%ix_data(-2:u%design%lat%n_ele_max))

! find number of datums at each element
do j = 1, size(u%data)
  data => u%data(j)
  if (.not. data%exists) cycle
  if (data%data_type(1:17) == 'multi_turn_orbit.') then
    ix_ele = -2
  elseif (data%ix_ele == -1) then
    ix_ele = -1
  elseif (data%ix_ele0 > data%ix_ele) then
    ix_ele = u%model%lat%n_ele_use
  else
    ix_ele = data%ix_ele
  endif
  n_data(ix_ele) = n_data(ix_ele) + 1
enddo
  
! allocate ix_ele array for each element
do j = lbound(u%ix_data, 1), ubound(u%ix_data, 1)
  if (associated(u%ix_data(j)%ix_datum)) deallocate (u%ix_data(j)%ix_datum)
  if (n_data(j) == 0) cycle
  allocate (u%ix_data(j)%ix_datum(n_data(j)))
enddo

! used for keeping track of current datum index in each ix_ele element
ix_next(:) = 1
  
! setup ix_ele array for each element
! This is the point where the datum is evaluated
! if ix_ele0 > ix_ele then there is "wrap around"
do j = 1, size(u%data)
  data => u%data(j)
  if (.not. data%exists) cycle
  if (data%data_type(1:17) == 'multi_turn_orbit.') then
    ix_ele = -2
  elseif (data%ix_ele == -1) then
    ix_ele = -1
  elseif (data%ix_ele0 > data%ix_ele) then
    ix_ele = u%model%lat%n_ele_use
  else
    ix_ele = data%ix_ele
  endif
  u%ix_data(ix_ele)%ix_datum(ix_next(ix_ele)) = j
  ix_next(ix_ele) = ix_next(ix_ele) + 1
enddo

end subroutine init_ix_data

end subroutine tao_init_global_and_universes
