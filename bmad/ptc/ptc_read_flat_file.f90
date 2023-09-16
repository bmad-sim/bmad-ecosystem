!+
! Subroutine ptc_read_flat_file (flat_file, err_flag, lat, create_end_marker, from_mad)
!
! Routine to read a PTC "flat file" (or two linked flat files) and optionally create 
! a Bmad lattice based upon the flat file(s).
!
! The information from the flat file is stored in the PTC universes M_u (and possibly M_t).
!
! Input:
!   flat_file(:)      -- character(*): Name(s) of PTC flat file(s).
!   create_end_marker -- logical, optional: Put a marker element named END at the end of the lattice brances?
!                          Default is True.
!   from_mad      -- logical, optional: If True, ignore PTC specific parameters like integrator_order.
!                      Default is False. True is used when the fibre has been created via MAD. In this
!                      case, the PTC specific parameters may not have good values.
!
! Output:
!   err_flag      -- logical: Set True if there is a problem.
!   lat           -- lat_struct, optional: If present then setup a Bmad lattice.
!-

subroutine ptc_read_flat_file (flat_file, err_flag, lat, create_end_marker, from_mad)

use ptc_interface_mod, except_dummy => ptc_read_flat_file
use madx_ptc_module

implicit none

type (lat_struct), optional, target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (fibre), pointer :: fib

integer iu, ios, i, ix_ele, n

logical err_flag
logical, optional :: create_end_marker, from_mad
logical old_style

character(*) flat_file(:)
character(*), parameter :: r_name = 'ptc_read_flat_file'
character(100) line

! Init ptc needed?

if (.not. associated (M_u)) call ptc_ini_no_append

! Old or new format?

err_flag = .true.

iu = lunget()
open(iu, file = flat_file(1), iostat = ios, action = 'read')
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // flat_file(1))
  return
endif

read (iu, '(a)') line
read (iu, '(a)') line
close (iu)

old_style = (index(line, 'high') == 0)

!

call set_ptc_quiet(11, set$, n)

if (size(flat_file) == 1) then
  call read_lattice_append(M_u, flat_file(1))

else  ! New style, two files
  call read_universe_database (M_u, flat_file(1), arpent = my_false)
  call read_universe_pointed (M_u, M_t, flat_file(2))
  call create_dna (M_u, M_t)
endif

call set_ptc_quiet(11, unset$, n)

err_flag = .false.
if (.not. present(lat)) return

! Fill in the lattice

if (size(flat_file) == 1) then
  call init_lat(lat, M_u%start%n + 1)
  branch => lat%branch(0)
  branch%param%geometry = open$

  fib => M_u%start%start
  branch%param%particle = species_of(1d9*fib%mass, nint(fib%charge))
  branch%param%default_tracking_species = ref_particle$

  ix_ele = 0
  do i = 1, M_u%start%n
    if (i /= 1) fib => fib%next
    ! $START and $END are marker elements inserted by mad-x
    if (index(fib%mag%name, '$START') /= 0) cycle
    if (index(fib%mag%name, '$END') /= 0) cycle

    call fibre_to_ele(fib, branch, ix_ele, err_flag, from_mad)
    if (err_flag) return

!    print *
!    print *, '!----------------------------------------------------------'
!    print *, i, '  ', fib%mag%name
!    call type_ptc_fibre (fib, .true.)
  enddo

  if (logic_option(.true., create_end_marker)) then
    ix_ele = ix_ele + 1
    ele => lat%ele(ix_ele)
    call init_ele (ele, marker$, -1, ix_ele, branch)
    ele%name = 'END'
  endif

  branch%ele(0)%value(E_tot_start$) = branch%ele(1)%value(E_tot_start$)
  branch%ele(0)%value(E_tot$)       = branch%ele(1)%value(E_tot$)      
  branch%ele(0)%value(p0c_start$)   = branch%ele(1)%value(p0c_start$)  
  branch%ele(0)%value(p0c$)         = branch%ele(1)%value(p0c$)        

  branch%n_ele_track = ix_ele
  branch%n_ele_max = ix_ele

  lat%input_file_name = flat_file(1)
  lat%use_name = 'from_flat_file'

  call lattice_bookkeeper(lat)
  call lat_sanity_check(lat, err_flag)

  branch%param%total_length = branch%ele(ix_ele)%value(l$)

else
  print *, 'Multiple flat file parsing not yet implemented!'
  stop
endif

end subroutine ptc_read_flat_file

