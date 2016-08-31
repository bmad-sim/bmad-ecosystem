!+
! Subroutine ptc_read_flat_file (flat_file, err_flag, lat)
!
! Routine to read a PTC "flat file" (or two linked flat files) and optionally create 
! a Bmad lattice based upon the flat file(s).
!
! The information from the flat file is stored in the PTC universes M_u (and possibly M_t).
!
! Input:
!   flat_file(:)  -- character(*): Name(s) of PTC flat file(s).
!
! Output:
!   err_flag      -- logical: Set True if there is a problem.
!   lat           -- lat_struct, optional: If present then setup a Bmad lattice.
!-

subroutine ptc_read_flat_file (flat_file, err_flag, lat)

use bmad, except_dummy => ptc_read_flat_file
use madx_ptc_module
use ptc_interface_mod

implicit none

type (lat_struct), optional :: lat
type (ele_struct), pointer :: ele
type (fibre), pointer :: fib

integer iu, ios, i, ix_ele

logical err_flag
logical old_style

character(*) flat_file(:)
character(*), parameter :: r_name = 'ptc_read_flat_file'
character(100) line

! Init ptc needed?

if (.not. associated (M_u)) call ptc_ini_no_append

! Old or new format?

err_flag = .true.

iu = lunget()
open(iu, file = flat_file(1), iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // flat_file(1))
  return
endif

read (iu, '(a)') line
read (iu, '(a)') line
close (iu)

old_style = (index(line, 'high') == 0)

!

if (old_style) then
  call read_and_append_virgin_general (M_u, flat_file(1))

elseif (size(flat_file) == 1) then
  call read_lattice_append(M_u, flat_file(1))

else  ! New style, two files
  call read_universe_database (M_u, flat_file(1), arpent = my_false)
  call read_universe_pointed (M_u, M_t, flat_file(2))
  call create_dna (M_u, M_t)
endif

err_flag = .false.
if (.not. present(lat)) return

! Fill in the lattice

if (size(flat_file) == 1) then
  call init_lat(lat, M_u%start%n)

  fib => M_u%start%start
  lat%param%particle = species_of(fib%mass, nint(fib%charge))

  ix_ele = 1
  do i = 1, M_u%start%n
    if (i /= 1) fib => fib%next
    ! $START and $END are marker elements inserted by mad-x
    if (index(fib%mag%name, '$START') /= 0) cycle
    if (index(fib%mag%name, '$END') /= 0) cycle

    call fibre_to_ele(fib, lat%ele, ix_ele, err_flag)
    if (err_flag) return

    print *
    print *, '!----------------------------------------------------------'
    print *, i, '  ', fib%mag%name
    call type_ptc_fibre (fib, .true.)

  enddo

else

endif

end subroutine ptc_read_flat_file

