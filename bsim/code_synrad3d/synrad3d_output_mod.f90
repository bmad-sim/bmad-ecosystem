module synrad3d_output_mod

use synrad3d_utils

implicit none

contains

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_write_hit_points (file_name, photon, wall_hit, lat, lots_of_digits)
!
! Routine to record the points at which a photon has been reflected.
!
! Input:
!   file_name       -- character(*): Name of the data file.
!   photon          -- sr3d_photon_track_struct: Photon 
!   wall_hit(0:)    -- sr3d_photon_wall_hit_struct: Array of hit point information.
!   lat             -- lat_struct: Lattice.
!   lots_of_digits  -- logical, optional: If true then print using more digits than what is
!                       printed if this argument is False (which is the default).
!   iunit           -- integer, optional: If not present or equal to 0 then create a new file using file_name.
!                       Otherwise  use this file unit number for writing instead of creating a new file.
!
! Output:
!   iunit           -- integer, optional: File unit number used.
!-

subroutine sr3d_write_hit_points (file_name, photon, wall_hit, lat, lots_of_digits, iunit)

type (sr3d_photon_track_struct), target :: photon
type (sr3d_photon_wall_hit_struct), pointer :: hit
type (sr3d_photon_wall_hit_struct), target :: wall_hit(0:)
type (lat_struct) lat

integer, optional :: iunit
integer iu, n, iu_hit_file

logical, optional :: lots_of_digits
logical new_file

character(*) file_name
character(100) fm, fm2
character(40) wall_name

! Open file

if (file_name == '') return

new_file = .true.

if (present(iunit)) then
  if (iunit == 0) then
    iunit = lunget()
    iu = iunit
    open (iu, file = file_name, recl = 500)
  else
    iu = iunit
    new_file = .false.
  endif

else
  iu = lunget()
  open (iu, file = file_name, recl = 500)
endif

!

if (logic_option(.false., lots_of_digits)) then
  fm  = '(6es25.15)'
  fm2 = '(i7, i4, f10.2, 5x, 3es25.15, i4, 2(5x, 3es25.15), 10x, 3f18.12, 5x, 3f16.10, 3x, a)' 
else
  fm  = '(6f12.6)'
  fm2 = '(i7, i4, f10.2, 5x, 2f10.6, f14.6, i4, 2(5x, 3f10.6), 10x, 3f10.6, 5x, 3f10.6, 3x, a)'
  if (new_file) then
    write (iu, '(a)') '#                                                             ix_    |            Before             |             After                 |             Wall Perpendicular                Cos(perp)'
    write (iu, '(a)') '# Index n_hit  Energy         X          Y            S       branch |   Vx        Vy         Vs     |       Vx        Vy        Vs      |            x         y        s             In        Out    Reflectivity'
  endif
endif

!

do n = 0, photon%n_wall_hit
  hit => wall_hit(n)
  if (n == 0) then
    wall_name = '<inside chamber>'
  else
    wall_name = lat%branch(hit%ix_branch)%wall3d(hit%ix_wall3d)%name
    if (wall_name == '') wall_name = '<default_subchamber>'
  endif

  write (iu, fm2) photon%ix_photon, n, hit%before_reflect%p0c, hit%after_reflect%vec(1:3:2), hit%after_reflect%s, &
          hit%ix_branch, hit%before_reflect%vec(2:6:2), hit%after_reflect%vec(2:6:2), &
          hit%dw_perp, hit%cos_perp_in, hit%cos_perp_out, hit%reflectivity, trim(wall_name)
enddo

if (.not. present(iunit)) then
  close (iu)
  print *, 'Written file of hit points: ', trim(file_name)
endif

end subroutine sr3d_write_hit_points

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_write_photon_start_file (file_name, photon, print_message)
!
! Routine to write a namelist file of a photon's starting position that can be used to reinit
! the photon. This is useful for debugging purposes.
!
! Input:
!   file_name       -- character(*): Name of the data file to create.
!   photon          -- sr3d_photon_track_struct: Photon whose initial coords are to be written.
!   print_message   -- logical, optional: Print file creation message? Default is True.
!-

subroutine sr3d_write_photon_start_file (file_name, photon, print_message)

type (sr3d_photon_track_struct) photon

integer iu
logical, optional :: print_message
character(*) file_name

!

iu = lunget()
open (iu, file = file_name, recl = 200)

write (iu, *) 'ix_photon           =', photon%ix_photon
write (iu, *) 'ix_photon_generated =', photon%ix_photon_generated
write (iu, *) 

write (iu, *) '&start'
write (iu, *) '  ran_state%ix                =', sr3d_params%ran_state%ix
write (iu, *) '  ran_state%iy                =', sr3d_params%ran_state%iy
write (iu, *) '  ran_state%am                =', sr3d_params%ran_state%am
write (iu, *) '  ran_state%h_saved           =', sr3d_params%ran_state%h_saved
write (iu, *) '  ran_state%gauss_converter   =', sr3d_params%ran_state%gauss_converter
write (iu, *) '  ran_state%gauss_sigma_cut   =', sr3d_params%ran_state%gauss_sigma_cut
write (iu, *) '  ran_state%seed              =', sr3d_params%ran_state%seed
write (iu, *) '  ran_state%number_stored     =', sr3d_params%ran_state%number_stored
write (iu, *) '  ran_state%engine            =', sr3d_params%ran_state%engine
write (iu, *)
write (iu, *) '  orbit%vec       =', photon%start%orb%vec
write (iu, *) '  orbit%p0c       =', photon%start%orb%p0c
write (iu, *) '  orbit%s         =', photon%start%orb%s
write (iu, *) '  orbit%direction =', photon%start%orb%direction
write (iu, *) '  ix_branch       =', photon%start%ix_branch
write (iu, *) '/'

close (iu)
if (logic_option(.true., print_message)) print *, 'Written photon init file: ', trim(file_name)

end subroutine sr3d_write_photon_start_file

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

subroutine sr3d_record_photon_position (action, photon)

type (sr3d_photon_track_struct), target, optional :: photon
type (coord_struct), pointer :: orb

integer ios
integer, save :: iu = 0
character(*) action
character(120) line

! Only tracks that pass the filter tests get into the official track file.

select case (action)

case ('START_RECORDING')
  ! Open a scratch file
  if (iu /= 0) call err_exit ! Should be zero
  iu = lunget()
  open (iu, status = 'scratch')

case ('ERASE_RECORDING')
  ! Close scratch file
  if (iu == 0) call err_exit ! Should be non-zero
  close (iu)
  iu = 0

case ('MOVE_TRACK_TO_FILE')
  ! Move info from scratch file to the official track file
  if (iu == 0) call err_exit ! Should be non-zero
  rewind (iu)
  do
    read (iu, '(a)', iostat = ios) line
    if (ios /= 0) exit
    write (sr3d_params%iu_photon_track, '(a)') trim(line)
  enddo
  close (iu)
  iu = 0

case ('RECORD_TRACK_POINT')
  ! Record a track point in the scratch file.
  if (iu == 0) call err_exit ! Should be non-zero
  orb => photon%now%orb
  write (iu, '(i8, i10, 2f11.6, f13.6, i4, 5x, 3f11.6)') &
     photon%ix_photon, photon%ix_photon_generated, orb%vec(1:3:2), orb%s, photon%now%ix_branch, orb%vec(2:6:2)

case default
  call err_exit   ! Should not be here
end select

end subroutine sr3d_record_photon_position

end module
