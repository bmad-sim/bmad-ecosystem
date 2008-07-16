module tao_read_beam_mod

use tao_mod

real(rp), private, save :: bunch_charge_set = 0

integer, private, save :: num_bunch_set = 0, num_particle_set = 0
integer, private, save :: iu

logical, private, save :: ascii_file

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_open_beam_file (file_name)
!
! Routine to open a beam file for reading.
! Call tao_close_beam_file when done.
!
! Modules needed:
!   use tao_read_beam_mod
!
! Input:
!   file_name -- Character(*): Name of the beam file.
!-

subroutine tao_open_beam_file (file_name)

implicit none

character(*) file_name
character(100) :: full_file_name
character(80) line

! Open file and determine whether the file is binary or ascii

call tao_open_file ('TAO_INIT_DIR', file_name, iu, full_file_name)
if (iu == 0) stop

ascii_file = .true.
read (iu, '(a80)') line
if (index(line, '!BINARY') /= 0) ascii_file = .false.

if (ascii_file) then
  rewind (iu)
else
  close (iu)
  open (iu, file = full_file_name, form = 'unformatted', status = 'old')
  read (iu) line(1:7)  ! read "!binary" line
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_close_beam_file ()
!
! Routine to close a beam file that was opened with open_beam_file.
!
! Modules needed:
!   use tao_read_beam_mod
!-

subroutine tao_close_beam_file ()

implicit none

close (iu)

! reset parameters

num_bunch_set = 0
num_particle_set = 0
bunch_charge_set = 0

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_read_beam_file_header (ix_ele, n_bunch, n_particle)
!
! Routine to read in the beam parameters from a file
!
! Modules needed:
!   use tao_read_beam_mod
!
! Output:
!   ix_ele     -- Integer: Index of the lattice element at which the beam is to
!                   be specified. 0 => beginning of lattice. -1 => End of file.
!   n_bunch    -- Integer: Number of bunches.
!   n_particle -- Integer: Number of particles per bunch.
!-


subroutine tao_read_beam_file_header (ix_ele, n_bunch, n_particle)

implicit none

real(rp) charge_bunch
integer n_bunch, n_particle, ix_ele, ios
character(20) :: r_name = 'read_beam_file_header'

! Read numbers

if (ascii_file) then
  read (iu, *, iostat = ios) ix_ele
  read (iu, *, iostat = ios) n_bunch
  read (iu, *, iostat = ios) n_particle
else
  read (iu, iostat = ios) ix_ele, n_bunch, n_particle
endif

if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING HEADER INFO!')
  call err_exit
endif

if (ios < 0) ix_ele = -1  ! End of file 

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_set_beam_params (n_bunch, n_particle, charge_bunch)
!
! Subroutine to set some beam parameters.
! These parameters will override the parameters set in the beam file.
! This subroutine must be called before read_beam.
!
! Modules needed:
!   use tao_read_beam_mod
!
! Input:
!   n_bunch       -- Integer: Number of bunches.
!   n_particle    -- Integer: Number of particles in a bunch.
!   charge_bunch  -- Real(rp): Bunch charge
!-

subroutine tao_set_beam_params (n_bunch, n_particle, charge_bunch)

implicit none

real(rp), optional :: charge_bunch
integer, optional :: n_bunch, n_particle

!

if (present(n_bunch))      num_bunch_set = n_bunch
if (present(n_particle))   num_particle_set = n_particle
if (present(charge_bunch)) bunch_charge_set = charge_bunch

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_read_beam (beam)

implicit none

type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (particle_struct), pointer :: p(:)

integer i, j, k, ix, ios
integer n_bunch, n_particle, n_particle_lines, ix_lost

real(rp) vec(6), sum_charge
complex(rp) spin(2)

character(16) :: r_name = 'read_beam'
character(200) line

! If an ASCII file: Count the number of particle lines.

if (ascii_file) then
  n_particle_lines = -3
  do 
    read (iu, '(a)', iostat = ios) line
    if (ios /= 0) then
      call out_io (s_fatal$, r_name, 'NO "END_BUNCH" MARKER FOUND!')
      call err_exit
    endif
    call string_trim (line, line, ix)
    if (line(1:9) == 'END_BUNCH') exit
    n_particle_lines = n_particle_lines + 1
  enddo  
  rewind (iu)
endif

! Parameters come from the beam file unless set_beam_params has been called.

call tao_read_beam_file_header (i, n_bunch, n_particle)
if (n_particle == 0) n_particle = n_particle_lines

if (num_bunch_set > 0) n_bunch = num_bunch_set
if (num_particle_set > 0) n_particle = num_particle_set

! Allocate space

call reallocate_beam (beam, n_bunch, n_particle)

! An ascii file, if it is generated by another program, may not include ix_lost or the spin.
! so add the default values

do i = 1, n_bunch
  bunch => beam%bunch(i)
  p => bunch%particle

  if (ascii_file) then
    read (iu, *) bunch%charge
    read (iu, *) bunch%z_center
    read (iu, *) bunch%t_center
    do j = 1, n_particle_lines 

      if (j > n_particle) cycle
      p(j)%charge = 0; p(j)%ix_lost = not_lost$; p(j)%r%spin = cmplx(0.0_rp, 0.0_rp)

      read (iu, '(a)') line

      ix = 0
      do k = 1, 6
        call string_trim(line(ix+1:), line, ix)
        read (line, *) p(j)%r%vec(k)
      enddo

      call string_trim(line(ix+1:), line, ix)
      if (ix == 0) cycle
      read (line, *) p(j)%charge

      call string_trim(line(ix+1:), line, ix)
      if (ix == 0) cycle
      read (line, *) p(j)%ix_lost

      call string_trim(line(ix+1:), line, ix)
      if (ix == 0) cycle
      read (line, *) p(j)%r%spin

    enddo
    read (iu, '(a)') line  ! END_BUNCH line

  else
    read (iu) bunch%charge, bunch%z_center, bunch%t_center, n_particle_lines
    do j = 1, n_particle_lines
      read (iu) p(j)%r%vec, p(j)%charge, p(j)%ix_lost, p(j)%r%spin
      if (j > n_particle) cycle
    enddo
  endif

  if (bunch_charge_set /= 0) bunch%charge = bunch_charge_set
  sum_charge = sum(p(:)%charge)
  if (bunch%charge == 0) then
    bunch%charge = sum_charge
  elseif (sum_charge == 0) then
    p%charge = bunch%charge / n_particle
  else
    p%charge = p%charge * bunch%charge / sum_charge
  endif
enddo

end subroutine

end module
