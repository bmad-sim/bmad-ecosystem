module tao_read_beam_mod

use tao_mod

type read_beam_common_struct
  real(rp) :: bunch_charge = 0
  integer :: n_bunch = 0, n_particle = 0
  integer :: iu
  character(100) :: file_name
  logical :: ascii_file
end type

type (read_beam_common_struct), private, save :: rb_com

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

call tao_open_file ('TAO_INIT_DIR', file_name, rb_com%iu, full_file_name)
if (rb_com%iu == 0) stop

rb_com%ascii_file = .true.
read (rb_com%iu, '(a80)') line
if (index(line, '!BINARY') /= 0) rb_com%ascii_file = .false.

if (rb_com%ascii_file) then
  rewind (rb_com%iu)
else
  close (rb_com%iu)
  open (rb_com%iu, file = full_file_name, form = 'unformatted', status = 'old')
endif

rb_com%file_name = full_file_name

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

close (rb_com%iu)

! reset parameters

rb_com%n_bunch = 0
rb_com%n_particle = 0
rb_com%bunch_charge = 0

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_read_beam_file_header (ix_ele, n_bunch, n_particle, err)
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
!   err        -- Logical: Set True if there is an error. False otherwise.
!-


subroutine tao_read_beam_file_header (ix_ele, n_bunch, n_particle, err)

implicit none

real(rp) charge_bunch
integer n_bunch, n_particle, ix_ele, ios
character(80) line
character(20) :: r_name = 'read_beam_file_header'
logical err

! Read numbers

err = .true.

if (rb_com%ascii_file) then
  read (rb_com%iu, *, iostat = ios) ix_ele
  read (rb_com%iu, *, iostat = ios) n_bunch
  read (rb_com%iu, *, iostat = ios) n_particle
else
  read (rb_com%iu) line(1:7)  ! read "!BINARY" line
  read (rb_com%iu, iostat = ios) ix_ele, n_bunch, n_particle
endif

if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING HEADER INFO!')
  return
endif

if (ios < 0) ix_ele = -1  ! End of file 

err = .false.

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

if (present(n_bunch))      rb_com%n_bunch = n_bunch
if (present(n_particle))   rb_com%n_particle = n_particle
if (present(charge_bunch)) rb_com%bunch_charge = charge_bunch

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_read_beam (beam, err)
!
! Subroutine to read in a beam definition file opened by tao_open_beam_file.
! 
! Modules needed:
!   use tao_read_beam_mod
!
! Output:
!   beam -- Beam_struct: Structure hoding the beam information.
!   err  -- Logical: Set True if there is an error. False otherwise.
!+ 


subroutine tao_read_beam (beam, err)

implicit none

type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (particle_struct), pointer :: p(:)

integer i, j, k, ix, ios
integer n_bunch, n_particle, n_particle_lines, ix_lost

real(rp) vec(6), sum_charge
complex(rp) spin(2)

character(16) :: r_name = 'read_beam'
character(200) line, line_in

logical err, error

! If an ASCII file: Count the number of particle lines.

err = .true.

! Parameters come from the beam file unless set_beam_params has been called.

rewind (rb_com%iu)
call tao_read_beam_file_header (i, n_bunch, n_particle, error); if (error) return
n_particle_lines = n_particle

if (rb_com%n_bunch > 0) n_bunch = rb_com%n_bunch
if (rb_com%n_particle > 0) n_particle = rb_com%n_particle

! Allocate space

call reallocate_beam (beam, n_bunch, n_particle)

! An ascii file, if it is generated by another program, may not include ix_lost or the spin.
! so add the default values

do i = 1, n_bunch
  bunch => beam%bunch(i)
  p => bunch%particle

  if (rb_com%ascii_file) then

    read (rb_com%iu, '(a)', iostat = ios) line
    if (ios /= 0 .or. index(line, 'BEGIN_BUNCH') == 0) then
      call out_io (s_error$, r_name, 'NO "BEGIN_BUNCH" MARKER FOUND IN: ' // rb_com%file_name)
      return
    endif

    read (rb_com%iu, *, iostat = ios) bunch%charge
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'ERROR READING BUNCH CHARGE IN: ' //  rb_com%file_name)
      return
    endif

    read (rb_com%iu, *, iostat = ios) bunch%z_center
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'ERROR READING BUNCH Z_CENTER IN: ' // rb_com%file_name)
      return
    endif

    read (rb_com%iu, *, iostat = ios) bunch%t_center
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'ERROR READING BUNCH T_CENTER IN: ' // rb_com%file_name)
      return
    endif

    ! particle coord loop

    j = 0
    do 

      read (rb_com%iu, '(a)', iostat = ios) line
      if (ios /= 0) then
        call out_io (s_error$, r_name, &
                      'ERROR READING PARTICLE COORDINATE LINE: ' // rb_com%file_name)
        return
      endif
      line_in = line ! save for error messages

      j = j + 1

      if (index(line, 'END_BUNCH') /= 0) exit
      if (j > n_particle) cycle

      p(j)%charge = 0; p(j)%ix_lost = not_lost$; p(j)%r%spin = cmplx(0.0_rp, 0.0_rp)

      ix = 0
      do k = 1, 6
        call string_trim(line(ix+1:), line, ix)
        read (line, *, iostat = ios) p(j)%r%vec(k)
        if (ios /= 0) then
          call out_io (s_error$, r_name, &
                        'ERROR READING PARTICLE COORDINATES IN: ' // rb_com%file_name, &
                        'BAD LINE: ' // trim(line_in))
          return
        endif
      enddo

      call string_trim(line(ix+1:), line, ix)
      if (ix == 0) cycle
      read (line, *, iostat = ios) p(j)%charge
      if (ios /= 0) then
        call out_io (s_error$, r_name, &
                        'ERROR READING PARTICLE CHARGE IN: ' // rb_com%file_name, &
                        'BAD LINE: ' // trim(line_in))
        return
      endif

      call string_trim(line(ix+1:), line, ix)
      if (ix == 0) cycle
      read (line, *, iostat = ios) p(j)%ix_lost
      if (ios /= 0) then
        call out_io (s_error$, r_name, &
                        'ERROR READING PARTICLE "IX_LOST" IN: ' // rb_com%file_name, &
                        'BAD LINE: ' // trim(line_in))
        return
      endif

      call string_trim(line(ix+1:), line, ix)
      if (ix == 0) cycle
      read (line, *, iostat = ios) p(j)%r%spin
      if (ios /= 0) then
        call out_io (s_error$, r_name, &
                        'ERROR READING PARTICLE SPIN IN: ' // rb_com%file_name, &
                        'BAD LINE: ' // trim(line_in))
        return
      endif

    enddo

  else
    read (rb_com%iu, iostat = ios) bunch%charge, &
                            bunch%z_center, bunch%t_center, n_particle_lines
    if (ios /= 0) then
      call out_io (s_error$, r_name, &
                        'ERROR READING BUNCH PARAMETERS IN: ' // rb_com%file_name)
      return
    endif

    do j = 1, n_particle_lines
      if (j > n_particle) exit
      read (rb_com%iu, iostat = ios) p(j)%r%vec, p(j)%charge, p(j)%ix_lost, p(j)%r%spin
      if (ios /= 0) then
        call out_io (s_error$, r_name, &
                        'ERROR READING PARTICLE COORDINATES IN: ' // rb_com%file_name)
        return
      endif
    enddo
  endif

  if (rb_com%bunch_charge /= 0) bunch%charge = rb_com%bunch_charge
  sum_charge = sum(p(:)%charge)
  if (bunch%charge == 0) then
    bunch%charge = sum_charge
  elseif (sum_charge == 0) then
    p%charge = bunch%charge / n_particle
  else
    p%charge = p%charge * bunch%charge / sum_charge
  endif
enddo

err = .false.

end subroutine

end module
