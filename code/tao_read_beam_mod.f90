module tao_read_beam_mod

use tao_mod

type read_beam_common_struct
  real(rp) :: bunch_charge = 0
  integer :: n_bunch = 0, n_particle = 0
  integer :: iu
  character(100) :: file_name
  character(8) file_type
end type

type (read_beam_common_struct), private, save :: rb_com

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_open_beam_file (file_name, err_flag)
!
! Routine to open a beam file for reading.
! Call tao_close_beam_file when done.
!
! Modules needed:
!   use tao_read_beam_mod
!
! Input:
!   file_name -- Character(*): Name of the beam file.
!   err_flag  -- Logical :: Set true if cannot open the file
!-

subroutine tao_open_beam_file (file_name, err_flag)

implicit none

logical err_flag

character(*) file_name
character(100) :: full_file_name
character(80) line
character(20), parameter :: r_name = 'tao_open_beam_file'

! Open file and determine whether the file is binary or ascii

err_flag = .true.

call tao_open_file (file_name, rb_com%iu, full_file_name, s_error$)
if (rb_com%iu == 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN BEAM FILE!')
  return
endif

read (rb_com%iu, '(a80)') line
if (index(line, '!BINARY') /= 0) then 
  rb_com%file_type = 'BIN:1'
elseif (index(line, '!BIN::2') /= 0) then 
  rb_com%file_type = 'BIN:2'
else
  rb_com%file_type = 'ASCII'
endif

if (rb_com%file_type == 'ASCII') then
  rewind (rb_com%iu)
else
  close (rb_com%iu)
  open (rb_com%iu, file = full_file_name, form = 'unformatted', status = 'old')
endif

rb_com%file_name = full_file_name
err_flag = .false.

end subroutine tao_open_beam_file

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

end subroutine tao_close_beam_file

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
character(*), parameter :: r_name = 'tao_read_beam_file_header'
logical err

! Read numbers

err = .true.

if (rb_com%file_type == 'ASCII') then
  read (rb_com%iu, *, iostat = ios, err = 8000) ix_ele
  read (rb_com%iu, *, iostat = ios, err = 8000) n_bunch
  read (rb_com%iu, *, iostat = ios, err = 8000) n_particle
else
  read (rb_com%iu) line(1:7)  ! read "!BINARY" line
  read (rb_com%iu, iostat = ios, err = 8000) ix_ele, n_bunch, n_particle
endif

8000 continue
if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING BEAM HEADER INFO IN FILE: ' // trim(rb_com%file_name))
  return
endif

if (ios < 0) ix_ele = -1  ! End of file 

err = .false.

end subroutine tao_read_beam_file_header

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

end subroutine tao_set_beam_params

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
type (coord_struct), pointer :: p(:)
type (coord_struct) orb_init

integer i, j, k, ix, ix_word, ios, ix_ele, species
integer n_bunch, n_particle, n_particle_lines, ix_lost

real(rp) vec(6), sum_charge
complex(rp) spin(2)

character(16) :: r_name = 'tao_read_beam'
character(300) line, line_in

logical err, error, in_parens

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
  p = orb_init   ! init with default params

  if (rb_com%file_type == 'ASCII') then

    read (rb_com%iu, '(a)', iostat = ios) line
    if (ios /= 0 .or. index(upcase(line), 'BEGIN_BUNCH') == 0) then
      call this_error_out ('NO "BEGIN_BUNCH" MARKER FOUND')
      return
    endif

    read (rb_com%iu, *, iostat = ios) line
    if (is_real(line, .true.)) then
      call this_error_out ('OLD STYLE FORMAT DOES NOT INCLUDE BUNCH SPECIES', 'PLEASE CORRECT.')
      return
    endif

    ix = species_id(line)
    if (ix == invalid$) then
      call this_error_out ('BAD SPECIES NAME: ' // trim(line))
      return
    endif
    bunch%particle%species = ix

    read (rb_com%iu, *, iostat = ios) bunch%charge_tot
    if (ios /= 0) then
      call this_error_out ('ERROR READING BUNCH')
      return
    endif

    read (rb_com%iu, *, iostat = ios) bunch%z_center
    if (ios /= 0) then
      call this_error_out ('ERROR READING BUNCH Z_CENTER')
      return
    endif

    read (rb_com%iu, *, iostat = ios) bunch%t_center
    if (ios /= 0) then
      call this_error_out ('ERROR READING BUNCH T_CENTER')
      return
    endif

    !----------------------------------------
    ! particle coord loop

    j = 0
    do 

      read (rb_com%iu, '(a)', iostat = ios) line
      if (ios /= 0) then
        call this_error_out ('ERROR READING PARTICLE COORDINATE LINE')
        return
      endif
      line_in = line ! save for error messages

      j = j + 1
      in_parens = .false.

      if (index(upcase(line), 'END_BUNCH') /= 0) exit
      if (j > n_particle) cycle

      p(j)%charge = 0; p(j)%state = alive$; p(j)%spin = cmplx(0.0_rp, 0.0_rp)

      call string_trim(line, line, ix_word)
      do k = 1, 6
        read (line, *, iostat = ios) p(j)%vec(k)
        if (ios /= 0) then
          call this_error_out ('ERROR READING PARTICLE COORDINATES', 'IN LINE: ' // trim(line_in))
          return
        endif
        if (.not. remove_first_number (line, ix_word, '', in_parens)) return
      enddo

      read (line, *, iostat = ios) p(j)%charge
      if (ios /= 0 .or. ix_word == 0) then
        call this_error_out ('ERROR READING PARTICLE CHARGE', 'IN LINE: ' // trim(line_in))
        return
      endif
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return

      if (ix_word == 0) goto 8000
      read (line, *, iostat = ios) p(j)%state
      if (ios /= 0) then
        call this_error_out ('ERROR READING PARTICLE "STATE"', 'IN LINE: ' // trim(line_in))
        return
      endif
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return

      if (ix_word == 0) goto 8000
      read (line, *, iostat = ios) p(j)%spin
      if (ios /= 0) then
        call this_error_out ('ERROR READING PARTICLE SPIN', 'IN LINE: ' // trim(line_in))
        return
      endif
      if (.not. remove_first_number (line, ix_word, '(x', in_parens)) return
      if (.not. remove_first_number (line, ix_word, 'x)(', in_parens)) return
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return
      if (.not. remove_first_number (line, ix_word, 'x)', in_parens)) return

      if (ix_word == 0) goto 8000
      read (line, *, iostat = ios) p(j)%ix_ele
      if (ios /= 0) then
        call this_error_out ('ERROR READING ELEMENT INDEX', 'IN LINE: ' // trim(line_in))
        return
      endif
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return

      if (ix_word == 0) goto 8000
      read (line, *, iostat = ios) p(j)%location
      if (ios /= 0) then
        call this_error_out ('ERROR READING PARTICLE LOCATION', 'IN LINE: ' // trim(line_in))
        return
      endif
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return

      8000 continue
      if (in_parens .or. ix_word /= 0) then
        call this_error_out ('UNMATCHED PARENTHESIS IN LINE: ' // trim(line_in))
        return
      endif

    enddo


  !------------------------------------------------------------------------------------
  ! Binary file

  else
    if (rb_com%file_type == 'BIN:1') then
      read (rb_com%iu, iostat = ios) bunch%charge_tot, bunch%z_center, bunch%t_center, n_particle_lines
      bunch%particle%species = electron$
      call this_error_out ('OLD STYLE BEAM0 FILE WITHOUT SPECIES INFO. ASSUMING ELECTRONS...') 
    else
      read (rb_com%iu, iostat = ios) species, bunch%charge_tot, bunch%z_center, bunch%t_center, n_particle_lines
      p%species = species
    endif

    if (ios /= 0) then
      call this_error_out ('ERROR READING BUNCH PARAMETERS')
      return
    endif

    do j = 1, n_particle_lines
      if (j > n_particle) exit
      read (rb_com%iu, iostat = ios) p(j)%vec, p(j)%charge, p(j)%state, p(j)%spin, &
                                     ix_ele, p(j)%location
      if (ios /= 0) then
        call this_error_out ('ERROR READING PARTICLE COORDINATES')
        return
      endif
    enddo
  endif

  if (rb_com%bunch_charge /= 0) bunch%charge_tot = rb_com%bunch_charge
  sum_charge = sum(p(:)%charge)
  if (bunch%charge_tot == 0) then
    bunch%charge_tot = sum_charge
  elseif (sum_charge == 0) then
    p%charge = bunch%charge_tot / n_particle
  else
    p%charge = p%charge * bunch%charge_tot / sum_charge
  endif
enddo

err = .false.

!---------------------------------------------------------------------------------------------------
contains

!+
! Function remove_first_number (line, ix_word, parens_expected, in_parens) result (pop_ok)
!
! Pop the leading number (which has just be read) off of line leaving the next number at the
! start of the line ready to be read.
!
! Input:
!   line            -- Character(*)
!   ix_word         -- integer: length of existing first word in line
!   parens_expected -- Chracter(*): Tells if parentheses may be expected.
!                       For example: 'x)(' means a ")' and then a '(' may be present after the number.
!
!
! Output
!   line          -- Character(*): Line with first word removed
!   ix_word       -- integer: length of new first word in line
!   in_parens     -- Logical: Inside a parenthesis '(...)' construct?
!   pop_ok        -- Logical: True if no errors encountered.
!-

function remove_first_number (line, ix_word, parens_expected, in_parens) result (pop_ok)

implicit none

integer ix_word
character(*) line, parens_expected
character(4) expect
logical in_parens, pop_ok


!

pop_ok = .false.
expect = parens_expected

! Remove leading '(' if present

if (expect(1:1) == '(') then
  if (line(1:1) == '(') then
    if (in_parens) then
      call this_error_out ('NESTED PARENTHESES FOUND', 'IN LINE: ' // trim(line_in))
      return
    endif
    in_parens = .true.
    call string_trim (line(2:), line, ix_word)
  endif

  expect = expect(2:)  
endif

! Remove word with possible trailing comma. But do not remove parenthesis.

ix = index(line(:ix_word), ',')
if (ix /= 0) ix_word = ix - 1

ix = index(line(:ix_word), '(')
if (ix /= 0) ix_word = ix - 1

ix = index(line(:ix_word), ')')
if (ix /= 0) ix_word = ix - 1

call string_trim (line(ix_word+1:), line, ix_word)
if (line(1:1) == ',') call string_trim(line(2:), line, ix_word)
if (expect(1:1) == 'x') expect = expect(2:)  

! Remove trailing ')' if present

if (expect(1:1) == ')') then
  if (line(1:1) == ')') then
    if (.not. in_parens) then
      call this_error_out ('MISMATCHED PARENTHESES ")" FOUND', 'IN LINE: ' // trim(line_in))
      return
    endif
    in_parens = .false.
    call string_trim (line(2:), line, ix_word)
    if (line(1:1) == ',') call string_trim(line(2:), line, ix_word)
  endif

  expect = expect(2:)  
endif

! Remove trailing '(' if present

if (expect(1:1) == '(') then
  if (line(1:1) == '(') then
    if (in_parens) then
      call this_error_out ('NESTED PARENTHESES FOUND', 'IN LINE: ' // trim(line_in))
      return
    endif
    in_parens = .true.
    call string_trim (line(2:), line, ix_word)
  endif

  expect = expect(2:)  
endif

!

if (line(1:1) == ',') then
  call this_error_out ('MISPLACED COMMA FOUND', 'IN LINE: ' // trim(line_in))
  return
endif

pop_ok = .true.

end function remove_first_number

!---------------------------------------------------------------------------------------------------
! contains

subroutine this_error_out (what, what2)

character(*) what
character(*), optional :: what2

call out_io (s_error$, r_name, 'IN FILE: ' // trim(rb_com%file_name), what, what2)

end subroutine this_error_out

end subroutine tao_read_beam

end module
