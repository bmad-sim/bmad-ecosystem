module beam_file_io

use spin_mod

implicit none

integer, parameter :: ascii$ = 1, binary$ = 2, hdf5$ = 3

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine write_beam_file (file_name, ele, beam, new_file, file_format)
!
! Routine to write a beam file.
!
! Input:
!   file_name     -- character(*): Name of file.
!   ele           -- ele_struct: Element at which beam is at.
!   beam          -- beam_struct: Beam to write
!   new_file      -- logical, optional: New file or append? Default = True.
!   file_format   -- logical, optional: binary$, ascii$, or hdf5$ (default).
!-

subroutine write_beam_file (file_name, ele, beam, new_file, file_format)

type (ele_struct) ele
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p

integer j, iu, ib, ip
integer, optional :: file_format

character(*) file_name
character(*), parameter :: r_name = 'write_beam_file'

logical, optional :: new_file

!

iu = lunget()

if (integer_option(hdf5$, file_format) == binary$) then
  if (logic_option(.true., new_file)) then
    open (iu, file = file_name, form = 'unformatted')
    write (iu) '!BIN::3'
  else
    open (iu, file = file_name, form = 'unformatted', access = 'append')
  endif
elseif (integer_option(hdf5$, file_format) == ascii$) then
  if (logic_option(.true., new_file)) then
    open (iu, file = file_name)
    write (iu, '(a)') '!ASCII::3'
  else
    open (iu, file = file_name, access = 'append')
  endif
else
  call out_io (s_fatal$, r_name, 'HDF5 FORMAT NOT YET SUPPORTED!')
  call err_exit
endif

!

if (integer_option(hdf5$, file_format) == binary$) then
  write (iu) ele%ix_ele, size(beam%bunch), size(beam%bunch(1)%particle)
  do ib = 1, size(beam%bunch)
    bunch => beam%bunch(ib)
    write (iu) bunch%particle(1)%species, bunch%charge_tot, bunch%z_center, bunch%t_center, size(bunch%particle)
    do ip = 1, size(bunch%particle)
      p => bunch%particle(ip)
      write (iu) p%vec, p%charge, p%state, p%spin, p%ix_ele, p%location
    enddo
  enddo
elseif (integer_option(hdf5$, file_format) == ascii$) then
  write (iu, *) ele%ix_ele, '  ! ix_ele' 
  write (iu, *) size(beam%bunch), '  ! n_bunch'
  write (iu, *) size(beam%bunch(1)%particle), '  ! n_particle'
  do ib = 1, size(beam%bunch)
    bunch => beam%bunch(ib)
    write (iu, *) 'BEGIN_BUNCH'
    write (iu, *) '  ', trim(species_name(bunch%particle(1)%species))
    write (iu, *) bunch%charge_tot, '  ! bunch_charge_tot'
    write (iu, *) bunch%z_center,   '  ! z_center'
    write (iu, *) bunch%t_center,   '  ! t_center'
    do ip = 1, size(bunch%particle)
      p => bunch%particle(ip)
      write (iu, '(6es19.10, es14.5, i6, 3es19.10, i6, i3)') &
            p%vec, p%charge, p%state, p%spin, p%ix_ele, p%location
    enddo
    write (iu, *) 'END_BUNCH'
  enddo
endif

close (iu)

end subroutine write_beam_file

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine read_beam_file (file_name, beam, beam_init, set_from_beam_init, err_flag)
!
! Subroutine to read in a beam definition file.
! If set_from_beam_init is True, the following components of beam_init are used to rescale the beam:
!     %n_bunch
!     %n_particle
!     %bunch_charge
! If set_from_beam_init is False, the above components of beam_init are set to the appropriate
! values found in the beam file.
! 
! Modules needed:
!   use beam_file_io
!
! Input:
!   file_name           -- character(*): Name of beam file.
!   set_from_beam_init  -- logical: See above.
!   beam_init           -- beam_init_struct: See above.
!
! Output:
!   beam        -- Beam_struct: Structure holding the beam information.
!   beam_init   -- beam_init_struct: See above.
!   err_flag    -- Logical: Set True if there is an error. False otherwise.
!+ 

subroutine read_beam_file (file_name, beam, beam_init, set_from_beam_init, err_flag)

type (beam_struct), target :: beam
type (beam_init_struct) beam_init
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p(:)
type (coord_struct) orb_init

integer i, j, k, ix, iu, ix_word, ios, ix_ele, species
integer n_bunch, n_particle, n_particle_lines, ix_lost

real(rp) vec(6), sum_charge, bunch_charge
complex(rp) spinor(2)

character(*) file_name
character(300) line, line_in
character(8) file_type
character(16) :: r_name = 'read_beam_file'

logical err_flag, error, in_parens, set_from_beam_init

! Open file and determine whether the file is binary or ascii

err_flag = .true.

iu = lunget()
open (iu, file = file_name, status = 'old', iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN BEAM FILE!')
  return
endif

read (iu, '(a80)') line
if (index(line, '!BINARY') /= 0) then
  call this_error_out ('OLD STYLE BEAM0 FILE NOT SUPPORTED...')
  return
elseif (index(line, '!BIN::2') /= 0) then 
  file_type = 'BIN::2'
elseif (index(line, '!BIN::3') /= 0) then
  file_type = 'BIN::3'
elseif (index(line, '!ASCII::3') /= 0) then
  file_type = 'ASCII::3'
else
  file_type = 'ASCII'
  rewind (iu)
endif

if (file_type(1:3) == 'BIN') then
  close (iu)
  open (iu, file = file_name, form = 'unformatted', status = 'old')
endif

! Read header info

if (file_type(1:5) == 'ASCII') then
  read (iu, *, iostat = ios, err = 9000) ix_ele
  read (iu, *, iostat = ios, err = 9000) n_bunch
  read (iu, *, iostat = ios, err = 9000) n_particle
else
  read (iu) line(1:7)  ! read "!BINARY" line
  read (iu, iostat = ios, err = 9000) ix_ele, n_bunch, n_particle
endif

! Set beam_init

n_particle_lines = n_particle

if (set_from_beam_init) then
  if (beam_init%n_bunch > 0) n_bunch = beam_init%n_bunch
  if (beam_init%n_particle > 0) n_particle = beam_init%n_particle
else
  beam_init%n_bunch = n_bunch
  beam_init%n_particle = n_particle
endif

! Allocate space

call reallocate_beam (beam, n_bunch, n_particle)

! An ascii file, if it is generated by another program, may not include ix_lost or the spin.
! so add the default values

do i = 1, n_bunch

  bunch => beam%bunch(i)
  p => bunch%particle
  p = orb_init   ! init with default params

  if (file_type(1:5) == 'ASCII') then

    read (iu, '(a)', iostat = ios) line
    if (ios /= 0 .or. index(upcase(line), 'BEGIN_BUNCH') == 0) then
      call this_error_out ('NO "BEGIN_BUNCH" MARKER FOUND')
      return
    endif

    read (iu, *, iostat = ios) line
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
    beam_init%species = bunch%particle(1)%species

    read (iu, *, iostat = ios) bunch%charge_tot
    if (ios /= 0) then
      call this_error_out ('ERROR READING BUNCH')
      return
    endif

    
    if (set_from_beam_init) then
      if (beam_init%bunch_charge /= 0) bunch_charge = beam_init%bunch_charge
    else
      bunch_charge = bunch%charge_tot
    endif

    read (iu, *, iostat = ios) bunch%z_center
    if (ios /= 0) then
      call this_error_out ('ERROR READING BUNCH Z_CENTER')
      return
    endif

    read (iu, *, iostat = ios) bunch%t_center
    if (ios /= 0) then
      call this_error_out ('ERROR READING BUNCH T_CENTER')
      return
    endif

    !----------------------------------------
    ! particle coord loop

    j = 0
    do 

      read (iu, '(a)', iostat = ios) line
      if (ios /= 0) then
        call this_error_out ('ERROR READING PARTICLE COORDINATE LINE')
        return
      endif
      line_in = line ! save for error messages

      j = j + 1
      in_parens = .false.

      if (index(upcase(line), 'END_BUNCH') /= 0) exit
      if (j > n_particle) cycle

      p(j)%charge = 0; p(j)%state = alive$; p(j)%spin = 0

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
      if (file_type == 'ASCII::3') then
        read (line, *, iostat = ios) p(j)%spin
        if (.not. remove_first_number (line, ix_word, '', in_parens)) return
        if (.not. remove_first_number (line, ix_word, '', in_parens)) return
        if (.not. remove_first_number (line, ix_word, '', in_parens)) return
      else
        read (line, *, iostat = ios) spinor
        p(j)%spin = spinor_to_vec(spinor)
        if (.not. remove_first_number (line, ix_word, '(x', in_parens)) return
        if (.not. remove_first_number (line, ix_word, 'x)(', in_parens)) return
        if (.not. remove_first_number (line, ix_word, '', in_parens)) return
        if (.not. remove_first_number (line, ix_word, 'x)', in_parens)) return
      endif

      if (ios /= 0) then
        call this_error_out ('ERROR READING PARTICLE SPIN', 'IN LINE: ' // trim(line_in))
        return
      endif


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
    read (iu, iostat = ios) species, bunch%charge_tot, bunch%z_center, bunch%t_center, n_particle_lines
    p%species = species

    if (ios /= 0) then
      call this_error_out ('ERROR READING BUNCH PARAMETERS')
      return
    endif

    do j = 1, n_particle_lines
      if (j > n_particle) exit
      if (file_type == 'BIN::3') then
        read (iu, iostat = ios) p(j)%vec, p(j)%charge, p(j)%state, p(j)%spin, ix_ele, p(j)%location
      else
        read (iu, iostat = ios) p(j)%vec, p(j)%charge, p(j)%state, spinor, ix_ele, p(j)%location
        p(j)%spin = spinor_to_vec(spinor)
      endif
      if (ios /= 0) then
        call this_error_out ('ERROR READING PARTICLE COORDINATES')
        return
      endif
    enddo
  endif

  if (bunch_charge /= 0) bunch%charge_tot = bunch_charge
  sum_charge = sum(p(:)%charge)
  if (bunch%charge_tot == 0) then
    bunch%charge_tot = sum_charge
  elseif (sum_charge == 0) then
    p%charge = bunch%charge_tot / n_particle
  else
    p%charge = p%charge * bunch%charge_tot / sum_charge
  endif
enddo

close (iu)
err_flag = .false.
return

!

9000 continue
if (ios > 0) then
  call out_io (s_error$, r_name, 'ERROR READING BEAM HEADER INFO IN FILE: ' // trim(file_name))
  close (iu)
  return
endif

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

call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), what, what2)

end subroutine this_error_out

end subroutine read_beam_file

end module
