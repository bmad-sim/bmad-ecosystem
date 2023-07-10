module beam_file_io

use equal_mod

implicit none

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine write_beam_file (file_name, beam, new_file, file_format, lat)
!
! Routine to write a beam file.
!
! A '.h5' suffix will be appended to the created file if hdf5$ format is used and file_name does not
! already have a '.h5' or '.hdf5' suffix.
!
! Input:
!   file_name     -- character(*): Name of file.
!   beam          -- beam_struct: Beam to write
!   new_file      -- logical, optional: New file or append? Default = True.
!   file_format   -- logical, optional: ascii$, or hdf5$ (default).
!   lat           -- lat_struct, optional: If present, lattice info will be writen to hdf5 files.
!
! Output:
!   file_name
!-

subroutine write_beam_file (file_name, beam, new_file, file_format, lat)

type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p
type (lat_struct), optional :: lat

integer j, iu, ib, ip, ix_ele, n, n0
integer, optional :: file_format

character(*) file_name
character(200) full_name
character(*), parameter :: r_name = 'write_beam_file'

logical, optional :: new_file
logical error, append

!

iu = lunget()
call fullfilename (file_name, full_name)

if (integer_option(hdf5$, file_format) == hdf5$) then
  n = len_trim(full_name)
  if (full_name(n-2:n) /= '.h5' .and. full_name(n-4:n) /= '.hdf5') full_name = trim(full_name) // '.h5'

  append = .not. logic_option(.true., new_file)
  call hdf5_write_beam(full_name, beam%bunch, append, error, lat)
  return
endif

!

if (logic_option(.true., new_file)) then
  open (iu, file = full_name)
  write (iu, '(a)') '!ASCII::3'
else
  open (iu, file = full_name, access = 'append')
endif

write (iu, *) beam%bunch(1)%particle(1)%ix_ele, '  ! ix_ele' 
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

close (iu)

end subroutine write_beam_file

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine read_beam_file (file_name, beam, beam_init, err_flag, ele, print_p0c_shift_warning, conserve_momentum)
!
! Subroutine to read in a beam definition file.
! If non_zero, the following components of beam_init are used to rescale the beam:
!     %n_bunch
!     %n_particle
!     %bunch_charge
!
! If the beam file has '.h5' or '.hdf5' suffix then the file is taken to be an HDF5 file.
! Otherwise the file is assumed to be ASCII.
!
! Input:
!   file_name   -- character(*): Name of beam file.
!   beam_init   -- beam_init_struct: See above.
!   ele         -- ele_struct, optional: Element with reference energy, etc.
!   print_p0c_shift_warning   -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!   shift_momentum            -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!
! Output:
!   beam        -- Beam_struct: Structure holding the beam information.
!   err_flag    -- Logical: Set True if there is an error. False otherwise.
!+ 

subroutine read_beam_file (file_name, beam, beam_init, err_flag, ele, print_p0c_shift_warning, conserve_momentum)

type (beam_struct), target :: beam
type (beam_init_struct) beam_init
type (ele_struct), optional :: ele
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p(:)
type (coord_struct), allocatable :: p_temp(:)
type (coord_struct) orb_init
type (pmd_header_struct) :: pmd_header

integer i, j, k, n, np, ix, iu, ix_word, ios, ix_ele, species
integer n_bunch, n_particle, n_particle_lines, ix_lost

real(rp) vec(6), sum_charge, bunch_charge
complex(rp) spinor(2)

character(*) file_name
character(200) full_name
character(300) line, line_in
character(8) file_type
character(*), parameter :: r_name = 'read_beam_file'

logical err_flag, error, in_parens, valid
logical, optional :: print_p0c_shift_warning, conserve_momentum

!

err_flag = .true.

call fullfilename(file_name, full_name, valid)
if (.not. valid) then
  call out_io (s_error$, r_name, 'NOT A VALID FILE NAME:' // file_name)
  return
endif

! HDF5 file

n = len_trim(full_name)
if (full_name(max(1,n-4):n) == '.hdf5' .or. full_name(max(1,n-2):n) == '.h5') then
  call hdf5_read_beam (full_name, beam, err_flag, ele, pmd_header, print_p0c_shift_warning, conserve_momentum)
  if (err_flag) return

  np = beam_init%n_particle
  if (np > 0) then
    do i = 1, size(beam%bunch)
      p => beam%bunch(i)%particle
      if (size(p) < beam_init%n_particle) then
        call out_io (s_warn$, r_name, &
                'Number of particles ' // int_str(size(p)) // ' defined in beam file: ' // full_name, &
                'less than number of particles wanted which is set by beam_init%n_particle: ' // int_str(np), &
                'The setting of beam_init%n_particle will be ignored.')
      endif

      np = min(size(p), np)
      call move_alloc (beam%bunch(i)%particle, p_temp)
      allocate (beam%bunch(i)%particle(np))
      beam%bunch(i)%particle = p_temp(1:np)
      deallocate (p_temp)
      call re_allocate (beam%bunch(i)%ix_z, np)
    enddo
  endif
  return
endif

! Open file and determine whether the file is binary or ascii

iu = lunget()
open (iu, file = full_name, status = 'old', iostat = ios, action = 'read')
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN BEAM FILE: ' // quote(file_name))
  return
endif

read (iu, '(a80)') line
if (index(line, '!BINARY') /= 0) then
  call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'OLD STYLE BEAM0 FILE NOT SUPPORTED...')
  return
elseif (index(line, '!BIN::2') /= 0) then 
  file_type = 'BIN::2'
elseif (index(line, '!BIN::3') /= 0) then
  file_type = 'BIN::3'
elseif (index(line, '!ASCII::3') /= 0) then
  file_type = 'ASCII::3'
else
  file_type = 'ASCII::3'
  rewind (iu)
endif

if (file_type(1:3) == 'BIN') then
  close (iu)
  open (iu, file = full_name, form = 'unformatted', status = 'old', action = 'read')
endif

! Read header info

if (file_type == 'ASCII::3') then
  read (iu, *, iostat = ios, err = 9000) ix_ele
  read (iu, *, iostat = ios, err = 9000) n_bunch
  read (iu, *, iostat = ios, err = 9000) n_particle
else
  read (iu) line(1:7)  ! read "!BINARY" line
  read (iu, iostat = ios, err = 9000) ix_ele, n_bunch, n_particle
endif

! Set beam_init

n_particle_lines = n_particle

if (beam_init%n_bunch > 0) n_bunch = beam_init%n_bunch
if (beam_init%n_particle > 0) n_particle = beam_init%n_particle

! Allocate space

call reallocate_beam (beam, n_bunch, n_particle)

! An ascii file, if it is generated by another program, may not include ix_lost or the spin.
! so add the default values

do i = 1, n_bunch
  bunch => beam%bunch(i)
  p => bunch%particle
  p = orb_init   ! init with default params

  if (file_type == 'ASCII::3') then

    read (iu, '(a)', iostat = ios) line
    if (ios /= 0 .or. index(upcase(line), 'BEGIN_BUNCH') == 0) then
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'NO "BEGIN_BUNCH" MARKER FOUND')
      return
    endif

    read (iu, *, iostat = ios) line
    if (is_real(line, .true.)) then
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'OLD STYLE FORMAT DOES NOT INCLUDE BUNCH SPECIES', 'PLEASE CORRECT.')
      return
    endif

    ix = species_id(line)
    if (ix == invalid$) then
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'BAD SPECIES NAME: ' // trim(line))
      return
    endif
    bunch%particle%species = ix
    beam_init%species = species_name(bunch%particle(1)%species)

    read (iu, *, iostat = ios) bunch%charge_tot
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING BUNCH')
      return
    endif
    if (beam_init%bunch_charge /= 0) bunch%charge_tot = beam_init%bunch_charge

    read (iu, *, iostat = ios) bunch%z_center
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING BUNCH Z_CENTER')
      return
    endif

    read (iu, *, iostat = ios) bunch%t_center
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING BUNCH T_CENTER')
      return
    endif

    !----------------------------------------
    ! particle coord loop

    j = 0
    do
      read (iu, '(a)', iostat = ios) line
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING PARTICLE COORDINATE LINE')
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
          call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING PARTICLE COORDINATES', 'IN LINE: ' // trim(line_in))
          return
        endif
        if (.not. remove_first_number (line, ix_word, '', in_parens)) return
      enddo

      if (ix_word == 0) goto 8000
      read (line, *, iostat = ios) p(j)%charge
      if (ios /= 0 .or. ix_word == 0) then
        call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING PARTICLE CHARGE', 'IN LINE: ' // trim(line_in))
        return
      endif
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return

      if (ix_word == 0) goto 8000
      read (line, *, iostat = ios) p(j)%state
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING PARTICLE "STATE"', 'IN LINE: ' // trim(line_in))
        return
      endif
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return

      if (ix_word == 0) goto 8000
      read (line, *, iostat = ios) p(j)%spin
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING PARTICLE SPIN', 'IN LINE: ' // trim(line_in))
        return
      endif

      if (.not. remove_first_number (line, ix_word, '', in_parens)) return
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return

      if (ix_word == 0) goto 8000
      read (line, *, iostat = ios) p(j)%ix_ele
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING ELEMENT INDEX', 'IN LINE: ' // trim(line_in))
        return
      endif
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return

      if (ix_word == 0) goto 8000
      read (line, *, iostat = ios) p(j)%location
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING PARTICLE LOCATION', 'IN LINE: ' // trim(line_in))
        return
      endif
      if (.not. remove_first_number (line, ix_word, '', in_parens)) return

      8000 continue
      if (in_parens .or. ix_word /= 0) then
        call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'UNMATCHED PARENTHESIS IN LINE: ' // trim(line_in))
        return
      endif
    enddo

  !------------------------------------------------------------------------------------
  ! Binary file

  else
    read (iu, iostat = ios) species, bunch%charge_tot, bunch%z_center, bunch%t_center, n_particle_lines
    p%species = species

    if (ios /= 0) then
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING BUNCH PARAMETERS')
      return
    endif

    do j = 1, n_particle_lines
      if (file_type == 'BIN::3') then
        read (iu, iostat = ios) p(j)%vec, p(j)%charge, p(j)%state, p(j)%spin, ix_ele, p(j)%location
      else
        read (iu, iostat = ios) p(j)%vec, p(j)%charge, p(j)%state, spinor, ix_ele, p(j)%location
        p(j)%spin = spinor_to_vec(spinor)
      endif
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING PARTICLE COORDINATES')
        return
      endif
      if (j == n_particle) exit
    enddo
  endif

  if (j < n_particle) then
    call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), &
            'NUMBER OF PARTICLES DEFINED IN THE FILE \i0\ IS LESS THAN THE NUMBER OF DESIRED PARTICLES \i0\.', &
            i_array = [j, n_particle])
    return
  endif

  sum_charge = sum(p(:)%charge)
  if (bunch%charge_tot == 0) then
    bunch%charge_tot = sum_charge
  elseif (sum_charge == 0) then
    p%charge = bunch%charge_tot / n_particle
  else
    p%charge = p%charge * bunch%charge_tot / sum_charge
  endif

  bunch%charge_live = sum(p%charge, (p%state == alive$))
  bunch%n_live = count(p%state == alive$)

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
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'NESTED PARENTHESES FOUND', 'IN LINE: ' // trim(line_in))
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
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'MISMATCHED PARENTHESES ")" FOUND', 'IN LINE: ' // trim(line_in))
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
      call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'NESTED PARENTHESES FOUND', 'IN LINE: ' // trim(line_in))
      return
    endif
    in_parens = .true.
    call string_trim (line(2:), line, ix_word)
  endif

  expect = expect(2:)  
endif

!

if (line(1:1) == ',') then
  call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'MISPLACED COMMA FOUND', 'IN LINE: ' // trim(line_in))
  return
endif

pop_ok = .true.

end function remove_first_number

end subroutine read_beam_file

end module
