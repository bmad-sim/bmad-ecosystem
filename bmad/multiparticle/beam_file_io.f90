module beam_file_io

use bmad_routine_interface

implicit none

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine write_beam_file (file_name, beam, new_file, file_format, lat, alive_only)
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
!   file_format   -- logical, optional: ascii$, or hdf5$ (default). old_ascii$ (deprecated) is still accepted.
!   lat           -- lat_struct, optional: If present, lattice info will be writen to hdf5 files.
!   alive_only    -- logical, optional: Only write live (includes pre_born) particles to the file? Default is False.
!-

subroutine write_beam_file (file_name, beam, new_file, file_format, lat, alive_only)

type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p
type (lat_struct), optional :: lat

integer j, iu, ib, ip, ix_ele, n, n0, format
integer, optional :: file_format

character(*) file_name
character(200) full_name
character(*), parameter :: r_name = 'write_beam_file'

logical, optional :: new_file, alive_only
logical error, append

!

call fullfilename (file_name, full_name)
n = len_trim(full_name)

if (present(file_format)) then
  format = file_format
elseif ((index(full_name, '.hdf5') == n-4 .and. n > 4) .or. (index(full_name, '.h5') == n-2 .and. n > 2)) then
  format = hdf5$
else
  format = ascii$
endif

if (format == ascii$) then
  call write_ascii_beam_file(full_name, beam, new_file, alive_only)
  return
endif

!

iu = lunget()

if (integer_option(hdf5$, format) == hdf5$) then
  n = len_trim(full_name)
  if (full_name(n-2:n) /= '.h5' .and. full_name(n-4:n) /= '.hdf5') full_name = trim(full_name) // '.h5'

  append = .not. logic_option(.true., new_file)
  call hdf5_write_beam(full_name, beam%bunch, append, error, lat, alive_only)
  return
endif

! Old_ASCII format ------------------------------

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
  write (iu, *) bunch%charge_tot, '  ! charge_tot'
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
! Subroutine write_ascii_beam_file (file_name, beam, new_file, alive_only)
!
! Routine to write a beam file in ASCII format (version 4).
!
! Input:
!   file_name     -- character(*): Name of file.
!   beam          -- beam_struct: Beam to write
!   new_file      -- logical, optional: New file or append? Default = True.
!   alive_only    -- logical, optional: Only write live (includes pre_born) particles to the file? Default is False.
!-

subroutine write_ascii_beam_file (file_name, beam, new_file, alive_only)

type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p

integer j, iu, ib, ic, ip, ix, ix_ele, n, n0, n_col, n_alive

character(*) file_name
character(20) col
character(200) name, col_out
character(600) line
character(12) colvec(20), cfmt(20)
character(*), parameter :: r_name = 'write_ascii_beam_file'

logical, optional :: new_file, alive_only
logical error, append, dt_ref0
logical spin0, field0, phase0, ix_branch0, loc0, dir0, tdir0, charge0, species0, s0, t0, p0c0

!

iu = lunget()
call fullfilename (file_name, name)

if (logic_option(.true., new_file)) then
  open (iu, file = name, recl = 600)
else
  open (iu, file = name, access = 'append', recl = 600)
endif

colvec = ''
colvec(1:17) = [character(12):: 'index', 'vec', 'p0c', 's_position', 'time', 'charge', 'spin', 'field', &
                 'phase', 'state', 'ix_ele', 'ix_branch', 'location', 'species', 'direction', 'time_dir', 'dt_ref']
spin0 = .true.; field0 = .true.; phase0 = .true.; ix_branch0 = .true.; loc0 = .true.; dir0 = .true.
tdir0 = .true.; charge0 = .true.; species0 = .true.; s0 = .true.; t0 = .true.; p0c0 = .true.; dt_ref0 = .true.
do ib = 1, size(beam%bunch)
  bunch => beam%bunch(ib)
  p => bunch%particle(1)
  if (logic_option(alive_only, .false.) .and. p%state /= alive$ .and. p%state /= pre_born$) cycle
  if (any(bunch%particle%spin(1) /= 0))              spin0 = .false.
  if (any(bunch%particle%spin(2) /= 0))              spin0 = .false.
  if (any(bunch%particle%spin(3) /= 0))              spin0 = .false.
  if (any(bunch%particle%field(1) /= 0))             field0 = .false.
  if (any(bunch%particle%field(2) /= 0))             field0 = .false.
  if (any(bunch%particle%phase(1) /= 0))             phase0 = .false.
  if (any(bunch%particle%phase(2) /= 0))             phase0 = .false.
  if (any(bunch%particle%ix_branch /= p%ix_branch))  ix_branch0 = .false.
  if (any(bunch%particle%location /= p%location))    loc0 = .false.
  if (any(bunch%particle%direction /= p%direction))  dir0 = .false.
  if (any(bunch%particle%time_dir /= p%time_dir))    tdir0 = .false.
  if (any(bunch%particle%charge /= p%charge))        charge0 = .false.
  if (any(bunch%particle%species /= p%species))      species0 = .false.
  if (any(bunch%particle%s /= p%s))                  s0 = .false.
  if (any(bunch%particle%t /= p%t))                  t0 = .false.
  if (any(bunch%particle%p0c /= p%p0c))              p0c0 = .false.
  if (any(bunch%particle%dt_ref /= 0))               dt_ref0 = .false.
enddo

!!! if (spin0)      call remove_col('spin', colvec)  ! Too confusing to not write out the spin even if the spin is zero.
if (field0)     call remove_col('field', colvec)
if (phase0)     call remove_col('phase', colvec)
if (dir0)       call remove_col('direction', colvec)
if (tdir0)      call remove_col('time_dir', colvec)
if (charge0)    call remove_col('charge', colvec)
if (species0)   call remove_col('species', colvec)
if (ix_branch0) call remove_col('ix_branch', colvec)
if (dt_ref0)    call remove_col('dt_ref', colvec)

!

do ib = 1, size(beam%bunch)
  bunch => beam%bunch(ib)

  ! Write header

  p => bunch%particle(1)
  n_alive = count(bunch%particle%state == alive$) + count(bunch%particle%state == pre_born$)
  write (iu, '(a, i6)') '# ix_bunch       =', ib
  write (iu, '(a, i6)') '# n_particle     =', size(bunch%particle)
  write (iu, '(a, i6)') '# n_alive        =', n_alive
  write (iu, '(a, i6)') '# n_dead         =', size(bunch%particle) - n_alive

  call headwrite_re(charge0, 'charge', colvec, p%charge)

  call headwrite_int(ix_branch0, 'ix_branch', colvec, p%ix_branch)
  call headwrite_int(dir0, 'direction', colvec, p%direction)
  call headwrite_int(tdir0, 'time_dir', colvec, p%time_dir)

  call headwrite_str(species0, 'species', colvec, species_name(p%species))
  call headwrite_str(loc0, 'location', colvec, location_name(p%location))

  ! Write table header

  line = ''

  do ic = 1, size(colvec)
    col = colvec(ic)
    if (col == '') cycle
    select case (col)
    case ('ix_branch', 'direction', 'time_dir');   write (line, '(a, a10)') trim(line), trim(col)
    case ('s_position', 'time', 'charge', 'p0c');  write (line, '(a, a23)') trim(line), trim(col)
    case ('index', 'ix_ele'); write (line, '(a, a8)') trim(line), trim(col)
    case ('vec');             write (line, '(a, 6a23)') trim(line), 'x', 'px', 'y', 'py', 'z', 'pz'
    case ('spin');            write (line, '(a, 6a23)') trim(line), 'spin_x', 'spin_y', 'spin_z'
    case ('field');           write (line, '(a, 6a23)') trim(line), 'field_x', 'field_y'
    case ('phase');           write (line, '(a, 6a23)') trim(line), 'phase_x', 'phase_y'
    case ('state');           write (line, '(a, a14)') trim(line), 'state'
    case ('species');         write (line, '(a, a14)') trim(line), 'species'
    case ('ix_user');         write (line, '(a, a16)') trim(line), 'ix_user' 
    case ('location');        write (line, '(a, a16)') trim(line), 'location'
    case ('dt_ref');          write (line, '(a, a23)') trim(line), 'dt_ref'
    case ('beta');            write (line, '(a, a23)') trim(line), 'beta'
    case ('r');               write (line, '(a, a23)') trim(line), 'r'
    case ('e_potential');     write (line, '(a, a23)') trim(line), 'e_potential'
    case default
      call err_exit
    end select
  enddo

  line(1:2) = '#!'
  write (iu, '(a)') trim(line)

  ! Write body.

  do ip = 1, size(bunch%particle)
    p => bunch%particle(ip)
    if (logic_option(alive_only, .false.) .and. p%state /= alive$ .and. p%state /= pre_born$) cycle
    line = ''
    do ic = 1, size(colvec)
      col = colvec(ic)
      if (col == '') cycle
      select case (col)
      case ('index');       write (line, '(a, i8)') trim(line), ip
      case ('ix_ele');      write (line, '(a, i8)') trim(line), p%ix_ele
      case ('ix_branch');   write (line, '(a, i10)') trim(line), p%ix_branch
      case ('direction');   write (line, '(a, i10)') trim(line), p%direction
      case ('time_dir');    write (line, '(a, i10)') trim(line), p%time_dir
      case ('s_position');  write (line, '(a,  es23.15)') trim(line), p%s
      case ('time');        write (line, '(a,  es23.15)') trim(line), p%t
      case ('charge');      write (line, '(a,  es23.15)') trim(line), p%charge
      case ('p0c');         write (line, '(a,  es23.15)') trim(line), p%p0c
      case ('dt_ref');      write (line, '(a,  es23.15)') trim(line), p%dt_ref
      case ('e_potential'); write (line, '(a,  es23.15)') trim(line), p%E_potential
      case ('beta');        write (line, '(a,  es23.15)') trim(line), p%beta  
      case ('r');           write (line, '(a,  es23.15)') trim(line), p%r
      case ('vec');         write (line, '(a, 6es23.15)') trim(line), p%vec
      case ('spin');        write (line, '(a, 3es23.15)') trim(line), p%spin
      case ('field');       write (line, '(a, 2es23.15)') trim(line), p%field
      case ('phase');       write (line, '(a, 2es23.15)') trim(line), p%phase
      case ('state');       write (line, '(a, a14)') trim(line), trim(state_name(p%state))
      case ('species');     write (line, '(a, a14)') trim(line), trim(species_name(p%species))
      case ('location');    write (line, '(a, a16)') trim(line), trim(location_name(p%location))
      case ('ix_user');     write (line, '(a, i16)') trim(line), p%ix_user
      end select
    enddo
    write (iu, '(a)') trim(line)
  enddo

enddo

close (iu)

!-----------------------------------------------------------------------------
contains

subroutine remove_col(who, colvec)

integer ix, n
character(*) who
character(12) colvec(:)

n = size(colvec)
ix = find_location(colvec, who)
colvec(ix:n-1) = colvec(ix+1:n)

end subroutine remove_col

!-----------------------------------------------------------------------------
! contains

subroutine headwrite_re(is_const, who, colvec, re_val)

character(*) who
character(12) colvec(:)
real(rp) re_val
logical is_const

!

if (any(colvec == who) .or. .not. is_const) return
write (iu, '(2a, t18, a, es24.16)') '# ', who, '=', re_val

end subroutine headwrite_re

!-----------------------------------------------------------------------------
! contains

subroutine headwrite_int(is_const, who, colvec, int_val)

character(*) who
character(12) colvec(:)
integer int_val
logical is_const

!

if (any(colvec == who) .or. .not. is_const) return
write (iu, '(2a, t18, a, i0)') '# ', who, '= ', int_val

end subroutine headwrite_int

!-----------------------------------------------------------------------------
! contains

subroutine headwrite_str(is_const, who, colvec, name)

character(*) who, name
character(12) colvec(:)
logical is_const

!

if (any(colvec == who) .or. .not. is_const) return
write (iu, '(2a, t18, 2a)') '# ', who, '= ', name

end subroutine headwrite_str

end subroutine write_ascii_beam_file

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine read_beam_file (file_name, beam, beam_init, err_flag, ele, print_mom_shift_warning, conserve_momentum)
!
! Subroutine to read in a beam definition file.
! If non_zero, the following components of beam_init are used to rescale the beam:
!     %n_bunch
!     %n_particle
!     %bunch_charge -> charge_tot
!     %species
!
! If the beam file has '.h5' or '.hdf5' suffix then the file is taken to be an HDF5 file.
! Otherwise the file is assumed to be ASCII.
!
! Input:
!   file_name   -- character(*): Name of beam file.
!   beam_init   -- beam_init_struct: See above.
!   ele         -- ele_struct, optional: Element with reference energy, etc.
!   print_mom_shift_warning   -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!   shift_momentum            -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!
! Output:
!   beam        -- Beam_struct: Structure holding the beam information.
!   err_flag    -- Logical: Set True if there is an error. False otherwise.
!+ 

subroutine read_beam_file (file_name, beam, beam_init, err_flag, ele, print_mom_shift_warning, conserve_momentum)

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

real(rp) vec(6), sum_charge, charge_tot
complex(rp) spinor(2)

character(*) file_name
character(200) full_name
character(300) line, line_in
character(8) file_type
character(*), parameter :: r_name = 'read_beam_file'

logical err_flag, error, in_parens, valid
logical, optional :: print_mom_shift_warning, conserve_momentum

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
  call hdf5_read_beam (full_name, beam, err_flag, ele, pmd_header, print_mom_shift_warning, conserve_momentum)
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

! Note: BIN type files no longer generated by Bmad (HDF5 is always used instead).

read (iu, '(a80)') line
if (index(line, '!BINARY') /= 0 .or. index(line, '!BIN::2') /= 0) then
  call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'OLD STYLE BINARY BEAM FILE NOT SUPPORTED...')
  return
elseif (index(line, '!BIN::3') /= 0) then
  file_type = 'BIN::3'
elseif (index(line, '!ASCII::3') /= 0) then
  file_type = 'ASCII::3'
else
  do
    if (len_trim(line) /= 0) exit
    read (iu, '(a)') line
  enddo

  if (line(1:1) == '#') then
    close (iu)
    call read_beam_ascii(file_name, beam, beam_init, err_flag)
    return
  endif

  file_type = 'ASCII::3'
  backspace (iu)
endif

! Old_ASCII format -------------------------------------------------------------------

if (file_type == 'BIN::3') then
  close (iu)
  open (iu, file = full_name, form = 'unformatted', status = 'old', action = 'read')
endif

! Read header info

if (file_type == 'ASCII::3') then
  read (iu, *, iostat = ios, err = 9000) ix_ele
  read (iu, *, iostat = ios, err = 9000) n_bunch
  read (iu, *, iostat = ios, err = 9000) n_particle

else
  read (iu) line(1:7)  ! Read "!BIN::" line
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
      read (iu, iostat = ios) p(j)%vec, p(j)%charge, p(j)%state, p(j)%spin, ix_ele, p(j)%location
      if (ios /= 0) then
        call out_io (s_error$, r_name, 'IN FILE: ' // trim(file_name), 'ERROR READING PARTICLE COORDINATES')
        return
      endif
      if (j == n_particle) exit
    enddo
  endif

  !-------------------------

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

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine read_beam_ascii (file_name, beam, beam_init, err_flag, ele, print_mom_shift_warning, conserve_momentum)
!
! Subroutine to read in a beam definition file.
! If non_zero, the following components of beam_init are used to rescale the beam:
!     %n_bunch
!     %n_particle
!     %charge_tot
!
! If the beam file has '.h5' or '.hdf5' suffix then the file is taken to be an HDF5 file.
! Otherwise the file is assumed to be ASCII.
!
! Input:
!   iu          -- integer: File unit number
!   file_name   -- character(*): Name of beam file.
!   beam_init   -- beam_init_struct: See above.
!   ele         -- ele_struct, optional: Element with reference energy, etc.
!   print_mom_shift_warning   -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!   shift_momentum            -- logical, optional: Default is True. See hdf5_read_beam doc. Only used when reading hdf5 file.
!
! Output:
!   beam        -- Beam_struct: Structure holding the beam information.
!   err_flag    -- Logical: Set True if there is an error. False otherwise.
!+ 

subroutine read_beam_ascii (file_name, beam, beam_init, err_flag)

type (beam_struct), target :: beam
type (beam_init_struct) beam_init
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p
type (coord_struct) p0

real(rp) charge_tot, tim

integer iu, ic, ip, np, ix, ios, n_particle, n_bunch, n_col, n_line

character(*) file_name
character(*), parameter :: r_name = 'read_beam_ascii'
character(20) col(30), acol
character(200) full_name
character(600) line, str, line_saved

logical err_flag, err, valid, beta_found

!

n_bunch = 0
err_flag = .true.
err = .false.
n_line = 0

call fullfilename(file_name, full_name, valid)
if (.not. valid) then
  call out_io (s_error$, r_name, 'NOT A VALID FILE NAME:' // file_name)
  return
endif

iu = lunget()
open (iu, file = full_name, status = 'old', iostat = ios, action = 'read')
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN BEAM FILE: ' // quote(file_name))
  return
endif

! bunch loop
do
  n_bunch = n_bunch + 1
  call reallocate_beam(beam, n_bunch, 1000, extend = .true.)
  bunch => beam%bunch(n_bunch)

  p0 = coord_struct()
  beta_found = .false.

  ! Read bunch header
  header_loop: do
    read (iu, '(a)', iostat = ios) line
    n_line = n_line + 1
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'CANNOT READ BUNCH HEADER IN FILE: ' // file_name, &
                                     'LINE NUMBER ' // int_str(n_line))
      return
    endif

    if (line == '') cycle

    if (line(1:2) == '#!') then
      call string_trim(line(3:), line, ix)
      n_col = 0
      do
        if (ix == 0) exit header_loop
        n_col = n_col + 1
        col(n_col) = line(:ix)
        call string_trim(line(ix+1:), line, ix)
      enddo
    endif

    if (line(1:1) /= '#') then
      call out_io (s_error$, r_name, 'FIRST CHARACTER IN HEADER LINE NOT A "#" CHARACTER: ' // quote(line(1:20)) // '...', &
                                     'IN FILE: ' // file_name, & 
                                     'LINE NUMBER ' // int_str(n_line))
      return
    endif

    if (index(line, '=') == 0) cycle
    call string_trim(line(2:), line, ix)

    select case (downcase(line(:ix)))
    case ('charge_tot');  bunch%charge_tot = read_param(line, n_line, err); if (err) return

    case ('location')
      str = read_string(line, n_line, err); if (err) return
      call read_switch(line(:ix), p0%location, str, n_line, err); if (err) return

    case ('state')
      str = unquote(read_string(line, n_line, err)); if (err) return
      call read_switch(line(:ix), p0%state, str, n_line, err); if (err) return

    case ('species')
      str = unquote(read_string(line, n_line, err)); if (err) return
      call read_switch(line(:ix), p0%species, str, n_line, err); if (err) return

    case ('r');                p0%r           = read_param(line, n_line, err); if (err) return
    case ('s', 's_position');  p0%s           = read_param(line, n_line, err); if (err) return
    case ('t', 'time');        p0%t           = read_param(line, n_line, err); if (err) return
    case ('p0c');              p0%p0c         = read_param(line, n_line, err); if (err) return
    case ('charge');           p0%charge      = read_param(line, n_line, err); if (err) return
    case ('dt_ref');           p0%dt_ref      = read_param(line, n_line, err); if (err) return
    case ('e_potential');      p0%E_potential = read_param(line, n_line, err); if (err) return
    case ('beta');             p0%beta        = read_param(line, n_line, err); if (err) return; beta_found = .true.
    case ('spin');             call read_params(line, n_line, p0%spin, err); if (err) return
    case ('field');            call read_params(line, n_line, p0%field, err); if (err) return
    case ('phase');            call read_params(line, n_line, p0%phase, err); if (err) return
    case ('time_dir');         p0%time_dir    = nint(read_param(line, n_line, err)); if (err) return
    case ('direction');        p0%direction   = nint(read_param(line, n_line, err)); if (err) return
    case ('ix_ele');           p0%ix_ele      = nint(read_param(line, n_line, err)); if (err) return
    case ('ix_branch');        p0%ix_branch   = nint(read_param(line, n_line, err)); if (err) return
    case ('ix_user');          p0%ix_user     = nint(read_param(line, n_line, err)); if (err) return
    end select
  enddo header_loop

  ! Read particle info

  ip = 0
  do 
    read (iu, '(a)', iostat = ios, end = 8000) line
    n_line = n_line + 1
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'CANNOT READ BEAM FILE TABLE IN FILE: ' // file_name, 'AT LINE NUMBER: ' // int_str(n_line))
      return
    endif

    if (line == '') cycle

    if (line(1:1) == '#') then
      call bunch_finalizer(bunch, ip, beam_init, beta_found)
      exit   ! Next bunch
    endif

    ip = ip + 1
    if (ip > size(bunch%particle)) call reallocate_bunch (bunch, 2*ip, .true.)

    p => bunch%particle(ip)
    p = p0

    call string_trim(line, line, ix)
    do ic = 1, n_col
      acol = col(ic)
      select case (downcase(acol))
      case ('x');                call read_component(acol, p%vec(1), line, n_line, ix, err); if (err) return
      case ('px');               call read_component(acol, p%vec(2), line, n_line, ix, err); if (err) return
      case ('y');                call read_component(acol, p%vec(3), line, n_line, ix, err); if (err) return
      case ('py');               call read_component(acol, p%vec(4), line, n_line, ix, err); if (err) return
      case ('z');                call read_component(acol, p%vec(5), line, n_line, ix, err); if (err) return
      case ('pz');               call read_component(acol, p%vec(6), line, n_line, ix, err); if (err) return
      case ('spin_x');           call read_component(acol, p%spin(1), line, n_line, ix, err); if (err) return
      case ('spin_y');           call read_component(acol, p%spin(2), line, n_line, ix, err); if (err) return
      case ('spin_z');           call read_component(acol, p%spin(3), line, n_line, ix, err); if (err) return
      case ('field_x');          call read_component(acol, p%field(1), line, n_line, ix, err); if (err) return
      case ('field_y');          call read_component(acol, p%field(2), line, n_line, ix, err); if (err) return
      case ('phase_x');          call read_component(acol, p%phase(1), line, n_line, ix, err); if (err) return
      case ('phase_y');          call read_component(acol, p%phase(2), line, n_line, ix, err); if (err) return
      case ('s', 's_position');  call read_component(acol, p%s, line, n_line, ix, err); if (err) return
      case ('t', 'time');        call read_component(acol, tim, line, n_line, ix, err); if (err) return; p%t = tim
      case ('charge');           call read_component(acol, p%charge, line, n_line, ix, err); if (err) return
      case ('dt_ref');           call read_component(acol, p%dt_ref, line, n_line, ix, err); if (err) return
      case ('r');                call read_component(acol, p%r, line, n_line, ix, err); if (err) return
      case ('p0c');              call read_component(acol, p%p0c, line, n_line, ix, err); if (err) return
      case ('E_potential');      call read_component(acol, p%E_potential, line, n_line, ix, err); if (err) return
      case ('beta');             call read_component(acol, p%beta, line, n_line, ix, err); if (err) return; beta_found = .true.
      case ('ix_ele');           call read_component_int(acol, p%ix_ele, line, n_line, ix, err); if (err) return
      case ('ix_branch');        call read_component_int(acol, p%ix_branch, line, n_line, ix, err); if (err) return
      case ('ix_user');          call read_component_int(acol, p%ix_user, line, n_line, ix, err); if (err) return
      case ('direction');        call read_component_int(acol, p%direction, line, n_line, ix, err); if (err) return
      case ('time_dir');         call read_component_int(acol, p%time_dir, line, n_line, ix, err); if (err) return
      case ('state');            call read_switch(acol, p%state, line, n_line, err, ix); if (err) return
      case ('species');          call read_switch(acol, p%species, line, n_line, err, ix); if (err) return
      case ('location');         call read_switch(acol, p%location, line, n_line, err, ix); if (err) return
      case ('index')
        ! Value is not used but check that it is an integer
        if (.not. is_integer(line(:ix))) then
          call out_io (s_error$, r_name, 'INDEX COLUMN ENTRY IS NOT AN INTEGER: ' // quote(line(:ix)), 'IN FILE: ' // file_name, &
                                         'AT LINE NUMBER: ' // int_str(n_line))
          return
        endif
        call string_trim(line(ix+1:), line, ix)     ! Ignore value
      case default
        err = .true.
        call out_io(s_error$, r_name, 'COLUMN NAME NOT RECOGNIZED: ' // col(ic), &
                                      'IN FILE: ' // file_name, &
                                      'AT LINE NUMBER: ' // int_str(n_line))
        return
      end select
    enddo
  enddo

  call bunch_finalizer(bunch, ip, beam_init, beta_found)
enddo

!

8000 continue
call bunch_finalizer(bunch, ip, beam_init, beta_found)
err_flag = .false.

!---------------------------------------------------------------------------------------------------
contains

subroutine bunch_finalizer (bunch, n_part, beam_init, beta_found)

type (bunch_struct), target :: bunch
type (beam_init_struct) beam_init
type (coord_struct), pointer :: p(:), pp

real(rp) sum_charge
integer ip, np, n_part
logical beta_found

! Particle number bookkeeping

np = beam_init%n_particle
if (np > 0) then
  if (np > n_part) then
    call out_io (s_warn$, r_name, &
            'Number of particles ' // int_str(n_part) // ' defined in beam file: ' // full_name, &
            'is less than the number of particles wanted which is set by beam_init%n_particle: ' // int_str(np), &
            'The setting of beam_init%n_particle will be ignored.')
  endif
  n_part = min(n_part, np)
endif

call reallocate_bunch(bunch, n_part, .true.)
p => bunch%particle

bunch%n_live = count(p%state == alive$)

! Charge bookkeeping

if (beam_init%bunch_charge /= 0) bunch%charge_tot = beam_init%bunch_charge

sum_charge = sum(p(:)%charge)
if (bunch%charge_tot == 0) then
  bunch%charge_tot = sum_charge
elseif (sum_charge == 0) then
  p%charge = bunch%charge_tot / n_part
else
  p%charge = p%charge * bunch%charge_tot / sum_charge
endif

bunch%charge_live = sum(p%charge, (p%state == alive$))

! Species

if (beam_init%species /= '') p%species = species_id(beam_init%species)

! Calc beta if able to

if (.not. beta_found) then
  do ip = 1, size(bunch%particle)
    pp => bunch%particle(ip)
    if (pp%species == not_set$ .or. pp%species == photon$ .or. pp%p0c <= 0) cycle
    call convert_pc_to ((1.0_rp + pp%vec(6)) * pp%p0c, pp%species, beta = pp%beta)
  enddo
endif

end subroutine bunch_finalizer

!---------------------------------------------------------------------------------------------------
! contains

subroutine read_switch(who, switch, str, n_line, err, ix_in)

integer switch, n_line, ix
integer, optional :: ix_in
character(*) who, str
logical err

!

if (present(ix_in)) then
  ix = ix_in
else
  call string_trim(str, str, ix)
endif

if (ix == 0) err = .true.

select case (who)
case ('location')
  call match_word(str(:ix), location_name(1:), switch)
  if (switch <= 0) then
    call out_io (s_error$, r_name, 'LOCATION NAME NOT RECOGNIZED: ' // quote(str(1:ix)), &
                                   'IN FILE: ' // file_name, 'AT LINE: ' // int_str(n_line))
    return
  endif

case ('state')
  call match_word(str(:ix), state_name, switch)
  if (switch <= 0) then
    call out_io (s_error$, r_name, 'PARTICLE STATE NAME NOT RECOGNIZED: ' // quote(str(1:ix)), &
                                   'IN FILE: ' // file_name, 'AT LINE: ' // int_str(n_line))
    return
  endif
  switch = switch - 1   ! Since state_name is zero based.

case ('species')
  switch = species_id(str(:ix), positron$)
  if (switch == invalid$) err = .true.
  if (err) then
    call out_io (s_error$, r_name, 'SPECIES NAME NOT RECOGNIZED: ' // quote(str(1:ix)), &
                                   'IN FILE: ' // file_name, 'AT LINE: ' // int_str(n_line))
    return
  endif
end select

if (present(ix_in)) then
  call string_trim(str(ix+1:), str, ix)
  ix_in = ix
endif

end subroutine read_switch

!---------------------------------------------------------------------------------------------------
! contains

subroutine read_component(name, component, line, n_line, ix, err)

real(rp) component
integer n_line, ix
character(*) name, line
logical err

!

if (ix == 0) then
  err = .true.
  call out_io (s_error$, r_name, 'LINE ENDED BEFORE COLUMN: ' // trim(name) // ' AT LINE NUMBER: ' // int_str(n_line), &
                                 'IN FILE: ' // file_name)
  return
endif

read (line, *, iostat = ios) component
err = (ios /= 0)
if (err) then
  call out_io (s_error$, r_name, 'ERROR READING REAL NUMBER: ' // quote(line(1:ix)), &
                                 'IN FILE: ' // file_name, 'IN COLUMN: ' // trim(name) // ' AT LINE NUMBER: ' // int_str(n_line))
  return
endif

call string_trim(line(ix+1:), line, ix)

end subroutine read_component

!---------------------------------------------------------------------------------------------------
! contains

subroutine read_component_int(name, component, line, n_line, ix, err)

integer component
integer n_line, ix
character(*) name, line
logical err

!

if (ix == 0) then
  err = .true.
  call out_io (s_error$, r_name, 'LINE ENDED BEFORE COLUMN: ' // trim(name) // &
                                 ' AT LINE NUMBER: ' // int_str(n_line), 'IN FILE: ' // file_name)
  return
endif

read (line, *, iostat = ios) component
if (err) then
  call out_io (s_error$, r_name, 'ERROR READING INTEGER NUMBER: ' // quote(line(1:ix)), &
                                 'IN FILE: ' // file_name, 'IN COLUMN: ' // trim(name) // ' AT LINE NUMBER: ' // int_str(n_line))
  return
endif
call string_trim(line(ix+1:), line, ix)

end subroutine read_component_int

!---------------------------------------------------------------------------------------------------
! contains

function read_param(line, n_line, err) result (param)
character(*) line
real(rp) param
integer n_line, ix, ios
logical err

!

ix = index(line, '=')
read(line(ix+1:), *, iostat = ios) param
err = (ios /= 0 .or. ix == 0) 
if (err) then
  call out_io (s_error$, r_name, 'ERROR PARSING PARAMETER: ' // line, &
                                 'IN FILE: ' // file_name, 'AT LINE NUMBER: ' // int_str(n_line))
endif

end function read_param

!---------------------------------------------------------------------------------------------------
! contains

subroutine read_params(line, n_line, param, err)
character(*) line
real(rp) param(:)
integer n_line, ix, ios
logical err

!

ix = index(line, '=')
read(line(ix+1:), *, iostat = ios) param
err = (ios /= 0 .or. ix == 0) 
if (err) then
  call out_io (s_error$, r_name, 'ERROR PARSING PARAMETER: ' // line, &
                                 'IN FILE: ' // file_name, 'AT LINE NUMBER: ' // int_str(n_line))
endif

end subroutine read_params

!---------------------------------------------------------------------------------------------------
! contains

function read_string(line, n_line, err) result (str)
character(*) line
character(200) str
integer n_line, ix, ios
logical err

!

ix = index(line, '=')
err = (ix == 0)
if (err) then
  call out_io (s_error$, r_name, 'ERROR PARSING PARAMETER: ' // line, &
                                 'IN FILE: ' // file_name, 'AT LINE NUMBER: ' // int_str(n_line))
endif

str = unquote(trim(line(ix+1:)))

end function read_string

end subroutine read_beam_ascii 

end module
