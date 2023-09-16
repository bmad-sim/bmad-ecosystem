!+
! Subroutine tao_hook_command (cmd_line, found)
!
! Put custom Tao commands here. These commands are searched before the standard
! tao commands are searched. This allows for the defining of new commands and 
! the overwriting of any standard tao command.
!
! This file is already set up so that it is rather simple to add a command. Just
! follow the directions. Keep in mind that you don't have to use the included 
! infrastructure if you feel like doing something really custom -- Just clear
! everything out of here and replace it with whatever you like. All that Tao 
! really cares about is that found is set to either TRUE or FALSE.
!
! Input:
!   cmd_line   -- Character(*): command line
!   found      -- Logical: Set True if the command is handled by this routine.
!                          If false, then the standard commands will be searched
!                          and Tao will spit out an error if the command is not
!                          found
!-

subroutine tao_hook_command (command_line, found)

use tao_command_mod, dummy => tao_hook_command

implicit none

type (tao_d2_data_array_struct), allocatable :: d2(:)
type (tao_d1_data_array_struct), allocatable, target :: d1_amp_arr(:), d1_phase_arr(:)
type (tao_d1_data_struct), pointer :: d1_phase, d1_amp
type (tao_data_struct), pointer :: datum_amp, datum_phase
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat

real(rp) r1, r2, design_tune_a, design_tune_b, data_tune_a, data_tune_b, offset, offset_best, rms, rms_best
real(rp), allocatable :: design(:), data(:)

integer i, ix, iu, ix_line, ix_cmd, which, ios
integer plane, n

character(*) :: command_line
character(140) :: cmd_line, line, file_name
character(140) :: cmd_word(12), ele_name
character(16) cmd_name, mode, data_or_ref
character(*), parameter :: r_name = 'tao_hook_command'

logical found
logical, allocatable :: ok(:)

!!!! put your list of hook commands in here. 

character(16) :: cmd_names(1) = [character(16):: 'pingread']  

logical quit_tao, err

! found will be set to TRUE if the command is found in here

found = .false.

! strip the command line of comments

call string_trim (command_line, cmd_line, ix_line)
ix = index(cmd_line, '!')
if (ix /= 0) cmd_line = cmd_line(:ix-1)        ! strip off comments

! blank line => nothing to do

if (cmd_line(1:1) == '') return

! match first word to a command name
! If not found then found = .false.

call match_word (cmd_line(:ix_line), cmd_names, ix_cmd, .true., .true., cmd_name)
if (ix_cmd == 0 .or. ix_line < 3) return   ! To be recognized, command must use at least 3 characters
if (ix_cmd < 0) then
  call out_io (s_error$, r_name, 'AMBIGUOUS HOOK COMMAND')
  found = .true.
  return
endif

found = .true.
call string_trim (cmd_line(ix_line+1:), cmd_line, ix_line)

!----------------------------------------------------------
! PUT YOUR CUSTOM COMMANDS IN THIS CASE CONSTRUCT
! Look at the file tao/code/tao_command.f90 to see how the standard commands are parsed.

select case (cmd_name)

!--------------------------------
! pingread - Read ping data from a file.
! Format: 
!   pingread <mode> <filename> <data_or_ref>
! where: 
!   <mode> = "a_mode" or "b_mode"
!   <data_or_ref> = "data", "reference", or blank (defaults to "data") 

case ('pingread')
 
  ! First parse the pingread command to extract the <mode>, <filename> and <data_or_ref> arguments.

  call tao_cmd_split (cmd_line, 3, cmd_word, .true., err); if (err) return

  call match_word(cmd_word(1), [character(16):: 'a_mode', 'b_mode'], ix, .true., .true., mode)
  if (ix < 1) then
    call out_io (s_fatal$, r_name, 'Syntax:  pingread <mode> <filename> <data_or_ref>', &
                                   '<mode> must be "a_mode" or "b_mode"')
    return
  endif

  if (cmd_word(3) == '') cmd_word(3) = 'data'
  call match_word(cmd_word(3), [character(16):: 'data', 'reference'], ix, .true., .true., data_or_ref)
  if (ix < 1) then
    call out_io (s_fatal$, r_name, 'Syntax:  pingread <mode> <filename> <data_or_ref>', &
                                   '<data_or_ref> must be "data", "reference", or blank')
    return
  endif

  call tao_open_file (cmd_word(2), iu, file_name, s_fatal$)
  if (iu == 0) return

  ! Read header from the data file

  do
    read (iu, '(a)', iostat = ios) line
    if (ios < 0) then
      call out_io (s_fatal$, r_name, 'END OF FILE ENCOUNTERED BEFORE DATA ENCOUNTERED: ' // file_name)
      return
    endif
    if (ios > 0) then
      call out_io (s_fatal$, r_name, 'ERROR READING FILE: ' // file_name)
      return
    endif
    if (line(1:2) == 'R:') exit
    call string_trim(line, line, ix)

    ! Read tunes.

    if (line(1:4) == 'Tune') then
      call tao_find_data (err, 'tune', d2_array = d2)
      if (size(d2) /= 1) then
        call out_io (s_fatal$, r_name, 'NO TUNE D2 DATA STRUCTURE DEFINED!')
        return
      endif

      u => s%u(d2(1)%d2%ix_uni)
      lat => u%design%lat
      design_tune_a = lat%ele(lat%n_ele_track)%a%phi / twopi
      design_tune_b = lat%ele(lat%n_ele_track)%b%phi / twopi

      ix = index(line, '(')
      line = line(ix+1:)
      ix = index(line, ')')
      line(ix:ix) = ' '
      read (line, *, iostat = ios) data_tune_a

      ix = index(line, '(')
      line = line(ix+1:)
      ix = index(line, ')')
      line(ix:ix) = ' '
      read (line, *, iostat = ios) data_tune_b

      if (data_or_ref == 'data') then
        d2(1)%d2%d1(1)%d(1)%meas_value = twopi * (data_tune_a + nint(design_tune_a))
        d2(1)%d2%d1(1)%d(1)%good_meas = .true.
        d2(1)%d2%d1(2)%d(1)%meas_value = twopi * (data_tune_b + nint(design_tune_b))
        d2(1)%d2%d1(2)%d(1)%good_meas = .true.
      else
        d2(1)%d2%d1(1)%d(1)%ref_value = twopi * (data_tune_a + nint(design_tune_a))
        d2(1)%d2%d1(1)%d(1)%good_ref = .true.
        d2(1)%d2%d1(2)%d(1)%ref_value = twopi * (data_tune_b + nint(design_tune_b))
        d2(1)%d2%d1(2)%d(1)%good_ref = .true.
      endif      
    endif

  enddo

  ! If beginning of data then setup pointers

  if (line(3:3) == 'H') then
    if (mode == 'a_mode') then
      call tao_find_data (err, 'ping_a.amp_x', d1_array = d1_amp_arr)
      call tao_find_data (err, 'ping_a.phase_x', d1_array = d1_phase_arr)
    else 
      call tao_find_data (err, 'ping_b.amp_x', d1_array = d1_amp_arr)
      call tao_find_data (err, 'ping_b.phase_x', d1_array = d1_phase_arr)
    endif
  elseif (line(3:3) == 'V') then
    if (mode == 'a_mode') then
      call tao_find_data (err, 'ping_a.amp_y', d1_array = d1_amp_arr)
      call tao_find_data (err, 'ping_a.phase_y', d1_array = d1_phase_arr)
    else 
      call tao_find_data (err, 'ping_b.amp_y', d1_array = d1_amp_arr)
      call tao_find_data (err, 'ping_b.phase_y', d1_array = d1_phase_arr)
    endif
  else
    call out_io (s_fatal$, r_name, 'CANNOT IDENTIFY BPM PLANE!')
    return    
  endif

  if (size(d1_amp_arr) /= 1 .or. size(d1_phase_arr) /= 1) then
    call out_io (s_fatal$, r_name, 'NO PING DATA STRUCTURE(S) DEFINED!')
    return
  endif

  d1_amp => d1_amp_arr(1)%d1
  d1_phase => d1_phase_arr(1)%d1
  allocate (data(size(d1_phase%d)), design(size(d1_phase%d)), ok(size(d1_phase%d)))
  ok = .false.

  ! Read data

  do
    if (line == '') cycle

    call tao_cmd_split (line, 4, cmd_word, .false., err)
    read (cmd_word(2), *) r1
    read (cmd_word(3), *) r2
    ele_name = cmd_word(1)
    datum_amp => tao_pointer_to_datum(d1_amp, ele_name(3:))
    datum_phase => tao_pointer_to_datum(d1_phase, ele_name(3:))
    if (.not. associated(datum_amp)) then
      call out_io (s_fatal$, r_name, 'ELEMENT NAME IN FILE DOES NOT HAVE A CORRESPONDING DATUM: ' // ele_name(3:))
      return
    endif

! Fix off-by-one problem in the data
!    if (mode == 'a_mode') then
!      if (line(3:3) == 'H') then
!        if (datum_phase%ix_d1 < 15 .or. any(datum_phase%ix_d1 == [100, 101, 102])) then
!          r1 = r1 - data_tune_a
!        endif
!      endif
!    else
!      if (line(3:3) == 'V') then
!        if (datum_phase%ix_d1 < 15 .or. any(datum_phase%ix_d1 == [100, 101, 102])) then
!          r1 = r1 - data_tune_b
!        endif
!      endif
!    endif

    if (data_or_ref == 'data') then
      datum_phase%good_meas = .true.
      datum_amp%meas_value = r2
      datum_amp%good_meas = .true.
    else
      datum_phase%good_ref = .true.
      datum_amp%ref_value = r2
      datum_amp%good_ref = .true.
    endif

    data(datum_phase%ix_d1) = r1
    design(datum_phase%ix_d1) = datum_phase%design_value / twopi
    ok(datum_phase%ix_d1) = .true.

    ! Next line

    read (iu, '(a)', iostat = ios) line
    if (ios < 0) exit
    if (ios > 0) then
      call out_io (s_fatal$, r_name, 'ERROR READING FILE: ' // file_name)
      return
    endif

  enddo

  ! Shift the measured phases until all of them are within +/- 180 deg of design.
  ! Must take into account that there is an overall arbitrary offset between
  ! the design and the measured. 

  n = count(ok)

  rms_best = 1e30

  do i = 1, 20
    offset = i / 20.0
    data = data + nint(design + offset - data)
    rms = sum((data - design - offset)**2, mask = ok)
    if (rms < rms_best) then
      offset_best = offset
      rms_best = rms
    endif
  enddo

  data = data + nint(design + offset_best - data)

  ! Put in data struct

  do i = 1, size(data)
    if (data_or_ref == 'data') then
      d1_phase%d(i)%meas_value = twopi * data(i)
    else
      d1_phase%d(i)%ref_value = twopi * data(i)
    endif
  enddo

end select

! Do the standard calculations and plotting after command execution.
! See the Tao manual section enititled "Lattice Calculation" in the 
! "Overview" chapter for details on what is performed here.

call tao_cmd_end_calc

end subroutine tao_hook_command
