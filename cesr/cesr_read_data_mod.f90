module cesr_read_data_mod

use cesr_basic_mod
use cesr_db_mod

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! 
!
!
!-

subroutine read_cesr_phase_data (phase_number, &
                      phase_cbar, raw_orbit, db, param, err)

implicit none

type (phase_cbar_data_struct) phase_cbar, pc_(0:120)
type (detector_struct) orbit_(0:120), raw_orbit
type (db_struct) db
type (cesr_data_params_struct) param

integer phase_number
integer ix, ios, iu

character(100) file_name
character(40) :: r_name = 'read_cesr_phase_data'

logical err

namelist / phase_cbar_data / pc_
namelist / raworbit / orbit_

!

! read a data file...
! first construct the file name

call calc_file_number ('CESR_MNT:[phase]phase.number', phase_number, ix, err)
if (err) return
call form_file_name_with_number ('PHASE', ix, file_name, err)
if (err) return

! open the file and read the contents

err = .true.

iu = lunget()
open (iu, file = file_name, type = 'old', readonly, iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // file_name)
  close (iu)
  return
endif

! call read_header (iu, data_or_ref, u%phase, u, graph, logic)
rewind (iu)

read (iu, nml = raworbit, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING RAWORBIT')
  return
endif

pc_(:)%ok_x = .false.
pc_(:)%ok_y = .false.
read (iu, nml = phase_cbar_data, iostat = ios)

if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING PHASE_CBAR_DATA')
  return
endif


err = .false.

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_cooked_data (data_file, all_data, param, err)
!
! Routine to read in the orbit, phase, cbar, or eta cooked data from a file. 
! 
! Modules needed:
!   use read_cesr_data_mod
!
! Input:
!   data_file -- Character(*): Name of the data file.
!
! Output:
!   all_data  -- Cesr_data_struct: Structure holding the data.
!   param     -- Cesr_data_params_struct: Parameters 
!   err       -- Logical: Set true if there is a read error.
!-

subroutine read_cesr_cooked_data (data_file, all_data, param, err)

implicit none

type (cesr_data_struct) d, all_data
type (cesr_data_params_struct) param

integer i, j, iu, ios

character(*) data_file
character(20) dat_name
character(40) :: r_name = 'read_cesr_data'

logical err

namelist / cesr_data / d

!--------------------------------------------------------------------

call read_cesr_cooked_data_parameters (data_file, param, err)
if (err) return

d%orbit%x  = real_garbage$
d%orbit%y  = real_garbage$
d%phase%x  = real_garbage$
d%phase%y  = real_garbage$
d%eta%x    = real_garbage$
d%eta%y    = real_garbage$
d%cbar11%val = real_garbage$
d%cbar12%val = real_garbage$
d%cbar21%val = real_garbage$
d%cbar22%val = real_garbage$

err = .true.

iu = lunget()
open (iu, file = data_file, status = 'old', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // data_file)
  return
endif

read (iu, nml = cesr_data)

call mark_good_data (d%orbit%y, d%orbit%good, d%orbit%x)
call mark_good_data (d%phase%y, d%phase%good, d%phase%x)
call mark_good_data (d%eta%y, d%eta%good, d%eta%x)
call mark_good_data (d%cbar11%val, d%cbar11%good)
call mark_good_data (d%cbar12%val, d%cbar12%good)
call mark_good_data (d%cbar21%val, d%cbar21%good)
call mark_good_data (d%cbar22%val, d%cbar22%good)

d%orbit%x = d%orbit%x / 1000.0     ! convert to m
d%orbit%y = d%orbit%y / 1000.0     ! convert to m

close (iu)
err = .false.

all_data = d

!--------------------------------------------------------
contains

subroutine mark_good_data (value, good, value2)

real(rp) value(:)
real(rp), optional :: value2(:)
logical good(:)

integer i

!

do i = lbound(value, 1), ubound(value, 1)
  if (value(i) == real_garbage$) then
    good(i) = .false.
    value(i) = 0
    if (present(value2)) value2(i) = 0
  else
    good(i) = .true.
  endif
end do

end subroutine
end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_cooked_data_parameters (data_file, param, err)
!
! Routine to read in the header information from a data file.
! 
! Modules needed:
!   use read_cesr_data_mod
!
! Input:
!   data_file -- Character(*): Name of the data file.
!
! Output:
!   param     -- Cesr_data_params_struct: Parameters 
!   err       -- Logical: Set true if there is a read error.
!-

subroutine read_cesr_cooked_data_parameters (data_file, param, err)

implicit none

type (cesr_data_params_struct) param

integer ios, iu

character(*) data_file
character(40) :: r_name = 'read_cesr_data_parameters'

namelist / data_parameters / param

logical err

! Init

param%var_attrib_name = ''
param%var_ele_name    = ''
param%data_type       = ''
param%dvar            = 0
param%csr_set         = 0
param%lattice         = ''
param%species         = 0
param%horiz_beta_freq = 0
param%vert_beta_freq  = 0
param%data_date       = ''
param%comment         = ''

! Read parameters from the file

err = .true.

iu = lunget()
open (iu, file = data_file, status = 'old', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // data_file)
  return
endif

read (iu, nml = data_parameters, iostat = ios)
if (ios /= 0) then
  call out_io (s_fatal$, r_name, &
          'ERROR READING "DATA_PARAMETERS" NAMELIST IN FILE: ' // data_file)
  rewind (iu)
  read (iu, nml = data_parameters)
endif

close (iu)
err = .false.

end subroutine

end module
