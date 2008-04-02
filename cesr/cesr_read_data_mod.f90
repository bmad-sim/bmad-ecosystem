module cesr_read_data_mod

use cesr_utils
use physical_constants

type cesr_data_params_struct
  character(20) data_date, data_type
  character(40) var_ele_name, var_attrib_name
  character(100) comment, lattice
  integer csr_set
  integer species
  real(rp) horiz_beta_freq, vert_beta_freq
  real(rp) dvar
end type

type cesr_cbar_datum_struct
  real(rp) val
  logical good
end type

type cesr_xy_datum_struct
  real(rp) x, y
  logical good
end type

type cesr_data_struct
  type (cesr_xy_datum_struct) orbit(0:120), phase(0:120), eta(0:120)
  type (cesr_cbar_datum_struct) cbar11(0:120), cbar12(0:120)
  type (cesr_cbar_datum_struct) cbar21(0:120), cbar22(0:120)
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_data_parameters (data_file, param, err)
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

subroutine read_cesr_data_parameters (data_file, param, err)

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

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_data (data_file, all_data, err)
!
! Routine to read in the orbit, phase, cbar, or eta data from a file. 
! Note: This routine does not read in butns files.
! 
! Modules needed:
!   use read_cesr_data_mod
!
! Input:
!   data_file -- Character(*): Name of the data file.
!
! Output:
!   all_data  -- Cesr_data_struct: Structure holding the data.
!   err       -- Logical: Set true if there is a read error.
!-

subroutine read_cesr_data (data_file, all_data, err)

implicit none

type (cesr_data_struct) d, all_data

integer i, j, iu, ios

character(*) data_file
character(20) dat_name
character(40) :: r_name = 'read_cesr_data'

logical err

namelist / cesr_data / d

!--------------------------------------------------------------------

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

end module
