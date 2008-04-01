module read_cesr_data_mod

use cesr_utils

type cesr_data_params_struct
  character(20) data_date, data_type
  character(40) var_ele_name, var_attrib_name
  character(100) comment, lattice
  integer csr_set
  integer species
  real(rp) horiz_beta_freq, vert_beta_freq
  real(rp) dvar
end type

type cesr_fake_orbit_struct
  real(rp) x_orbit(0:120), y_orbit(0:120)
  logical  good(0:120)
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_data_parameters (data_file, data_params, err)
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
!   data_params -- Cesr_data_params_struct: Parameters 
!   err         -- Logical: Set true if there is a read error.
!-

subroutine read_cesr_data_parameters (data_file, data_params, err)

implicit none

type (cesr_data_params_struct) data_params

real(rp) dvar

integer ios, iu, csr_set, species

real(rp) horiz_beta_freq, vert_beta_freq

character(*) data_file
character(20) data_date, data_type
character(40) var_ele_name, var_attrib_name
character(100) comment, line, lattice
character(40) :: r_name = 'read_cesr_data_parameters'

namelist / data_parameters / data_date, lattice, species, &
          horiz_beta_freq, vert_beta_freq, csr_set, data_type, comment, &
          var_ele_name, var_attrib_name, dvar

logical err

! Init

data_params%var_attrib_name = ''
data_params%var_ele_name    = ''
data_params%data_type       = ''
data_params%dvar            = 0
data_params%csr_set         = 0
data_params%lattice         = ''
data_params%species         = 0
data_params%horiz_beta_freq = 0
data_params%vert_beta_freq  = 0
data_params%data_date       = ''
data_params%comment         = ''

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
  return
endif

close (iu)

! transfer them to the data_parameters structure

data_params%var_attrib_name = var_attrib_name
data_params%var_ele_name    = var_ele_name
data_params%data_type       = data_type
data_params%dvar            = dvar
data_params%csr_set         = csr_set
data_params%lattice         = lattice
data_params%species         = species
data_params%horiz_beta_freq = horiz_beta_freq
data_params%vert_beta_freq  = vert_beta_freq
data_params%data_date       = data_date
data_params%comment         = comment

err = .false.

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_fake_orbit_data (data_file, orbit, err)
!
! Routine to read in the orbit from a data file. 
! Note: This routine does not read in butns files.
! 
! Modules needed:
!   use read_cesr_data_mod
!
! Input:
!   data_file -- Character(*): Name of the data file.
!
! Output:
!   orbit     -- Cesr_fake_orbit_struct: Structure holding the orbit data.
!   err       -- Logical: Set true if there is a read error.
!-

subroutine read_cesr_fake_orbit_data (data_file, orbit, err)

implicit none

type (cesr_fake_orbit_struct) orbit

integer i, j, iu, ios

real(rp) x_orbit, y_orbit

character(*) data_file
character(100) line
character(20) dat_name
character(40) :: r_name = 'tc_read_fake_orbit_data'

logical err

!--------------------------------------------------------------------

orbit%x_orbit = 0
orbit%y_orbit = 0
orbit%good = .false.

iu = lunget()
open (iu, file = data_file, status = 'old', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // data_file)
  return
endif

call skip_header(iu, err)
if (err) return

do while (.true.)
  read (iu, *, iostat = ios) i, x_orbit, y_orbit
  if (ios < 0) exit   ! end-of-file
  if (ios > 0) then   ! read error
    backspace (iu)
    read (iu, '(a)') line
    call out_io (s_error$, r_name, 'ERROR READING ORBIT DATA IN: ' // data_file, &
                                   'BAD LINE: ' // line)
    err = .true.
    close (iu)
    return
  endif
  j = i
  if (i == 100) j = 0
  orbit%x_orbit(i) = x_orbit / 1000.0     ! convert to m
  orbit%y_orbit(i) = y_orbit / 1000.0     ! convert to m
  orbit%good(i) = .true.
enddo

close (iu)

err = .false.

end subroutine

end module
