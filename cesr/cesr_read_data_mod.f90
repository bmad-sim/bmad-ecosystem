module cesr_read_data_mod

use cesr_basic_mod
use cesr_db_mod

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine read_cesr_data_parameters (iu, all_dat, err)

implicit none

type (cesr_all_data_struct) all_dat

integer i, iu, ios
integer save_set, species
integer orbit_num1, orbit_num2, orbit_num
integer csr_quad_cur(98), csr_qadd_cur(n_qadd_maxx), csr_sext_cur(98)
integer csr_sqewquad(98), csr_sqewsext(n_csr_sqewsext_maxx), csr_horz_cur(98)
integer csr_hbnd_cur(n_hbnd_maxx), csr_vert_cur(98)
integer csr_hsp_volt(n_sep_maxx), csr_vsp_volt(n_sep_maxx), csr_hsp_vval(12)
integer csr_octu_cur(n_oct_maxx), scir_quadcur(n_scir_quad_maxx)
integer scir_skqucur(n_scir_quad_maxx), scir_vertcur(n_scir_quad_maxx)
integer ir_sksxcur(1), scir_pos_stp(n_scir_cam_maxx)
integer scir_enc_cnt(n_scir_cam_maxx), scir_pos_rd(3*n_scir_cam_maxx)
integer nir_shuntcur(n_nir_shuntcur_maxx)
integer scwig_contrl(100), scw_cur_read(100)  !<---- FIX THIS!
integer scsol_contrl(5*n_sc_sol_maxx)

character(80) comment
character(40) route_name
character(20) data_date, file_type
character(40) lattice
character(40) :: r_name = 'read_cesr_data_parameters'

logical err

namelist / data_parameters / data_date, lattice, comment, &
          file_type, save_set, route_name

namelist / data_base / &
     csr_quad_cur, csr_qadd_cur, csr_sext_cur, csr_sqewquad, csr_sqewsext, &
     csr_horz_cur, csr_hbnd_cur, csr_vert_cur, csr_hsp_volt, csr_vsp_volt, &
     csr_hsp_vval, csr_octu_cur, scir_quadcur, scir_skqucur, scir_vertcur, &
     ir_sksxcur, scir_pos_stp, scir_enc_cnt, scir_pos_rd , scwig_contrl, &
     scw_cur_read, nir_shuntcur, scsol_contrl

! Init

err = .true.
rewind (iu)

route_name = ''
nir_shuntcur = 0

read (iu, nml = data_parameters, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING "DATA_PARAMETERS" NAMELIST')
  return
endif

all_dat%param%data_date  = data_date 
all_dat%param%lattice    = lattice
all_dat%param%comment    = comment
all_dat%param%csr_set    = save_set
all_dat%param%route_name = route_name

if (file_type(:5) == 'PHASE') then
  all_dat%param%data_type = 'PHASE'
elseif (file_type(:10) == 'DISPERSION') then
  all_dat%param%data_type = 'ETA'
else
  all_dat%param%data_type = file_type
endif

read (iu, nml = data_base, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING "DATA_BASE" NAMELIST IN DATA FILE.')
  return
endif

all_dat%db%csr_quad_cur%cu_now = csr_quad_cur
all_dat%db%csr_qadd_cur%cu_now = csr_qadd_cur
all_dat%db%nir_shuntcur%cu_now = nir_shuntcur
all_dat%db%csr_sext_cur%cu_now = csr_sext_cur
all_dat%db%csr_sqewquad%cu_now = csr_sqewquad
all_dat%db%csr_scsolcur%cu_now = scsol_contrl
all_dat%db%csr_sqewsext%cu_now = csr_sqewsext
all_dat%db%csr_horz_cur%cu_now = csr_horz_cur
all_dat%db%csr_hbnd_cur%cu_now = csr_hbnd_cur
all_dat%db%csr_vert_cur%cu_now = csr_vert_cur
call hsp_vval_to_volt (csr_hsp_vval, all_dat%db%csr_hsp_volt%cu_now)
all_dat%db%csr_vsp_volt%cu_now = csr_vsp_volt
all_dat%db%csr_octu_cur%cu_now = csr_octu_cur
all_dat%db%scir_quadcur%cu_now = scir_quadcur
all_dat%db%scir_skqucur%cu_now = scir_skqucur
all_dat%db%scir_vertcur%cu_now = scir_vertcur
all_dat%db%ir_sksxcur%cu_now = ir_sksxcur
all_dat%db%scir_pos_stp(:)%cu_now = scir_pos_stp
all_dat%db%scir_pos_rd(:)%cu_now  = scir_pos_rd
all_dat%db%scir_enc_cnt(:)%cu_now = scir_enc_cnt

all_dat%db%csr_quad_cur%valid_cu_now = .true.
all_dat%db%csr_qadd_cur%valid_cu_now = .true.
all_dat%db%nir_shuntcur%valid_cu_now = .true.
all_dat%db%csr_sext_cur%valid_cu_now = .true.
all_dat%db%csr_sqewquad%valid_cu_now = .true.
all_dat%db%csr_scsolcur%valid_cu_now = .true.
all_dat%db%csr_sqewsext%valid_cu_now = .true.
all_dat%db%csr_horz_cur%valid_cu_now = .true.
all_dat%db%csr_hbnd_cur%valid_cu_now = .true.
all_dat%db%csr_vert_cur%valid_cu_now = .true.
all_dat%db%csr_hsp_volt%valid_cu_now = .true.
all_dat%db%csr_vsp_volt%valid_cu_now = .true.
all_dat%db%csr_octu_cur%valid_cu_now = .true.
all_dat%db%scir_quadcur%valid_cu_now = .true.
all_dat%db%scir_skqucur%valid_cu_now = .true.
all_dat%db%scir_vertcur%valid_cu_now = .true.
all_dat%db%ir_sksxcur%valid_cu_now = .true.
all_dat%db%scir_pos_stp(:)%valid_cu_now = .true.
all_dat%db%scir_pos_rd(:)%valid_cu_now  = .true.
all_dat%db%scir_enc_cnt(:)%valid_cu_now = .true.

err = .false.

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_orbit_data (orbit_number, data, err)
!
! Routine to read the data in a raw orbit file.
!
! Modules needed:
!   use cesr_read_data_mod
!
! Input:
!   orbit_number -- Integer: Number of the orbit data file.
!
! Output:
!   data   -- All_data_struct: Holds the data along with 
!                     data base and measurement parameters.
!-

subroutine read_cesr_orbit_data (orbit_number, data, err)

implicit none

type (butns_struct) butns
type (cesr_all_data_struct) data

character(40) :: r_name = 'read_cesr_orbit_data'

integer i, ix_in, ix_set, orbit_number

logical err_flag, read_ok, err

! 

call cesr_all_data_struct_init (data)

err_flag = .true.
call calc_file_number ('CESR_MNT:[orbit]next_butnum.num', orbit_number, ix_set, err)
if (err) call err_exit
if (orbit_number < 1)  ix_set = ix_set - 1  ! Number in file is 1 + current number

call read_butns_file (ix_set, .true., butns, data%db, read_ok, .true.)
if (.not. read_ok) call err_exit
data%param%ix_data_set = ix_set

data%param%lattice   = butns%lattice
data%param%data_date = butns%date
data%param%comment   = butns%comment(1)
data%param%csr_set   = butns%save_set

data%orbit_x%value = butns%det%x_orb
data%orbit_y%value = butns%det%y_orb
data%orbit_x%good  = butns%det%ok
data%orbit_y%good  = butns%det%ok

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_dispersion_data (eta_number, data, err)
!
! Routine to read the data in a raw dispersion file.
!
! Modules needed:
!   use cesr_read_data_mod
!
! Input:
!   eta_number -- Integer: Number of the eta data file.
!
! Output:
!   data   -- All_data_struct: Holds the data along with 
!                     data base and measurement parameters.
!-

subroutine read_cesr_dispersion_data (eta_number, data, err)

implicit none

type (cesr_all_data_struct) data
type (cesr_xy_data_struct) eta_(0:120)

integer eta_number
integer ix, iu, ios

logical err

character(40) :: r_name = 'read_cesr_dispersion_data'
character(100) file_name

namelist / dispersion_data / eta_

! Init

call cesr_all_data_struct_init (data)

! first construct the file name

call calc_file_number ('$CESR_MNT/eta/eta.number', eta_number, ix, err)
if (err) return
call form_file_name_with_number ('ETA', ix, file_name, err)
if (err) return

! open the file and read the contents

err = .true.

iu = lunget()
open (iu, file = file_name, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // file_name)
  return
endif

! Read header info

call read_cesr_data_parameters (iu, data, err)
if (err) then
  close (iu)
  return
endif

! Read phase parameters

rewind (iu)
eta_%good = .false.
eta_%x = 0
eta_%y = 0
read (iu, nml = dispersion_data)
close (iu)

data%eta_x%value = eta_%x
data%eta_y%value = eta_%y
data%eta_x%good  = eta_%good
data%eta_y%good  = eta_%good

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_phase_data (phase_number, data, err)
!
! Routine to read the data in a raw phase/coupling file.
!
! Modules needed:
!   use cesr_read_data_mod
!
! Input:
!   phase_number -- Integer: Number of the phase data file.
!
! Output:
!   data   -- Cesr_all_data_struct: Holds everything from the data
!                     file: phase, cbar, and raw orbit data along with 
!                     data base and measurement parameters.
!-

subroutine read_cesr_phase_data (phase_number, data, err)

implicit none

type (cesr_all_data_struct) data
type (phase_cbar_data_struct) pc_(0:120)
type (detector_struct) orbit_(0:120)
type (db_struct) db
type (raw_det_struct) :: h_(0:120), v_(0:120)

real(rp) horiz_beta_freq, vert_beta_freq
real(rp) horiz_reflection_shake, vert_reflection_shake

integer phase_number
integer ix, ios, iu, species

character(100) file_name
character(40) :: r_name = 'read_cesr_phase_data'

logical err

namelist / phase_parameters / species, horiz_beta_freq, vert_beta_freq, &
          horiz_reflection_shake, vert_reflection_shake
namelist / phase_cbar_data / pc_
namelist / rawdata / h_, v_
namelist / raworbit / orbit_

! Init

call cesr_all_data_struct_init (data)

! read a data file...
! first construct the file name

call calc_file_number ('$CESR_MNT/phase/phase.number', phase_number, ix, err)
if (err) return
call form_file_name_with_number ('PHASE', ix, file_name, err)
if (err) return

! open the file and read the contents

err = .true.

iu = lunget()
open (iu, file = file_name, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // file_name)
  return
endif

! Read header info

call read_cesr_data_parameters (iu, data, err)
if (err) then
  close (iu)
  return
endif

! Read phase parameters

rewind (iu)

read (iu, nml = phase_parameters, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'WARNING: ERROR READING "PHASE_PARAMETERS" NAMELIST')
  rewind(iu)
endif

data%param%horiz_beta_freq = horiz_beta_freq
data%param%vert_beta_freq  = vert_beta_freq
data%param%species = species

! read the raw orbit

read (iu, nml = raworbit, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING RAWORBIT')
  return
endif

data%raw_orbit = orbit_

! read the phase and cbar data

pc_(:)%ok_x = .false.
pc_(:)%ok_y = .false.
read (iu, nml = phase_cbar_data, iostat = ios)

if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING PHASE_CBAR_DATA')
  return
endif

data%phase_x%value = pc_%x_phase * twopi / 360
data%phase_x%good  = pc_%ok_x
data%phase_y%value = pc_%y_phase * twopi / 360
data%phase_y%good  = pc_%ok_y

data%cbar11_y%value = pc_%y_cbar11
data%cbar11_y%good  = pc_%ok_y
data%cbar12_x%value = pc_%x_cbar12
data%cbar12_x%good  = pc_%ok_x
data%cbar12_y%value = pc_%y_cbar12
data%cbar12_y%good  = pc_%ok_y
data%cbar22_x%value = pc_%x_cbar22
data%cbar22_x%good  = pc_%ok_x

err = .false.

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_cooked_data (data_file, all_data, err)
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
!   all_data  -- Cesr_all_data_struct: Structure holding the data.
!   err       -- Logical: Set true if there is a read error.
!-

subroutine read_cesr_cooked_data (data_file, all_data, err)

implicit none

type xy_data_input_struct
  real(rp) x, y
end type

type (cesr_all_data_struct) all_data
type (xy_data_input_struct) orbit(0:120), phase(0:120), eta(0:120)
real(rp) cbar11(0:120), cbar12_x(0:120), cbar12_y(0:120), cbar22(0:120)

integer i, j, iu, ios

character(*) data_file
character(20) dat_name
character(40) :: r_name = 'read_cesr_cooked_data'

logical err

namelist / cesr_data / orbit, phase, eta, cbar11, cbar12_x, cbar12_y, cbar22

!--------------------------------------------------------------------

call cesr_all_data_struct_init (all_data)

call read_cesr_cooked_data_parameters (data_file, all_data%param, err)
if (err) return

orbit%x  = real_garbage$
orbit%y  = real_garbage$
phase%x  = real_garbage$
phase%y  = real_garbage$
eta%x    = real_garbage$
eta%y    = real_garbage$
cbar11   = real_garbage$
cbar12_x = real_garbage$
cbar12_y = real_garbage$
cbar22   = real_garbage$

err = .true.

iu = lunget()
open (iu, file = data_file, status = 'old', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // data_file)
  return
endif

read (iu, nml = cesr_data, iostat = ios)
close (iu)
if (ios /= 0) then       ! abort on read error
  call out_io (s_error$, r_name, 'ERROR READING DATA FROM: ' // data_file)
  return
endif

call transfer_data (orbit%x,  all_data%orbit_x)
call transfer_data (orbit%y,  all_data%orbit_y)
call transfer_data (phase%x,  all_data%phase_x)
call transfer_data (phase%y,  all_data%phase_y)
call transfer_data (eta%x,    all_data%eta_x)
call transfer_data (eta%y,    all_data%eta_y)
call transfer_data (cbar11,   all_data%cbar11_y)
call transfer_data (cbar12_x, all_data%cbar12_x)
call transfer_data (cbar12_y, all_data%cbar12_y)
call transfer_data (cbar22,   all_data%cbar22_x)

err = .false.

!--------------------------------------------------------
contains

subroutine transfer_data (from, to)

real(rp) from(:)
type (cesr_data1_struct) to(:)

integer i

!

do i = lbound(from, 1), ubound(from, 1)
  if (from(i) == real_garbage$) then
    to(i)%good = .false.
    to(i)%value = 0
  else
    to(i)%good = .true.
    to(i)%value = from(i)
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
character(40) :: r_name = 'read_cesr_cooked_data_parameters'

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
