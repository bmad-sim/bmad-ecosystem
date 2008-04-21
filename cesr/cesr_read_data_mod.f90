module cesr_read_data_mod

use cesr_basic_mod
use cesr_db_mod

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine read_cesr_data_parameters (iu, pc_dat, err)

implicit none

type (all_phase_cbar_data_struct) pc_dat

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

!

err = .true.
rewind (iu)

route_name = ''
nir_shuntcur = 0

read (iu, nml = data_parameters, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING "DATA_PARAMETERS" NAMELIST')
  return
endif

pc_dat%param%data_date  = data_date 
pc_dat%param%lattice    = lattice
pc_dat%param%comment    = comment
pc_dat%param%csr_set    = save_set
pc_dat%param%route_name = route_name

if (file_type(:5) == 'PHASE') then
  pc_dat%param%data_type = 'PHASE'
elseif (file_type(:10) == 'DISPERSION') then
  pc_dat%param%data_type = 'ETA'
else
  pc_dat%param%data_type = file_type
endif

read (iu, nml = data_base, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING "DATA_BASE" NAMELIST IN DATA FILE.')
  return
endif

pc_dat%db%csr_quad_cur%cu_now = csr_quad_cur
pc_dat%db%csr_qadd_cur%cu_now = csr_qadd_cur
pc_dat%db%nir_shuntcur%cu_now = nir_shuntcur
pc_dat%db%csr_sext_cur%cu_now = csr_sext_cur
pc_dat%db%csr_sqewquad%cu_now = csr_sqewquad
pc_dat%db%csr_scsolcur%cu_now = scsol_contrl
pc_dat%db%csr_sqewsext%cu_now = csr_sqewsext
pc_dat%db%csr_horz_cur%cu_now = csr_horz_cur
pc_dat%db%csr_hbnd_cur%cu_now = csr_hbnd_cur
pc_dat%db%csr_vert_cur%cu_now = csr_vert_cur
call hsp_vval_to_volt (csr_hsp_vval, pc_dat%db%csr_hsp_volt%cu_now)
pc_dat%db%csr_vsp_volt%cu_now = csr_vsp_volt
pc_dat%db%csr_octu_cur%cu_now = csr_octu_cur
pc_dat%db%scir_quadcur%cu_now = scir_quadcur
pc_dat%db%scir_skqucur%cu_now = scir_skqucur
pc_dat%db%scir_vertcur%cu_now = scir_vertcur
pc_dat%db%ir_sksxcur%cu_now = ir_sksxcur
pc_dat%db%scir_pos_stp(:)%cu_now = scir_pos_stp
pc_dat%db%scir_pos_rd(:)%cu_now  = scir_pos_rd
pc_dat%db%scir_enc_cnt(:)%cu_now = scir_enc_cnt

err = .false.

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
!   data   -- All_phase_cbar_data_struct: Holds everything from the data
!                     file: phase, cbar, and raw orbit data along with 
!                     data base and measurement parameters.
!-

subroutine read_cesr_phase_data (phase_number, data, err)

implicit none

type (all_phase_cbar_data_struct) data
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

! read a data file...
! first construct the file name

call calc_file_number ('CESR_MNT:[phase]phase.number', phase_number, ix, err)
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

data%phase_cbar = pc_

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
