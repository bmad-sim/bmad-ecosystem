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
integer save_set
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
! Subroutine number_to_cesr_data_file (number, who, file_name, err, print_err)
!
! Routine to form the file name for a phase, orbit, eta, or ac_dispersion
! data file. 
! If number is negative then the number of the data file will be:
!   number_of_last_data_set + number
!
! Modules needed:
!   use cesr_read_data_mod
!
! Input:
!   number    -- Integer: Number of the orbit data file.
!   who       -- Character(*): 'phase', 'orbit', 'eta', or 'ac_eta'  
!   print_err -- Logical, optional: Print an error message if there is an error?
!                 Default is True.
!
! Output:
!   file_name -- Character(*): Full file name including directory spec.
!   err       -- Logical: Set True if there is an error. False otherwise.
!-

subroutine number_to_cesr_data_file (number, who, file_name, err, print_err)

implicit none

integer number, ix_set

character(*) who, file_name
character(40) :: r_name = 'number_to_cesr_data_file'
logical err
logical, optional :: print_err

!

file_name = ''
call number_to_data_file_number (number, who, ix_set, err, print_err)
if (err) return
call form_file_name_with_number (who, ix_set, file_name, err)

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine number_to_data_file_number (num_in, who, num_out, err, print_err)
!
! Routine to form the file name for a phase, orbit, eta, or ac_dispersion
! data file. 
!   num_out = num_in                             if num_in > 0
!           = number_of_last_data_set + num_in   otherwise
!
! Modules needed:
!   use cesr_read_data_mod
!
! Input:
!   num_in    -- Integer: Number of the orbit data file.
!   who       -- Character(*): 'phase', 'orbit', 'eta', or 'ac_eta'  
!   print_err -- Logical, optional: Print an error message if there is an error?
!                 Default is True.
!
! Output:
!   file_name -- Character(*): Full file name including directory spec.
!                   Set to '' if there is an error.
!   err       -- Logical: Set True if there is an error. False otherwise.
!-

subroutine number_to_data_file_number (num_in, who, num_out, err, print_err)

implicit none

integer num_in, num_out

character(*) who
character(40) :: r_name = 'number_to_data_file_number'
logical err
logical, optional :: print_err

!

select case (who)

case ('phase')
  call calc_file_number ('$CESR_MNT/phase/phase.number', num_in, num_out, err)

case ('orbit')
  call calc_file_number ('CESR_MNT:[orbit]next_butnum.num', num_in, num_out, err)
  if (err) return
  if (num_in < 1)  num_out = num_out - 1  ! Number in file is 1 + current number

case ('eta')
  call calc_file_number ('$CESR_MNT/eta/eta.number', num_in, num_out, err)
  if (err) return

case ('ac_eta')

case default
  if (logic_option(.true., print_err)) &
                call out_io (s_error$, r_name, 'BAD "WHO": ' // who)
  return
end select




end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_orbit_data (file_name, all_dat, err)
!
! Routine to read the data in a raw orbit file.
! Note: You can use number_to_cesr_data_file to convert from an orbit number
! to a data file name.
!
! Modules needed:
!   use cesr_read_data_mod
!
! Input:
!   file_name -- Character(*): Name of the orbit data file. 
!
! Output:
!   all_dat   -- cesr_all_data_struct: Holds the data along with 
!                     data base and measurement parameters.
!-

subroutine read_cesr_orbit_data (file_name, all_dat, err)

implicit none

type (butns_struct) butns
type (cesr_all_data_struct) all_dat

character(*) file_name
character(100) dir, name
character(40) :: r_name = 'read_cesr_orbit_data'

integer i, ix, ix_in

logical err_flag, read_ok, err, is_rel

! 

call cesr_all_data_struct_init (all_dat)

err_flag = .true.

call read_butns_file (file_name, .true., butns, all_dat%db, read_ok, .true.)
if (.not. read_ok) call err_exit

ix = splitfilename (file_name, dir, name, is_rel)
all_dat%param%file_name = name
ix = index(name, '.')
if (is_integer(name(ix+1:))) read (name(ix+1:), *) all_dat%param%ix_data_set 

all_dat%param%lattice   = butns%lattice
all_dat%param%data_date = butns%date
all_dat%param%comment   = butns%comment(1)
all_dat%param%csr_set   = butns%save_set

all_dat%orbit_x%value = butns%det%x_orb
all_dat%orbit_y%value = butns%det%y_orb
all_dat%orbit_x%good  = butns%det%ok
all_dat%orbit_y%good  = butns%det%ok

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_dispersion_data (file_name, all_dat, err)
!
! Routine to read the data in a raw dispersion file.
! Note: You can use number_to_cesr_data_file to convert from an orbit number
! to a data file name.
!
! Modules needed:
!   use cesr_read_data_mod
!
! Input:
!   file_name -- Character(*): Name of the orbit data file. 
!
! Output:
!   all_dat   -- Cesr_all_data_struct: Holds the data along with 
!                     data base and measurement parameters.
!-

subroutine read_cesr_dispersion_data (file_name, all_dat, err)

implicit none

type (cesr_all_data_struct) all_dat
type (cesr_xy_data_struct) eta_(0:120)

integer ix, iu, ios

logical err

character(*) file_name
character(100) dir, name
character(40) :: r_name = 'read_cesr_dispersion_data'

logical is_rel

namelist / dispersion_data / eta_

! Init

call cesr_all_data_struct_init (all_dat)

! open the file and read the contents

err = .true.

iu = lunget()
open (iu, file = file_name, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // file_name)
  return
endif

! Read header info

call read_cesr_data_parameters (iu, all_dat, err)
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

ix = splitfilename (file_name, dir, name, is_rel)
all_dat%param%file_name = name
ix = index(name, '.')
if (is_integer(name(ix+1:))) read (name(ix+1:), *) all_dat%param%ix_data_set 

all_dat%eta_x%value = eta_%x
all_dat%eta_y%value = eta_%y
all_dat%eta_x%good  = eta_%good
all_dat%eta_y%good  = eta_%good

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_phase_data (file_name, all_dat, err)
!
! Routine to read the data in a raw phase/coupling file.
! Note: You can use number_to_cesr_data_file to convert from an orbit number
! to a data file name.
!
! Modules needed:
!   use cesr_read_data_mod
!
! Input:
!   file_name -- Character(*): Name of the orbit data file. 
!
! Output:
!   all_dat   -- Cesr_all_data_struct: Holds everything from the data
!                     file: phase, cbar, and raw orbit data along with 
!                     data base and measurement parameters.
!-

subroutine read_cesr_phase_data (file_name, all_dat, err)

implicit none

type (cesr_all_data_struct) all_dat
type (phase_cbar_data_struct) pc_(0:120)
type (detector_struct) orbit_(0:120)
type (db_struct) db
type (raw_det_struct) :: h_(0:120), v_(0:120)

real(rp) horiz_beta_freq, vert_beta_freq
real(rp) horiz_reflection_shake, vert_reflection_shake

integer ix, ios, iu, species

character(*) file_name
character(100) name, dir
character(40) :: r_name = 'read_cesr_phase_data'

logical err, is_rel

namelist / phase_parameters / species, horiz_beta_freq, vert_beta_freq, &
          horiz_reflection_shake, vert_reflection_shake
namelist / phase_cbar_data / pc_
namelist / rawdata / h_, v_
namelist / raworbit / orbit_

! Init

call cesr_all_data_struct_init (all_dat)

! Open the file and read the contents

err = .true.

iu = lunget()
open (iu, file = file_name, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // file_name)
  return
endif

! Read header info

call read_cesr_data_parameters (iu, all_dat, err)
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

ix = splitfilename (file_name, dir, name, is_rel)
all_dat%param%file_name = name
ix = index(name, '.')
if (is_integer(name(ix+1:))) read (name(ix+1:), *) all_dat%param%ix_data_set 

all_dat%param%horiz_beta_freq = horiz_beta_freq
all_dat%param%vert_beta_freq  = vert_beta_freq
all_dat%param%species = species

! Read in raw data

read (iu, nml = rawdata)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'WARNING: ERROR READING "RAWDATA" NAMELIST')
  rewind(iu)
endif

all_dat%raw_phase%a = h_
all_dat%raw_phase%b = v_

! read the raw orbit

read (iu, nml = raworbit, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING RAWORBIT')
  return
endif

all_dat%raw_orbit = orbit_

! read the phase and cbar data

pc_(:)%ok_x = .false.
pc_(:)%ok_y = .false.
read (iu, nml = phase_cbar_data, iostat = ios)
close (iu)

if (ios /= 0) then
  call out_io (s_error$, r_name, 'ERROR READING PHASE_CBAR_DATA')
  return
endif

all_dat%phase_a%value = pc_%x_phase * twopi / 360
all_dat%phase_a%good  = pc_%ok_x
all_dat%phase_b%value = pc_%y_phase * twopi / 360
all_dat%phase_b%good  = pc_%ok_y

all_dat%cbar11_b%value = pc_%y_cbar11
all_dat%cbar11_b%good  = pc_%ok_y
all_dat%cbar12_a%value = pc_%x_cbar12
all_dat%cbar12_a%good  = pc_%ok_x
all_dat%cbar12_b%value = pc_%y_cbar12
all_dat%cbar12_b%good  = pc_%ok_y
all_dat%cbar22_a%value = pc_%x_cbar22
all_dat%cbar22_a%good  = pc_%ok_x

err = .false.

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine read_cesr_cooked_data (file_name, all_dat, err)
!
! Routine to read in the orbit, phase, cbar, or eta cooked data from a file. 
! 
! Modules needed:
!   use read_cesr_data_mod
!
! Input:
!   file_name -- Character(*): Name of the data file.
!
! Output:
!   all_data  -- Cesr_all_data_struct: Structure holding the data.
!   err       -- Logical: Set true if there is a read error.
!-

subroutine read_cesr_cooked_data (file_name, all_dat, err)

implicit none

type xy_data_input_struct
  real(rp) x, y
end type

type (cesr_all_data_struct) all_dat
type (xy_data_input_struct) orbit(0:120), phase(0:120), eta(0:120)
real(rp) cbar11(0:120), cbar12_x(0:120), cbar12_y(0:120), cbar22(0:120)

integer i, j, iu, ios

character(*) file_name
character(20) dat_name
character(40) :: r_name = 'read_cesr_cooked_data'

logical err

namelist / cesr_data / orbit, phase, eta, cbar11, cbar12_x, cbar12_y, cbar22

!--------------------------------------------------------------------

call cesr_all_data_struct_init (all_dat)

call read_cesr_cooked_data_parameters (file_name, all_dat%param, err)
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
open (iu, file = file_name, status = 'old', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // file_name)
  return
endif

read (iu, nml = cesr_data, iostat = ios)
close (iu)
if (ios /= 0) then       ! abort on read error
  call out_io (s_error$, r_name, 'ERROR READING DATA FROM: ' // file_name)
  return
endif

call transfer_data (orbit%x,  all_dat%orbit_x)
call transfer_data (orbit%y,  all_dat%orbit_y)
call transfer_data (phase%x,  all_dat%phase_a)
call transfer_data (phase%y,  all_dat%phase_b)
call transfer_data (eta%x,    all_dat%eta_x)
call transfer_data (eta%y,    all_dat%eta_y)
call transfer_data (cbar11,   all_dat%cbar11_b)
call transfer_data (cbar12_x, all_dat%cbar12_a)
call transfer_data (cbar12_y, all_dat%cbar12_b)
call transfer_data (cbar22,   all_dat%cbar22_a)

all_dat%param%file_name = file_name

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
! Subroutine read_cesr_cooked_data_parameters (file_name, param, err)
!
! Routine to read in the header information from a data file.
! 
! Modules needed:
!   use read_cesr_data_mod
!
! Input:
!   file_name -- Character(*): Name of the data file.
!
! Output:
!   param     -- Cesr_data_params_struct: Parameters 
!   err       -- Logical: Set true if there is a read error.
!-

subroutine read_cesr_cooked_data_parameters (file_name, param, err)

implicit none

type (cesr_data_params_struct) param

integer ios, iu

character(*) file_name
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
open (iu, file = file_name, status = 'old', iostat = ios)
if (ios /= 0) then       ! abort on open error
  call out_io (s_error$, r_name, 'ERROR OPENING: ' // file_name)
  return
endif

read (iu, nml = data_parameters, iostat = ios)
if (ios /= 0) then
  call out_io (s_fatal$, r_name, &
          'ERROR READING "DATA_PARAMETERS" NAMELIST IN FILE: ' // file_name)
  rewind (iu)
  read (iu, nml = data_parameters)
endif

close (iu)
err = .false.

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine read_butns_file (file_name, nonlinear_calc, butns, db, read_ok, type_err)
!
! Subroutine to read the raw information from a BUTNS.nnnnn file and convert
! it to orbit x-y positions.
!
! Modules Needed:
!   use cesr_read_data_mod
!
! Input:
!   file_name      -- Character(*): Name of the orbit data file. 
!   nonlinear_calc -- Logical: Calculate orbit using Rich Helms nonlinear 
!                       routines? Recomendation: True.
!   type_err       -- Logical: If True then open error message will be printed.
!                         If False then open error message will be supressed.
!
! Output:
!   butns -- Butns_struct: Orbit information.
!       %lattice    -- Character(40): Lattice name.
!       %save_set   -- Integer: Save set number.
!       %date       -- Character(20): Date orbit was taken
!       %turn       -- Integer: Turn number for injection orbits. 0 otherwise.
!       %comment(5) -- Character(72): Comment.
!       %det(0:99)%ok     -- Logical: Was there a valid orbit reading?
!       %det(0:99)%amp(4) -- Integer: raw button numbers.
!       %det(0:99)%x_orb  -- Real(rp): Horizontal orbit in meters.
!       %det(0:99)%y_orb  -- Real(rp): Horizontal orbit in meters.
!   db    -- Db_struct: Structure holding the steering settings.
!     %csr_horz_cur(i)%cu_now -- CU settings for CSR HORZ CUR
!     %csr_hbnd_cur(i)%cu_now -- CU settings for CSR HBND CUR
!     %csr_vert_cur(i)%cu_now -- CU settings for CSR VERT CUR
!     %csr_hsp_volt(i)%cu_now -- CU settings for CSR HSP VVAL
!     %csr_vsp_volt(i)%cu_now -- CU settings for CSR VSP VOLT
!     %scir_vertcur(i)%cu_now -- CU settings for SCIR VERTCUR
!     %scir_pos_stp(i)%cu_now -- CU settings for SCIR POS STP
!     %scir_enc_cnt(i)%cu_now -- CU settings for SCIR ENC CNT
!     %scir_pos_rd(i)%cu_now  -- CU settings for SCIR POS RD
!   read_ok -- Logical: Set True if butns file was successfuly parsed.
!
! Note: orbit numbers are from 0 to 99
!
! Note: db%csr_hsp_volt is actually obtained from the node CSR HSP VVAL which
! records the actual voltage as opposed to the command. Since CSR HSP VVAL
! is a readback node the values put in db%csr_hsp_volt will not be exactly the
! actual commands (and if the separators have been turned off they will not
! even be approximately the same).
!-

subroutine read_butns_file (file_name, nonlinear_calc, butns, db, &
                                                          read_ok, type_err)

  implicit none

  type (db_struct) db
  type (butns_struct) butns

  integer vec(120), det_type(120)
  integer i, ix, j, iu, ios, raw(4, 120)
  integer n_node, n_ele

  real x_orbit(120), y_orbit(120), rdummy
  character(*) file_name

  character(130) line_in

  logical nonlinear_calc, read_ok, type_err, err_flag, is_ok(120)

! init

  butns%comment = ' '
  read_ok = .false.

! read header line in the data file (to get lattice name)
! and sterlat strengths

  iu = lunget()
  open(unit = iu, file = file_name, status = 'old', action = 'READ', iostat = ios)
  if (ios /= 0) then
    if (type_err) print *, &
          'ERROR IN READ_BUTNS_FILE: ERROR OPENING: ', trim(file_name)
    return
  endif

  read (iu, '(a)') line_in
  butns%lattice = line_in(61:)
  call string_trim (butns%lattice, butns%lattice, ix)

  butns%date = line_in(30:)                     ! get date
  read (line_in(54:), *) butns%save_set

  do
    read (iu, '(a)', iostat = ios) line_in
    if (ios /= 0) then
      print *, 'ERROR IN ORBIT_READ: ERROR READING STEERINGS IN ORBIT FILE'
      goto 1000
    endif
    ix = index(line_in, 'END BUTNS')
    if (ix /= 0) exit
    ix = index(line_in, 'TURN=')
    if (ix /= 0) then
      read (line_in(ix+5:), *, iostat = ios) butns%turn
      if (ios /= 0) then
        print *, 'ERROR IN READ_BUTNS_FILE: ERROR READING TURN NUMBER.'
        print *, '      IN FILE: ', trim(file_name)
      endif
    endif
  enddo

  read (line_in(ix+10:), *, iostat = ios) n_node
  if (ios /= 0) then
    print *, 'ERROR IN ORBIT_READ: ERROR READING STEERING NODE NUMBER IN ORBIT FILE'
    goto 1000
  endif

! read data base valuse stored in the orbit file

  do i = 1, n_node

    read (iu, '(a)', iostat = ios) line_in
    if (ios /= 0) then
      print *, 'ERROR IN ORBIT_READ: ERROR READING STEERING NODE IN ORBIT FILE'
      exit
    endif

    read (line_in(14:), *, iostat = ios) n_ele

    if (line_in(2:13) == 'CSR COMMENTS') then
      do j = 1, n_ele
        read (iu, '(1x, a)', iostat = ios) butns%comment(j)
        if (ios /= 0) then
          print *, 'ERROR IN ORBIT_READ: ERROR READING COMMENT #', j
          exit
        endif
      enddo
      cycle
    endif

    read (iu, '(12x, 10i6)') vec(1:n_ele)

    if (line_in(2:13) == 'CSR HORZ CUR') then
      db%csr_horz_cur(1:n_ele)%cu_now = vec(1:n_ele)
      db%csr_horz_cur(1:n_ele)%valid_cu_now = .true.
    elseif (line_in(2:13) == 'CSR VERT CUR') then
      db%csr_vert_cur(1:n_ele)%cu_now = vec(1:n_ele)
      db%csr_vert_cur(1:n_ele)%valid_cu_now = .true.
    elseif (line_in(2:13) == 'CSR HBND CUR') then
      db%csr_hbnd_cur(1:n_ele)%cu_now = vec(1:n_ele)
      db%csr_hbnd_cur(1:n_ele)%valid_cu_now = .true.
    elseif (line_in(2:13) == 'CSR HSP VVAL') then
      call hsp_vval_to_volt (vec, vec)
      n_ele = size(db%csr_hsp_volt)
      db%csr_hsp_volt(1:n_ele)%cu_now = vec(1:n_ele)
      db%csr_hsp_volt(1:n_ele)%valid_cu_now = .true.
    elseif (line_in(2:13) == 'CSR VSP VOLT') then
      db%csr_vsp_volt(1:n_ele)%cu_now = vec(1:n_ele)
      db%csr_vsp_volt(1:n_ele)%valid_cu_now = .true.
    elseif (line_in(2:13) == 'SCIR VERTCUR') then
      db%scir_vertcur(1:n_ele)%cu_now = vec(1:n_ele)
      db%scir_vertcur(1:n_ele)%valid_cu_now = .true.
    elseif (line_in(2:13) == 'SCIR POS STP') then
      db%scir_pos_stp(1:n_ele)%cu_now = vec(1:n_ele)
      db%scir_pos_stp(1:n_ele)%valid_cu_now = .true.
    elseif (line_in(2:13) == 'SCIR ENC CNT') then
      db%scir_enc_cnt(1:n_ele)%cu_now = vec(1:n_ele)
      db%scir_enc_cnt(1:n_ele)%valid_cu_now = .true.
    elseif (line_in(2:13) == 'SCIR POS RD ') then
      db%scir_pos_rd(1:n_ele)%cu_now = vec(1:n_ele)
      db%scir_pos_rd(1:n_ele)%valid_cu_now = .true.
    else
      print *, 'ERROR IN ORBIT_READ: UNKNOWN NODE IN ORBIT FILE: ', line_in(2:13)
      goto 1000
    endif

  enddo

! close

  1000 continue
  close (unit = iu)

! read in the raw data

  call butfilget (raw, file_name, rdummy, det_type)      ! read in raw data
  if (nonlinear_calc) then
    call nonlin_butcon (raw, 1, 100, y_orbit, x_orbit)
  else
    call butcon (raw, 1, 100, y_orbit, x_orbit)
  endif

  is_ok = .false.
  call det_ok (raw, 1, 100, det_type, is_ok)

  do i = 1, 100
    j = i
    if (i == 100) j = 0
    butns%det(j)%amp = raw(1:4, i)
    butns%det(j)%x_orb = x_orbit(i) / 1000.0   ! convert to m
    butns%det(j)%y_orb = y_orbit(i) / 1000.0   ! convert to m
    butns%det(j)%ok = is_ok(i)
  end do

  read_ok = .true.

end subroutine

end module
