!+
! Module sodom_mod
! 
! Routines used by the sodom program
!+

module sodom2_mod
use bmad
use sim_utils
use mode3_mod 
!use twiss_and_track_mod
use beam_mod
use f95_lapack
!use twiss_and_track_mod
!use ptc_map_with_radiation_mod
!use high_energy_space_charge_mod
!use superimpose_mod, only: add_superimpose
!use radiation_mod
!use expression_mod
!use mode3_mod, only: make_N
!use random_mod
!use s_fitting_new, only: probe, internal_state, track_probe, assignment(=), operator(+), default, spin0

implicit none

!type ltt_column_struct
!  character(120) :: param = ''
!  character(40) :: header_str = ''
!  character(20) :: format = ''
!end type
integer, parameter :: master_rank$  = 0
integer, parameter :: a_mode$ = 1
integer, parameter :: b_mode$ = 2
integer, parameter :: c_mode$ = 3
integer, parameter :: angle_max$ = 40
integer, private :: i_loop

! User settable parameters

type sodom2_params_struct
  character(40) :: ele_eval = ''
  character(200) :: lat_file = ''
  character(200) :: output_file = 'sodom2.out' ! to store n_axes for inputted action/angles
  integer :: mode = -1 ! mode to excite orbital motion
  real(rp) :: J = -1
  real(rp) :: angles_to_eval(angle_max$) = [0.0_rp, (-10.0_rp, i_loop = 2, angle_max$)]  ! outside of [0,2*pi] do not evaluate
  integer :: n_particle = -1
end type

! Common vars

type sodom2_com_struct
  type (lat_struct) :: lat
  !type (coord_struct), allocatable :: bmad_closed_orb(:)
  !type (normal_modes_struct) modes
  !type (beam_init_struct) beam_init
  type (bunch_struct) :: bunch
  integer :: ix_ele_eval
  integer :: ix_branch = 0                   ! Lattice branch being tracked.
  real(rp) :: time_start = 0
  logical :: debug = .false.
  !integer :: n_particle      ! Num particles per bunch. Needed with MPI.
  integer :: mpi_rank = master_rank$
  !integer :: mpi_run_index = 0                 ! Run index
  !integer :: mpi_ix0_particle = 0              ! Index of first particle
  logical :: using_mpi = .false.
  character(200) :: master_input_file = ''
  real(rp) :: Q_2pi
  real(rp) :: ADST = 0
  real(rp) :: n_axis(3, angle_max$) = 0
  real(rp) :: n_mat(6,6)
  complex(rp), allocatable :: U11(:)
  complex(rp), allocatable :: U12(:)
  complex(rp), allocatable :: U21(:)
  complex(rp), allocatable :: U22(:)
  complex(rp), allocatable :: M(:,:)
  complex(rp), allocatable :: eig_vec(:,:)
  complex(rp), allocatable :: eig_val(:)
  integer :: idx_max_0 = 0 ! Index of eigenvector with max |psi_0|
  !character(40) :: ps_fmt = '(2i7, 8es16.8, 3x, 3f10.6, 4x, a)'
end type

integer, parameter :: new$ = 0,  valid$ = 1, written$ = 2

type (sodom2_params_struct), pointer, save :: sodom2_params_global   ! Needed for track1_preprocess and track1_bunch_hook
type (sodom2_com_struct),    pointer, save :: sodom2_com_global      ! Needed for track1_preprocess and track1_bunch_hook

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_read_params(sodom, sodom2_com)

type (sodom2_params_struct), target :: sodom
type (sodom2_com_struct), target :: sodom2_com

! will need to generate a beam_init 
! type (beam_init_struct) beam_init

integer i, ix
character(200) arg
character(40) m_name
character(*), parameter :: r_name = 'sodom2_read_params'

namelist / params / sodom ! input file is fortran namelist format, can easily load into this struct

sodom2_params_global => sodom
sodom2_com_global => sodom2_com
sodom2_com%master_input_file = ''

! Parse command line

i = 0
do while (i < command_argument_count())
  i = i + 1
  call get_command_argument(i, arg)
  call match_word (arg, ['-debug'], ix, .false., .true., m_name)
  select case (m_name)
  case ('-debug')
    sodom2_com%debug = .true.
  case default
    if (sodom2_com%master_input_file /= '') then
      print '(2a)', 'Extra stuff on the command line: ', quote(arg)
      print '(a)',  'Stopping here.'
      stop
    endif
    sodom2_com%master_input_file = arg
  end select
end do

if (sodom2_com%master_input_file == '') sodom2_com%master_input_file = 'sodom.init'

! Read parameters

if (.not. sodom2_com%using_mpi .or. sodom2_com%mpi_rank == master_rank$) then
  print '(2a)', 'Initialization file: ', trim(sodom2_com%master_input_file)
endif

open (1, file = sodom2_com%master_input_file, status = 'old', action = 'read')
read (1, nml = params)
close (1)

call bmad_parser (sodom%lat_file, sodom2_com%lat)

print *,'Mode = ', sodom%mode

! Check:
if (sodom%J < 0) then
  call out_io (s_fatal$, r_name, 'sodom%J MUST BE GREATER THAN ZERO. STOPPING HERE.')
  stop
endif

if (sodom%mode < 0 .or. sodom%mode > 3) then
  call out_io (s_fatal$, r_name, 'sodom%mode MUST BE 1, 2, OR 3. STOPPING HERE.')
  stop
endif

if (sodom%n_particle < 1) then
  call out_io (s_fatal$, r_name, 'sodom%n_particle MUST BE GREATER THAN 1. STOPPING HERE.')
  stop
endif

!i = 1
!do while(sodom%angles_to_eval(i) /= -10)
!  if (sodom%angles_to_eval(i) > 2*acos(-1.0_rp) .or. sodom%angles_to_eval(i) < 0 ) then
!    call out_io(s_fatal$, r_name, 'Angles in sodom%angles_to_eval MUST BE WITHIN [0, 2*pi]. STOPPING HERE.')
!  endif
!  i = i + 1
!enddo

end subroutine sodom2_read_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_init_params(sodom, sodom2_com)

type (sodom2_params_struct), target :: sodom
type (sodom2_com_struct), target :: sodom2_com
type (lat_struct) :: lat
!type (branch_struct):: branch
type (ele_pointer_struct), allocatable :: eles(:)

integer n_loc, i
logical err

lat = sodom2_com%lat
bmad_com%auto_bookkeeper = .false.
bmad_com%spin_tracking_on = .true.

! make sure RF is on so 6D action-angle variables can be calculated
print '(2a)', 'Note: RF is always set ON for SODOM-2'
call set_on_off (rfcavity$, lat, on$)

if (sodom%ele_eval /= '') then
  call lat_ele_locator (sodom%ele_eval, lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    print '(2a)', 'Evaluate element not found: ', trim(sodom%ele_eval)
    stop
  endif
  if (n_loc > 1) then
    print '(2a)', 'Multiple elements found with ele_eval name: ', trim(sodom%ele_eval)
    print '(a)', 'Will stop here.'
    stop
  endif
  !branch => pointer_to_branch(eles(1)%ele)
  sodom2_com%ix_branch = eles(1)%ele%ix_branch
  eles(1)%ele%type = '@START_ELE' ! from ltt - why?
  sodom2_com%ix_ele_eval = eles(1)%ele%ix_ele
else
  sodom2_com%ix_branch = 0
  sodom2_com%ix_ele_eval = lat%ele(0)%ix_ele
endif

!branch => lat%branch(sodom2_com%ix_branch)

! set tracking method for each element to linear$
do i = 1, lat%branch(sodom2_com%ix_branch)%n_ele_track
  call set_ele_attribute(lat%branch(sodom2_com%ix_branch)%ele(i), "tracking_method = linear", err, .true., .true.)! , .false.)
enddo

call lattice_bookkeeper(lat)

call fullfilename (sodom%output_file, sodom%output_file)

! allocate SU2 quaternion discrete function
allocate(sodom2_com%U11(sodom%n_particle))
allocate(sodom2_com%U12(sodom%n_particle))
allocate(sodom2_com%U21(sodom%n_particle))
allocate(sodom2_com%U22(sodom%n_particle))

! allocate 2N x 2N matrix
allocate(sodom2_com%M(2*sodom%n_particle, 2*sodom%n_particle))

! allocate eig vecs/vals
allocate(sodom2_com%eig_vec(2*sodom%n_particle, 2*sodom%n_particle))
allocate(sodom2_com%eig_val(2*sodom%n_particle))

end subroutine sodom2_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_init_bunch(sodom, sodom2_com)

type (sodom2_params_struct), target :: sodom
type (sodom2_com_struct), target :: sodom2_com
type (bunch_struct) :: bunch
type (lat_struct) :: lat

real(rp) m(6,6), vec0(6), phi, J1, J2, PI, vec(6)
real(rp) :: jvec(6) = 0
real(rp) :: spin1(3) = (/1, 0, 0/)
real(rp) :: spin2(3) = (/0, 1, 0/)
real(rp) :: spin3(3) = (/0, 0, 1/)
complex(rp) :: eig_val(6)
complex(rp) :: eig_vec(6,6)

integer ix_branch
logical error
integer i

PI = acos(-1.0_rp)

lat = sodom2_com%lat
! allocate memory for bunch
call reallocate_bunch(bunch, 3*sodom%n_particle)

call transfer_matrix_calc (lat, m, ix1 = sodom2_com%ix_ele_eval, ix_branch = sodom2_com%ix_branch, one_turn = .true.) !, sodom2_com%ele_eval%ix_ele, sodom2_com%ele_eval%ix_branch, .true.)
call mat_eigen(m, eig_val, eig_vec, error, .false.)
!call twiss3_at_start(lat, error, sodom2_com%ix_branch, tune3)
!sodom2_com%Q_2pi = tune3(sodom%mode)*2*PI
sodom2_com%Q_2pi = real(log(eig_val(sodom%mode*2))*(0.0_rp, 1.0_rp))

call make_N(m, sodom2_com%n_mat, error)
if (error) then
  print *, 'Cannot compute 1-turn matrix eigen vectors. Stopping here.'
  stop
endif

do i = 1, sodom%n_particle
  phi = (2d0 * PI * (i - 1)) / sodom%n_particle
  J1 = sqrt(2*sodom%J)*cos(phi)
  J2 = -sqrt(2*sodom%J)*sin(phi)
  jvec(2*sodom%mode-1) = J1
  jvec(2*sodom%mode) = J2

  vec = matmul(sodom2_com%n_mat, jvec)
 

  ! initialize 3 particles w/ orthogonal spins
  call init_coord(bunch%particle(3*i-2), vec, sodom2_com%lat%ele(sodom2_com%ix_ele_eval), downstream_end$, sodom2_com%lat%param%particle, 1, spin = spin1)
  call init_coord(bunch%particle(3*i-1), vec, sodom2_com%lat%ele(sodom2_com%ix_ele_eval), downstream_end$, sodom2_com%lat%param%particle, 1, spin = spin2)
  call init_coord(bunch%particle(3*i), vec, sodom2_com%lat%ele(sodom2_com%ix_ele_eval), downstream_end$, sodom2_com%lat%param%particle, 1, spin = spin3)
enddo

sodom2_com%bunch = bunch

end subroutine sodom2_init_bunch

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_track_bunch(sodom, sodom2_com)

type (sodom2_params_struct), target :: sodom
type (sodom2_com_struct), target :: sodom2_com
logical err_flag
call track_bunch(sodom2_com%lat, sodom2_com%bunch, sodom2_com%lat%ele(sodom2_com%ix_ele_eval),sodom2_com%lat%ele(sodom2_com%ix_ele_eval), err_flag )!sodom2_com%ele_eval, sodom2_com%ele_eval, err_flag)

end subroutine sodom2_track_bunch

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine quat_to_SU2(q,U)

real(rp) q(0:3)
complex(rp) U(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)
complex(rp) :: eye_2(2,2) = reshape((/ (1,0), (0,0), (0,0), (1,0) /), [2,2])
sigma_1 = reshape((/ (0,0), (1,0), (1,0), (0,0) /), shape(sigma_2))
sigma_2 = reshape((/ (0,0), (0,1), (0,-1), (0,0) /), shape(sigma_2))
sigma_3 = reshape((/ (1,0), (0,0), (0,0), (-1,0) /), shape(sigma_2))

U = q(0)*eye_2+(0,-1)*(q(1)*sigma_1+q(2)*sigma_2+q(3)*sigma_3)

end subroutine quat_to_SU2

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_construct_quaternions(sodom, sodom2_com)
type (sodom2_params_struct), target :: sodom
type (sodom2_com_struct), target :: sodom2_com
complex(rp) U(2,2)

real(rp) R(3,3), q(0:3)
integer i

do i = 1, sodom%n_particle
  R(:,1) = sodom2_com%bunch%particle(3*i-2)%spin
  R(:,2) = sodom2_com%bunch%particle(3*i-1)%spin
  R(:,3) = sodom2_com%bunch%particle(3*i)%spin

  q = w_mat_to_quat(R)
  call quat_to_SU2(q,U)
  sodom2_com%U11(i) = U(1,1)
  sodom2_com%U12(i) = U(1,2)
  sodom2_com%U21(i) = U(2,1)
  sodom2_com%U22(i) = U(2,2)
enddo


end subroutine sodom2_construct_quaternions

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_construct_mat(sodom, sodom2_com)
type (sodom2_params_struct), target :: sodom
type (sodom2_com_struct), target :: sodom2_com

integer j, k, n, idx

logical err_flag

call fft_1d(sodom2_com%U11, -1)
call fft_1d(sodom2_com%U12, -1)
call fft_1d(sodom2_com%U21, -1)
call fft_1d(sodom2_com%U22, -1)

! center harmonics around 0:
sodom2_com%U11 = cshift(sodom2_com%U11, -(sodom%n_particle-1)/2)/sodom%n_particle
sodom2_com%U21 = cshift(sodom2_com%U21, -(sodom%n_particle-1)/2)/sodom%n_particle
sodom2_com%U12 = cshift(sodom2_com%U12, -(sodom%n_particle-1)/2)/sodom%n_particle
sodom2_com%U22 = cshift(sodom2_com%U22, -(sodom%n_particle-1)/2)/sodom%n_particle

do j = -(sodom%n_particle - 1)/2, (sodom%n_particle - 1)/2 
  do k = -(sodom%n_particle - 1)/2, (sodom%n_particle - 1)/2 
    
    ! enforce circular indexing of harmonics
    idx = modulo((j-k)+(sodom%n_particle-1)/2, sodom%n_particle)+1
    !M(2*j+n_particles:2*j+n_particles+1,2*k+n_particles:2*k+n_particles+1) = exp(-1i*2*pi*j*Q)*A(j-k);
    if (idx == 0) then
      print *, 'j = ', j, ', k = ', k
    endif
    sodom2_com%M(2*j+sodom%n_particle,2*k+sodom%n_particle) = exp(-(0,1)*j*sodom2_com%Q_2pi)*sodom2_com%U11(idx)
    sodom2_com%M(2*j+sodom%n_particle,2*k+sodom%n_particle+1) = exp(-(0,1)*j*sodom2_com%Q_2pi)*sodom2_com%U12(idx)
    sodom2_com%M(2*j+sodom%n_particle+1,2*k+sodom%n_particle) = exp(-(0,1)*j*sodom2_com%Q_2pi)*sodom2_com%U21(idx)
    sodom2_com%M(2*j+sodom%n_particle+1,2*k+sodom%n_particle+1) = exp(-(0,1)*j*sodom2_com%Q_2pi)*sodom2_com%U22(idx)
  enddo
enddo

end subroutine sodom2_construct_mat

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_deallocate_memory(sodom2_com)

type (sodom2_com_struct), target :: sodom2_com

deallocate(sodom2_com%U11)
deallocate(sodom2_com%U12)
deallocate(sodom2_com%U21)
deallocate(sodom2_com%U22)
deallocate(sodom2_com%M)
deallocate(sodom2_com%eig_vec)
deallocate(sodom2_com%eig_val)
call reallocate_bunch(sodom2_com%bunch, 0)

end subroutine sodom2_deallocate_memory

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_eig(sodom2_com)
type (sodom2_com_struct), target :: sodom2_com

logical err

call la_geev(sodom2_com%M, sodom2_com%eig_val, VL = sodom2_com%eig_vec) !, err, .true.)


end subroutine sodom2_eig

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_calc_ADST(sodom, sodom2_com)
type (sodom2_params_struct), target :: sodom
type (sodom2_com_struct), target :: sodom2_com
integer j, n_particle, max_idx
real(rp) :: PI
real(rp), parameter::inv_pi = 1.0_rp / acos(-1.0_rp)

PI = acos(-1.0_rp)

n_particle = sodom%n_particle

! determine eigenvector with largest |psi_0|
max_idx = 1
do j = 1, 2*n_particle
    if (real(sodom2_com%eig_vec(n_particle,j)*conjg(sodom2_com%eig_vec(n_particle,j))) > real(sodom2_com%eig_vec(n_particle, max_idx)*conjg(sodom2_com%eig_vec(n_particle, max_idx)))) then
        max_idx = j
    endif
enddo

sodom2_com%idx_max_0 = max_idx
sodom2_com%ADST = real(log(sodom2_com%eig_val(max_idx))*(0.0_rp, inv_pi))
print *,'ADST = ', sodom2_com%ADST

end subroutine sodom2_calc_ADST

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_calc_n(sodom, sodom2_com)
type (sodom2_params_struct), target :: sodom
type (sodom2_com_struct), target :: sodom2_com
complex(rp) sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)
complex(rp) :: psi_n(2)
real(rp) :: phi
integer :: i = 1
integer j


sigma_1 = reshape((/ (0,0), (1,0), (1,0), (0,0) /), shape(sigma_1))
sigma_2 = reshape((/ (0,0), (0,1), (0,-1), (0,0) /), shape(sigma_2))
sigma_3 = reshape((/ (1,0), (0,0), (0,0), (-1,0) /), shape(sigma_3))

do while(sodom%angles_to_eval(i) /= -10)
  phi = sodom%angles_to_eval(i)
  psi_n = (/ (0,0), (0,0) /)
  do j = -(sodom%n_particle-1)/2, (sodom%n_particle-1)/2
    psi_n = psi_n + sodom2_com%eig_vec(2*j+sodom%n_particle:2*j+sodom%n_particle+1, sodom2_com%idx_max_0)*exp((0,1)*j*phi);
  enddo
  
  sodom2_com%n_axis(1,i) = dot_product(psi_n, matmul(sigma_1, psi_n))
  sodom2_com%n_axis(2,i) = dot_product(psi_n, matmul(sigma_2, psi_n))
  sodom2_com%n_axis(3,i) = dot_product(psi_n, matmul(sigma_3, psi_n))
  sodom2_com%n_axis(1:3,i) = sodom2_com%n_axis(1:3,i)/norm2(sodom2_com%n_axis(1:3,i))

  !print *, "phi = ", sodom%angles_to_eval(i),":"
  !print *, "n_x = ", sodom2_com%n_axis(1,i)
  !print *, "n_y = ", sodom2_com%n_axis(2,i)
  !print *, "n_z = ", sodom2_com%n_axis(3,i)
  !print *, "-----------------------------------------------"
  i = i + 1
enddo


end subroutine sodom2_calc_n

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_write(sodom, sodom2_com)
type (sodom2_params_struct), target :: sodom
type (sodom2_com_struct), target :: sodom2_com
integer iu
integer :: ia = 1
real(rp) phi, J1, J2
real(rp) :: jvec(6) = 0

iu = lunget()

open (iu, file = sodom%output_file, recl = 300)
write (iu,  '(3a)')     	'# lat_file                            	= ', sodom%lat_file !quote(sodom2_com%lat)
write (iu, '(a, a)')   		'# ele_eval                   		= ', sodom%ele_eval
write (iu,  '(a, i8)')  	'# n_particle                          	= ', sodom%n_particle
write (iu,  '(a, i8)')  	'# mode                          	= ', sodom%mode
write (iu,  '(a, 8es16.8)')  	'# J                          		= ', sodom%J
write (iu,  '(a)') 	'# -----------------------------------------------------------'
write (iu,  '(a, 8es16.8)') 	'# ADST                        		= ', sodom2_com%ADST
write (iu,  '(a)') 	'# -----------------------------------------------------------'
write (iu, '(a)')  '## phi [rad]	 |          x              px               y              py               z              pz        |         nx              ny              nz'

do while(sodom%angles_to_eval(ia) /= -10)
  jvec = 0.0
  phi = sodom%angles_to_eval(ia)
  J1 = sqrt(2*sodom%J)*cos(phi)
  J2 = -sqrt(2*sodom%J)*sin(phi)
  jvec(2*sodom%mode-1) = J1
  jvec(2*sodom%mode) = J2

  jvec = matmul(sodom2_com%n_mat, jvec)

  write (iu, '(1es16.8, 3x, 6es16.8, 3x, 3es16.8)') phi, jvec, sodom2_com%n_axis(1,ia), sodom2_com%n_axis(2,ia), sodom2_com%n_axis(3,ia)
  ia = ia + 1
enddo

close(iu)

end subroutine sodom2_write

end module
