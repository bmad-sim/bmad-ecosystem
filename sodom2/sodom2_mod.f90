!+
! Module sodom_mod
! 
! Routines used by the sodom program
!+

module sodom2_mod
use bmad
use sim_utils
use mode3_mod 
use beam_mod
use f95_lapack

implicit none

integer, parameter :: master_rank$  = 0
integer, parameter :: a_mode$ = 1
integer, parameter :: b_mode$ = 2
integer, parameter :: c_mode$ = 3
integer, parameter :: angle_max$ = 1000
integer, private :: i_loop

! User settable parameters

type sodom2_params_struct
  character(40) :: ele_eval = ''
  character(200) :: lat_file = ''
  character(200) :: output_file = 'sodom2.out' ! to store n_axes for inputted action/angles
  integer :: mode = -1 ! mode to excite orbital motion
  real(rp) :: J = -1
  real(rp) :: angles_to_eval(angle_max$) = [(-10.0_rp, i_loop = 1, angle_max$)]
  integer :: n_particle = -1
  logical :: rf_on = .true.		! By default, RF is assumed on
end type

! Common vars

type sodom2_com_struct
  type (lat_struct) :: lat
  type (bunch_struct) :: bunch
  integer :: ix_ele_eval
  integer :: ix_branch = 0                   ! Lattice branch being tracked.
  real(rp) :: time_start = 0
  logical :: debug = .false.
  integer :: mpi_rank = master_rank$
  logical :: using_mpi = .false.
  character(200) :: master_input_file = ''
  real(rp) :: Q_2pi(3)
  real(rp) :: ADST = 0
  real(rp) :: n_mat(6,6) = 0
  complex(rp), allocatable :: U11(:)
  complex(rp), allocatable :: U12(:)
  complex(rp), allocatable :: U21(:)
  complex(rp), allocatable :: U22(:)
  complex(rp), allocatable :: M(:,:)
  complex(rp), allocatable :: eig_vec(:,:)
  complex(rp), allocatable :: eig_val(:)
  integer :: idx_max_0 = 0		     ! Index of eigenvector with max |psi_0|
end type


contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_read_params(sodom2, sodom2_com)
type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
integer i, ix
character(200) arg
character(40) m_name
character(*), parameter :: r_name = 'sodom2_read_params'

namelist / params / sodom2 ! input file is fortran namelist format, can easily load into this struct

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

if (sodom2_com%master_input_file == '') sodom2_com%master_input_file = 'sodom2.init'

! Read parameters

if (.not. sodom2_com%using_mpi .or. sodom2_com%mpi_rank == master_rank$) then
  print '(2a)', 'Initialization file: ', trim(sodom2_com%master_input_file)
endif

open (1, file = sodom2_com%master_input_file, status = 'old', action = 'read')
read (1, nml = params)
close (1)

call bmad_parser (sodom2%lat_file, sodom2_com%lat)

print *,'Mode = ', sodom2%mode

! Check:
if (sodom2%J < 0) then
  call out_io (s_fatal$, r_name, 'sodom2%J MUST BE GREATER THAN ZERO. STOPPING HERE.')
  stop
endif

if (sodom2%mode < 0 .or. sodom2%mode > 3) then
  call out_io (s_fatal$, r_name, 'sodom2%mode MUST BE 1, 2, OR 3. STOPPING HERE.')
  stop
endif

if (sodom2%n_particle < 1) then
  call out_io (s_fatal$, r_name, 'sodom2%n_particle MUST BE GREATER THAN 1. STOPPING HERE.')
  stop
endif

end subroutine sodom2_read_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_init_params(sodom2, sodom2_com)

type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
type (lat_struct) :: lat
!type (branch_struct):: branch
type (ele_pointer_struct), allocatable :: eles(:)

integer n_loc, i
logical err

lat = sodom2_com%lat
bmad_com%auto_bookkeeper = .false.
bmad_com%spin_tracking_on = .true.

! Check if RF is on:
! print '(2a)', 'Note: RF is always set ON for SODOM-2'
! call set_on_off (rfcavity$, lat, on$)

if (sodom2%ele_eval /= '') then
  call lat_ele_locator (sodom2%ele_eval, lat, eles, n_loc, err)
  if (err .or. n_loc == 0) then
    print '(2a)', 'Evaluate element not found: ', trim(sodom2%ele_eval)
    stop
  endif
  if (n_loc > 1) then
    print '(2a)', 'Multiple elements found with ele_eval name: ', trim(sodom2%ele_eval)
    print '(a)', 'Will stop here.'
    stop
  endif
  sodom2_com%ix_branch = eles(1)%ele%ix_branch
  eles(1)%ele%type = '@START_ELE'
  sodom2_com%ix_ele_eval = eles(1)%ele%ix_ele
else
  sodom2_com%ix_branch = 0
  sodom2_com%ix_ele_eval = lat%ele(0)%ix_ele
endif

! set tracking method for each element to linear$
do i = 1, lat%branch(sodom2_com%ix_branch)%n_ele_track
  call set_ele_attribute(lat%branch(sodom2_com%ix_branch)%ele(i), "tracking_method = linear", err, .true., .true.)! , .false.)
enddo

call lattice_bookkeeper(lat)

call fullfilename (sodom2%output_file, sodom2%output_file)

! allocate SU2 quaternion discrete function
allocate(sodom2_com%U11(sodom2%n_particle))
allocate(sodom2_com%U12(sodom2%n_particle))
allocate(sodom2_com%U21(sodom2%n_particle))
allocate(sodom2_com%U22(sodom2%n_particle))

! allocate 2N x 2N matrix
allocate(sodom2_com%M(2*sodom2%n_particle, 2*sodom2%n_particle))

! allocate eig vecs/vals
allocate(sodom2_com%eig_vec(2*sodom2%n_particle, 2*sodom2%n_particle))
allocate(sodom2_com%eig_val(2*sodom2%n_particle))

end subroutine sodom2_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_init_bunch(sodom2, sodom2_com)

type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
type (bunch_struct) :: bunch
type (lat_struct) :: lat

real(rp) m(6,6), vec0(6), vec(6)
real(rp) :: jvec(6) = 0

real(rp) m_4(4,4)
complex(rp) :: eig_val4(4) = 0
complex(rp) :: eig_vec4(4,4) = 0

real(rp) phi, J1, J2, PI
real(rp) :: spin1(3) = (/1, 0, 0/)
real(rp) :: spin2(3) = (/0, 1, 0/)
real(rp) :: spin3(3) = (/0, 0, 1/)
complex(rp) :: eig_val(6) = 0
complex(rp) :: eig_vec(6,6) = 0
real(rp) :: Qr(6,6) = 0
real(rp) :: Qr2(2,2) = 0
real(rp) :: Qi(6,6) = 0
real(rp) :: Qi2(2,2) = 0
real(rp) :: c = 1.0_rp/sqrt(2.0_rp)
real(rp) :: mat_test(4,4) = 0

logical error, rf_on
integer i

PI = acos(-1.0_rp)
Qr2 = c*reshape((/ 1, 1, 0, 0 /), shape(Qr2))
Qr(1:2,1:2) = Qr2
Qr(3:4,3:4) = Qr2
Qr(5:6,5:6) = Qr2
Qi2 = c*reshape((/ 0, 0, 1, -1 /), shape(Qi2))
Qi(1:2,1:2) = Qi2
Qi(3:4,3:4) = Qi2
Qi(5:6,5:6) = Qi2


lat = sodom2_com%lat
! allocate memory for bunch
call reallocate_bunch(bunch, 3*sodom2%n_particle)

call transfer_matrix_calc (lat, m, ix1 = sodom2_com%ix_ele_eval, ix_branch = sodom2_com%ix_branch, one_turn = .true.) !, sodom2_com%ele_eval%ix_ele, sodom2_com%ele_eval%ix_branch, .true.)
rf_on = rf_is_on(lat%branch(sodom2_com%ix_branch))

if (rf_on == .true.) then
  print *, "RF is ON"
  call mat_eigen(m, eig_val, eig_vec, error, .false.)
  if (error) then
    print *, 'Cannot compute 1-turn matrix eigen vectors. Stopping here.'
    stop
  endif

  ! mat_eigen uses 1,3,5 for tunes:
  sodom2_com%Q_2pi(1) = real(log(eig_val(1))*(0.0_rp, 1.0_rp))
  sodom2_com%Q_2pi(2) = real(log(eig_val(3))*(0.0_rp, 1.0_rp))
  sodom2_com%Q_2pi(3) = real(log(eig_val(5))*(0.0_rp, 1.0_rp))

  eig_vec = transpose(eig_vec)
  sodom2_com%n_mat = matmul(real(eig_vec), Qr) - matmul(aimag(eig_vec), Qi)
 
  ! Apply rotation matrix to make N12, N34, N56 equal to 0:
  sodom2_com%n_mat = matmul(sodom2_com%n_mat, Rot3(MyTan(sodom2_com%n_mat(1,2), sodom2_com%n_mat(1,1)), MyTan(sodom2_com%n_mat(3,4), sodom2_com%n_mat(3,3)),  MyTan(sodom2_com%n_mat(5,6), sodom2_com%n_mat(5,5)) ))
else
  print *, "RF is OFF"
  if (sodom2%mode == 3) then
    print *, 'Cannot excite c-mode with RF off. Stopping here.'
    stop
  endif
  m_4 = m(1:4,1:4)
  call mat_eigen(m_4, eig_val4, eig_vec4, error, .false.)
  if (error) then
    print *, 'Cannot compute 1-turn matrix eigen vectors. Stopping here.'
    stop
  endif

  sodom2_com%Q_2pi(1) = real(log(eig_val4(1))*(0.0_rp, 1.0_rp))
  sodom2_com%Q_2pi(2) = real(log(eig_val4(3))*(0.0_rp, 1.0_rp))
  
  eig_vec4 = transpose(eig_vec4)
  sodom2_com%n_mat(1:4,1:4) = matmul(real(eig_vec4), Qr(1:4,1:4)) - matmul(aimag(eig_vec4), Qi(1:4,1:4))
  
  ! Apply rotation matrix to make N12, N34 equal to 0:
  sodom2_com%n_mat(1:4,1:4) = matmul(sodom2_com%n_mat(1:4,1:4), Rot2(MyTan(sodom2_com%n_mat(1,2), sodom2_com%n_mat(1,1)), MyTan(sodom2_com%n_mat(3,4), sodom2_com%n_mat(3,3)) ))
endif



do i = 1, sodom2%n_particle
  phi = (2d0 * PI * (i - 1)) / sodom2%n_particle
  J1 = sqrt(2*sodom2%J)*cos(phi)
  J2 = -sqrt(2*sodom2%J)*sin(phi)
  jvec(2*sodom2%mode-1) = J1
  jvec(2*sodom2%mode) = J2

  vec = matmul(sodom2_com%n_mat, jvec)
 

  ! initialize 3 particles w/ orthogonal spins
  call init_coord(bunch%particle(3*i-2), vec, sodom2_com%lat%ele(sodom2_com%ix_ele_eval), downstream_end$, sodom2_com%lat%param%particle, 1, spin = spin1)
  call init_coord(bunch%particle(3*i-1), vec, sodom2_com%lat%ele(sodom2_com%ix_ele_eval), downstream_end$, sodom2_com%lat%param%particle, 1, spin = spin2)
  call init_coord(bunch%particle(3*i), vec, sodom2_com%lat%ele(sodom2_com%ix_ele_eval), downstream_end$, sodom2_com%lat%param%particle, 1, spin = spin3)
enddo

sodom2_com%bunch = bunch

contains

function Rot3(a, b, c) result(mat)
  real(rp) a, b, c
  real(rp) mat(6,6)

  mat = 0.0d0
  mat(1:2,1:2) = Rot1(a)
  mat(3:4,3:4) = Rot1(b)
  mat(5:6,5:6) = Rot1(c)
end function

function Rot2(a, b) result(mat)
  real(rp) a, b
  real(rp) mat(4,4)

  mat = 0.0d0
  mat(1:2,1:2) = Rot1(a)
  mat(3:4,3:4) = Rot1(b)
end function

function Rot1(theta) result(mat)
  real(rp) theta
  real(rp) mat(2,2)

  mat(1,1) = cos(theta)
  mat(2,2) = mat(1,1)
  mat(1,2) = -sin(theta)
  mat(2,1) = -mat(1,2)
end function

end subroutine sodom2_init_bunch

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_track_bunch(sodom2, sodom2_com)

type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
logical err_flag
call track_bunch(sodom2_com%lat, sodom2_com%bunch, sodom2_com%lat%ele(sodom2_com%ix_ele_eval),sodom2_com%lat%ele(sodom2_com%ix_ele_eval), err_flag )

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

subroutine sodom2_construct_quaternions(sodom2, sodom2_com)
type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
complex(rp) U(2,2)

real(rp) R(3,3), q(0:3)
integer i

do i = 1, sodom2%n_particle
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

subroutine sodom2_construct_mat(sodom2, sodom2_com)
type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com

integer j, k, n, idx

logical err_flag

call fft_1d(sodom2_com%U11, -1)
call fft_1d(sodom2_com%U12, -1)
call fft_1d(sodom2_com%U21, -1)
call fft_1d(sodom2_com%U22, -1)

! center harmonics around 0:
sodom2_com%U11 = cshift(sodom2_com%U11, -(sodom2%n_particle-1)/2)/sodom2%n_particle
sodom2_com%U21 = cshift(sodom2_com%U21, -(sodom2%n_particle-1)/2)/sodom2%n_particle
sodom2_com%U12 = cshift(sodom2_com%U12, -(sodom2%n_particle-1)/2)/sodom2%n_particle
sodom2_com%U22 = cshift(sodom2_com%U22, -(sodom2%n_particle-1)/2)/sodom2%n_particle

do j = -(sodom2%n_particle - 1)/2, (sodom2%n_particle - 1)/2 
  do k = -(sodom2%n_particle - 1)/2, (sodom2%n_particle - 1)/2 
    
    ! enforce circular indexing of harmonics
    idx = modulo((j-k)+(sodom2%n_particle-1)/2, sodom2%n_particle)+1
    !M(2*j+n_particles:2*j+n_particles+1,2*k+n_particles:2*k+n_particles+1) = exp(-1i*2*pi*j*Q)*A(j-k);
    if (idx == 0) then
      print *, 'j = ', j, ', k = ', k
    endif
    sodom2_com%M(2*j+sodom2%n_particle,2*k+sodom2%n_particle) = exp(-(0,1)*j*sodom2_com%Q_2pi(sodom2%mode))*sodom2_com%U11(idx)
    sodom2_com%M(2*j+sodom2%n_particle,2*k+sodom2%n_particle+1) = exp(-(0,1)*j*sodom2_com%Q_2pi(sodom2%mode))*sodom2_com%U12(idx)
    sodom2_com%M(2*j+sodom2%n_particle+1,2*k+sodom2%n_particle) = exp(-(0,1)*j*sodom2_com%Q_2pi(sodom2%mode))*sodom2_com%U21(idx)
    sodom2_com%M(2*j+sodom2%n_particle+1,2*k+sodom2%n_particle+1) = exp(-(0,1)*j*sodom2_com%Q_2pi(sodom2%mode))*sodom2_com%U22(idx)
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

subroutine sodom2_calc_ADST(sodom2, sodom2_com)
type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
integer j, n_particle, max_idx
real(rp) :: PI
real(rp), parameter::inv_pi = 1.0_rp / acos(-1.0_rp)

PI = acos(-1.0_rp)

n_particle = sodom2%n_particle

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

subroutine sodom2_write(sodom2, sodom2_com)
type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
complex(rp) sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)
complex(rp) :: psi_n(2)
real(rp) :: n_axis(3)
integer iu, i
integer j
integer :: ia = 1
real(rp) phi, J1, J2
real(rp) :: jvec(6) = 0


sigma_1 = reshape((/ (0,0), (1,0), (1,0), (0,0) /), shape(sigma_1))
sigma_2 = reshape((/ (0,0), (0,1), (0,-1), (0,0) /), shape(sigma_2))
sigma_3 = reshape((/ (1,0), (0,0), (0,0), (-1,0) /), shape(sigma_3))

iu = lunget()

open (iu, file = sodom2%output_file, recl = 300)
write (iu,  '(3a)')     	'# lat_file                            	= ', sodom2%lat_file !quote(sodom2_com%lat)
write (iu, '(a, a)')   		'# ele_eval                   		= ', sodom2%ele_eval
write (iu,  '(a, i8)')  	'# n_particle                          	= ', sodom2%n_particle
write (iu,  '(a, i8)')  	'# mode                          	= ', sodom2%mode
write (iu,  '(a, 8es16.8)')  	'# J                          		= ', sodom2%J
write (iu,  '(a)') 	'# -----------------------------------------------------------'
write (iu,  '(a, 8es16.8)') 	'# ADST                        		= ', sodom2_com%ADST
write (iu,  '(a)') 	'# -----------------------------------------------------------'
write (iu, '(a)')  '## phi [rad]	 |          x              px               y              py               z              pz        |         nx              ny              nz'

! if no angles_to_eval given, just write the ISF for each particle used in tracking:
if (sodom2%angles_to_eval(1) == -10) then
  do i = 1, sodom2%n_particle
    jvec = 0.0
    phi = (2d0 * PI * (i - 1)) / sodom2%n_particle
    J1 = sqrt(2*sodom2%J)*cos(phi)
    J2 = -sqrt(2*sodom2%J)*sin(phi)
    jvec(2*sodom2%mode-1) = J1
    jvec(2*sodom2%mode) = J2

    jvec = matmul(sodom2_com%n_mat, jvec)
    psi_n = (/ (0,0), (0,0) /)
    do j = -(sodom2%n_particle-1)/2, (sodom2%n_particle-1)/2
      psi_n = psi_n + sodom2_com%eig_vec(2*j+sodom2%n_particle:2*j+sodom2%n_particle+1, sodom2_com%idx_max_0)*exp((0,1)*j*phi);
    enddo
  
    n_axis(1) = dot_product(psi_n, matmul(sigma_1, psi_n))
    n_axis(2) = dot_product(psi_n, matmul(sigma_2, psi_n))
    n_axis(3) = dot_product(psi_n, matmul(sigma_3, psi_n))
    n_axis(1:3) = n_axis(1:3)/norm2(n_axis(1:3))


    write (iu, '(1es16.8, 3x, 6es16.8, 3x, 3es16.8)') phi, jvec, n_axis(1), n_axis(2), n_axis(3)
    ia = ia + 1
 enddo
else
  do while(sodom2%angles_to_eval(ia) /= -10 .and. ia <= angle_max$)
    jvec = 0.0
    phi = sodom2%angles_to_eval(ia)
    J1 = sqrt(2*sodom2%J)*cos(phi)
    J2 = -sqrt(2*sodom2%J)*sin(phi)
    jvec(2*sodom2%mode-1) = J1
    jvec(2*sodom2%mode) = J2

    jvec = matmul(sodom2_com%n_mat, jvec)

    psi_n = (/ (0,0), (0,0) /)
    do j = -(sodom2%n_particle-1)/2, (sodom2%n_particle-1)/2
      psi_n = psi_n + sodom2_com%eig_vec(2*j+sodom2%n_particle:2*j+sodom2%n_particle+1, sodom2_com%idx_max_0)*exp((0,1)*j*phi);
    enddo
  
    n_axis(1) = dot_product(psi_n, matmul(sigma_1, psi_n))
    n_axis(2) = dot_product(psi_n, matmul(sigma_2, psi_n))
    n_axis(3) = dot_product(psi_n, matmul(sigma_3, psi_n))
    n_axis(1:3) = n_axis(1:3)/norm2(n_axis(1:3))


    write (iu, '(1es16.8, 3x, 6es16.8, 3x, 3es16.8)') phi, jvec, n_axis(1), n_axis(2), n_axis(3)
    ia = ia + 1
  enddo
endif



close(iu)

end subroutine sodom2_write

end module
