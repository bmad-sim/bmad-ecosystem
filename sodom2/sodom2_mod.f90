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
use, intrinsic :: iso_c_binding

implicit none
include 'fftw3.f03'

integer, parameter :: master_rank$  = 0
integer, parameter :: a_mode$ = 1
integer, parameter :: b_mode$ = 2
integer, parameter :: c_mode$ = 3
real(rp), parameter :: J_tol$ = 1e-19
!integer, parameter :: angle_max$ = 1000
integer, private :: i_loop

! User settable parameters

type sodom2_params_struct
  character(40) :: ele_eval = ''
  character(200) :: lat_file = ''
  character(200) :: particle_output_file = 'sodom2.out'
  character(200) :: n_axis_output_file = 'n_axis.out' 
  real(rp) :: J(3) = [-1, -1, -1] 
  integer :: n_samples(3) = [1, 1, 1]
  logical :: linear_tracking = .true.
  logical :: add_closed_orbit_to_particle_output = .false.
  logical :: write_as_beam_init = .true.
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
  real(rp) :: Q_2pi(3) = 0
  real(rp) :: ADST = 0
  real(rp) :: n_mat(6,6) = 0
  integer :: n(3) = [1, 1, 1]			! n_samples in each plane
  complex(rp), allocatable :: U11(:,:,:)
  complex(rp), allocatable :: U12(:,:,:)
  complex(rp), allocatable :: U21(:,:,:)
  complex(rp), allocatable :: U22(:,:,:)
  complex(rp), allocatable :: M(:,:)
  complex(rp), allocatable :: eig_vec(:,:)
  complex(rp), allocatable :: eig_val(:)
  type (coord_struct), allocatable :: closed_orb(:)
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

! Check:
if (sodom2%J(1) < 0 .and. sodom2%J(2) < 0 .and. sodom2%J(3) < 0) then
  call out_io (s_fatal$, r_name, 'ATLEAST ONE sodom2%J MUST BE SPECIFIED. STOPPING HERE.')
  stop
endif

if (sodom2%J(1) < 0) then
  print *, 'sodom2%J(1) not specified. Assuming sodom2%J(1) = 0'
  sodom2%J(1) = 0
endif

if (sodom2%J(2) < 0) then
  print *, 'sodom2%J(2) not specified. Assuming sodom2%J(2) = 0'
  sodom2%J(2) = 0
endif

if (sodom2%J(3) < 0) then
  print *, 'sodom2%J(3) not specified. Assuming sodom2%J(3) = 0'
  sodom2%J(3) = 0
endif

do i = 1,3
  if (sodom2%n_samples(i) < 1) then
    call out_io (s_fatal$, r_name, 'sodom2%n_samples MUST BE GREATER THAN 1. STOPPING HERE.')
    stop
  endif

  if (mod(sodom2%n_samples(i),2) == 0) then
    sodom2%n_samples(i) = sodom2%n_samples(i) + 1
    print '(a,I0,a,I0)', ' Setting n_samples(',i,') = ', sodom2%n_samples(i)
  endif
enddo

end subroutine sodom2_read_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_init_params(sodom2, sodom2_com)

type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
type (lat_struct) :: lat
type (ele_pointer_struct), allocatable :: eles(:)
integer N

integer n_loc, i
logical err
logical :: J(3) = .false.

lat = sodom2_com%lat
bmad_com%auto_bookkeeper = .false.
bmad_com%spin_tracking_on = .true.

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

if (sodom2%linear_tracking) then
! set tracking method for each element to linear$
  do i = 1, lat%branch(sodom2_com%ix_branch)%n_ele_track
    call set_ele_attribute(lat%branch(sodom2_com%ix_branch)%ele(i), "tracking_method = linear", err, .true., .true.)! , .false.)
  enddo
!else
! set tracking method to Bmad standard
!  do i = 1, lat%branch(sodom2_com%ix_branch)%n_ele_track
!    call set_ele_attribute(lat%branch(sodom2_com%ix_branch)%ele(i), "tracking_method = bmad_standard", err, .true., .true.)! , .false.)
!  enddo
endif

call lattice_bookkeeper(lat)

call fullfilename (sodom2%particle_output_file, sodom2%particle_output_file)
call fullfilename (sodom2%n_axis_output_file, sodom2%n_axis_output_file)

! Determine number of particles necessary + memory allocation:
do i =1,3
  if (sodom2%J(i) > J_tol$) then
    J(i) = .true.
    sodom2_com%n(i) = sodom2%n_samples(i)
  endif
enddo

print '(a, 3l4)', ' Oscillating modes: ', [J(1), J(2), J(3)]

N = sodom2_com%n(1)*sodom2_com%n(2)*sodom2_com%n(3)

! allocate SU2 quaternion discrete function
allocate(sodom2_com%U11(sodom2_com%n(1),sodom2_com%n(2),sodom2_com%n(3)))
allocate(sodom2_com%U12(sodom2_com%n(1),sodom2_com%n(2),sodom2_com%n(3)))
allocate(sodom2_com%U21(sodom2_com%n(1),sodom2_com%n(2),sodom2_com%n(3)))
allocate(sodom2_com%U22(sodom2_com%n(1),sodom2_com%n(2),sodom2_com%n(3)))

! allocate 2N x 2N matrix
allocate(sodom2_com%M(2*N, 2*N))

! allocate eig vecs/vals
allocate(sodom2_com%eig_vec(2*N,2*N))
allocate(sodom2_com%eig_val(2*N))

end subroutine sodom2_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_init_bunch(sodom2, sodom2_com)

type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
type (bunch_struct) :: bunch
type (lat_struct) :: lat
type (coord_struct) :: closed_orbit

real(rp) m(6,6), vec0(6), vec(6)
real(rp) :: jvec(6) = 0

real(rp) m_4(4,4)
complex(rp) :: eig_val4(4) = 0
complex(rp) :: eig_vec4(4,4) = 0

real(rp) PI, J1_a, J1_b, J2_a, J2_b, J3_a, J3_b, phi1, phi2, phi3
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
integer i, j, k, N

PI = acos(-1.0_rp)
Qr2 = c*reshape((/ 1, 1, 0, 0 /), shape(Qr2))
Qr(1:2,1:2) = Qr2
Qr(3:4,3:4) = Qr2
Qr(5:6,5:6) = Qr2
Qi2 = c*reshape((/ 0, 0, 1, -1 /), shape(Qi2))
Qi(1:2,1:2) = Qi2
Qi(3:4,3:4) = Qi2
Qi(5:6,5:6) = Qi2

N = sodom2_com%n(1)*sodom2_com%n(2)*sodom2_com%n(3)
lat = sodom2_com%lat
! allocate memory for bunch
call reallocate_bunch(bunch, 3*N)


call transfer_matrix_calc (lat, m, ix1 = sodom2_com%ix_ele_eval, ix_branch = sodom2_com%ix_branch, one_turn = .true.) !, sodom2_com%ele_eval%ix_ele, sodom2_com%ele_eval%ix_branch, .true.)
rf_on = rf_is_on(lat%branch(sodom2_com%ix_branch))

if (rf_on) then
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
  if (sodom2_com%n(3) /= 1) then
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

! Calc closed orbit:
call closed_orbit_calc(lat, sodom2_com%closed_orb, ix_branch = sodom2_com%ix_branch, print_err = .true.)
!call twiss_and_track (lat, sodom2_com%bmad_closed_orb, ix_branch = sodom2_com%ix_branch)
!print '(a, 6es16.8)', 'Closed orbit: ', sodom2_com%bmad_closed_orb(sodom2_com%ix_ele_eval)%vec

do i = 1, sodom2_com%n(1)
  phi1 = (2.0_rp * PI * (i - 1)) / sodom2_com%n(1)
  J1_a = sqrt(2*sodom2%J(1))*cos(phi1)
  J1_b = -sqrt(2*sodom2%J(1))*sin(phi1)
  jvec(1) = J1_a
  jvec(2) = J1_b

  do j =1, sodom2_com%n(2)
    phi2 = (2.0_rp * PI * (j - 1)) / sodom2_com%n(2)
    J2_a = sqrt(2*sodom2%J(2))*cos(phi2)
    J2_b = -sqrt(2*sodom2%J(2))*sin(phi2)
    jvec(3) = J2_a
    jvec(4) = J2_b

    do k =1, sodom2_com%n(3)
      phi3 = (2.0_rp * PI * (k - 1)) / sodom2_com%n(3)
      J3_a = sqrt(2*sodom2%J(3))*cos(phi3)
      J3_b = -sqrt(2*sodom2%J(3))*sin(phi3)
      jvec(5) = J3_a
      jvec(6) = J3_b

      vec = matmul(sodom2_com%n_mat, jvec)
      vec = vec + sodom2_com%closed_orb(sodom2_com%ix_ele_eval)%vec
 
      ! initialize 3 particles w/ orthogonal spins
      call init_coord(bunch%particle(3*(k + (j-1)*sodom2_com%n(3) + (i-1)*sodom2_com%n(2)*sodom2_com%n(3))-2), vec, sodom2_com%lat%ele(sodom2_com%ix_ele_eval), downstream_end$, sodom2_com%lat%param%particle, 1, spin = spin1)
      call init_coord(bunch%particle(3*(k + (j-1)*sodom2_com%n(3) + (i-1)*sodom2_com%n(2)*sodom2_com%n(3))-1), vec, sodom2_com%lat%ele(sodom2_com%ix_ele_eval), downstream_end$, sodom2_com%lat%param%particle, 1, spin = spin2)
      call init_coord(bunch%particle(3*(k + (j-1)*sodom2_com%n(3) + (i-1)*sodom2_com%n(2)*sodom2_com%n(3))), vec, sodom2_com%lat%ele(sodom2_com%ix_ele_eval), downstream_end$, sodom2_com%lat%param%particle, 1, spin = spin3)
    enddo
  enddo
enddo
!call reallocate_coord(closed_orbit,0)
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
integer i,j,k

do i = 1, sodom2_com%n(1)
  do j =1, sodom2_com%n(2)
    do k =1, sodom2_com%n(3)
      R(:,1) = sodom2_com%bunch%particle(3*(k + (j-1)*sodom2_com%n(3) + (i-1)*sodom2_com%n(2)*sodom2_com%n(3))-2)%spin
      R(:,2) = sodom2_com%bunch%particle(3*(k + (j-1)*sodom2_com%n(3) + (i-1)*sodom2_com%n(2)*sodom2_com%n(3))-1)%spin
      R(:,3) = sodom2_com%bunch%particle(3*(k + (j-1)*sodom2_com%n(3) + (i-1)*sodom2_com%n(2)*sodom2_com%n(3)))%spin

      q = w_mat_to_quat(R)
      call quat_to_SU2(q,U)
      sodom2_com%U11(i,j,k) = U(1,1)
      sodom2_com%U12(i,j,k) = U(1,2)
      sodom2_com%U21(i,j,k) = U(2,1)
      sodom2_com%U22(i,j,k) = U(2,2)
    enddo
  enddo
enddo


end subroutine sodom2_construct_quaternions

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_fft(arr)
integer(8) plan
complex(rp) arr(:,:,:)
call dfftw_plan_dft_3d(plan,size(arr,1), size(arr,2), size(arr,3), arr, arr, FFTW_FORWARD, FFTW_ESTIMATE)
call dfftw_execute_dft(plan, arr, arr)
call dfftw_destroy_plan(plan)

end subroutine sodom2_fft

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_construct_mat(sodom2, sodom2_com)
type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com

integer i, j1, j2, j3, k1, k2, k3, idx,  row_idx, col_idx, N, n1, n2, n3, idx1, idx2, idx3
logical err_flag
real(rp) jdotQ_2pi

n1 = sodom2_com%n(1)
n2 = sodom2_com%n(2)
n3 = sodom2_com%n(3)

call sodom2_fft(sodom2_com%U11)
call sodom2_fft(sodom2_com%U12)
call sodom2_fft(sodom2_com%U21)
call sodom2_fft(sodom2_com%U22)

! center harmonics around 0:
do i = 1, 3
  sodom2_com%U11 = cshift(sodom2_com%U11, -(sodom2_com%n(i)-1)/2, i)
  sodom2_com%U21 = cshift(sodom2_com%U21, -(sodom2_com%n(i)-1)/2, i)
  sodom2_com%U12 = cshift(sodom2_com%U12, -(sodom2_com%n(i)-1)/2, i)
  sodom2_com%U22 = cshift(sodom2_com%U22, -(sodom2_com%n(i)-1)/2, i)
enddo

N = n1*n2*n3

sodom2_com%U11 = sodom2_com%U11/N
sodom2_com%U21 = sodom2_com%U21/N
sodom2_com%U12 = sodom2_com%U12/N
sodom2_com%U22 = sodom2_com%U22/N


do j1 = -(n1 - 1)/2, (n1 - 1)/2 
  do j2 = -(n2 - 1)/2, (n2 - 1)/2
    do j3 = -(n3 - 1)/2, (n3 - 1)/2
      do k1 = -(n1 - 1)/2, (n1 - 1)/2
        do k2 = -(n2 - 1)/2, (n2 - 1)/2
          do k3 = -(n3 - 1)/2, (n3 - 1)/2

            ! enforce circular indexing of harmonics
            idx1 = modulo((j1-k1)+(n1-1)/2, n1)+1
	    idx2 = modulo((j2-k2)+(n2-1)/2, n2)+1
	    idx3 = modulo((j3-k3)+(n3-1)/2, n3)+1
   	    
 	    ! Each row_idx, col_idx correspond to a 2x2 SU(2) matrix
	    row_idx = (j3 + (n3-1)/2 + 1) + (j2 + (n2-1)/2)*n3 + (j1 + (n1-1)/2)*n3*n2 
	    col_idx = (k3 + (n3-1)/2 + 1) + (k2 + (n2-1)/2)*n3 + (k1 + (n1-1)/2)*n3*n2 

	    jdotQ_2pi = j1*sodom2_com%Q_2pi(1)+j2*sodom2_com%Q_2pi(2)+j3*sodom2_com%Q_2pi(3)

	    sodom2_com%M(2*row_idx - 1, 2*col_idx - 1) = exp(-(0,1)*jdotQ_2pi)*sodom2_com%U11(idx1,idx2,idx3)
	    sodom2_com%M(2*row_idx - 1, 2*col_idx    ) = exp(-(0,1)*jdotQ_2pi)*sodom2_com%U12(idx1,idx2,idx3)
	    sodom2_com%M(2*row_idx    , 2*col_idx - 1) = exp(-(0,1)*jdotQ_2pi)*sodom2_com%U21(idx1,idx2,idx3)
	    sodom2_com%M(2*row_idx    , 2*col_idx    ) = exp(-(0,1)*jdotQ_2pi)*sodom2_com%U22(idx1,idx2,idx3)
          enddo
        enddo
      enddo
    enddo
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
call reallocate_coord(sodom2_com%closed_orb, 0)

end subroutine sodom2_deallocate_memory

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_eig(sodom2, sodom2_com)
type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com

call la_geev(sodom2_com%M, sodom2_com%eig_val, VL = sodom2_com%eig_vec) !, err, .true.)
call sodom2_determine_ADST(sodom2_com)

end subroutine sodom2_eig

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_determine_ADST(sodom2_com)
type (sodom2_com_struct), target :: sodom2_com
integer j, max_idx, N
real(rp) PI
real(rp), parameter::inv_pi = 1.0_rp / acos(-1.0_rp)

N = sodom2_com%n(1)*sodom2_com%n(2)*sodom2_com%n(3)

PI = acos(-1.0_rp)

! determine eigenvector with largest |psi_0|
max_idx = 1
do j = 1, 2*N
    if (real(sodom2_com%eig_vec(N,j)*conjg(sodom2_com%eig_vec(N,j))) > real(sodom2_com%eig_vec(N, max_idx)*conjg(sodom2_com%eig_vec(N, max_idx)))) then
        max_idx = j
    endif
enddo

sodom2_com%idx_max_0 = max_idx
sodom2_com%ADST = real(log(sodom2_com%eig_val(max_idx))*(0.0_rp, inv_pi))

end subroutine sodom2_determine_ADST

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_write_n(sodom2, sodom2_com)
type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
complex(rp) :: psi_n(2)
real(rp) :: n_axis(3)
real(rp) PI
integer iu, j1, j2, j3, n1, n2, n3, row_idx

n1 = sodom2_com%n(1)
n2 = sodom2_com%n(2)
n3 = sodom2_com%n(3)

PI = acos(-1.0_rp)

iu = lunget()

open (iu, file = sodom2%n_axis_output_file, recl = 300)
write (iu,  '(3a)')     	'# lat_file                            	= ', sodom2%lat_file !quote(sodom2_com%lat)
write (iu, '(a, a)')   		'# ele_eval                   		= ', sodom2%ele_eval
write (iu,  '(a, 3i8)')  	'# n_samples                         	= ', sodom2%n_samples
write (iu,  '(a, 8es16.8)')  	'# J                          		= ', sodom2%J
write (iu,  '(a)') 	'# -----------------------------------------------------------'
write (iu,  '(a, 8es16.8)') 	'# Orbital Tunes                  	= ', sodom2_com%Q_2pi/(2.0_rp*PI)
write (iu,  '(a)') 	'# -----------------------------------------------------------'
write (iu,  '(a, 8es16.8)') 	'# ADST                        		= ', sodom2_com%ADST
write (iu,  '(a)') 	'# -----------------------------------------------------------'
write (iu, '(a)')  '##    j1      j2      j3  |      Psi_n1                          Psi_n2'
write (iu, '(a)')  '##  		          |         Re             Im               Re             Im'
do j1 = -(n1 - 1)/2, (n1 - 1)/2 
  do j2 = -(n2 - 1)/2, (n2 - 1)/2
    do j3 = -(n3 - 1)/2, (n3 - 1)/2
      row_idx = (j3 + (n3-1)/2 + 1) + (j2 + (n2-1)/2)*n3 + (j1 + (n1-1)/2)*n3*n2 
      psi_n = sodom2_com%eig_vec(2*row_idx-1:2*row_idx, sodom2_com%idx_max_0)
      write (iu, '(3i8, 3x, 4es16.8)') j1, j2, j3, real(psi_n(1)), aimag(psi_n(1)), real(psi_n(2)), aimag(psi_n(2))
    enddo
  enddo
enddo

close(iu)

end subroutine sodom2_write_n

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine sodom2_write_particles(sodom2, sodom2_com)
type (sodom2_params_struct), target :: sodom2
type (sodom2_com_struct), target :: sodom2_com
complex(rp) sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)
complex(rp) :: psi_n(2)
real(rp) :: n_axis(3)
integer iu, i1, i2, i3, j1, j2, j3, N, n1, n2, n3, row_idx
integer :: ia = 1
real(rp) J1_a, J1_b, J2_a, J2_b, J3_a, J3_b, phi1, phi2, phi3, jdotphi
real(rp) :: jvec(6) = 0
real(rp) :: vec(6) = 0


sigma_1 = reshape((/ (0,0), (1,0), (1,0), (0,0) /), shape(sigma_1))
sigma_2 = reshape((/ (0,0), (0,1), (0,-1), (0,0) /), shape(sigma_2))
sigma_3 = reshape((/ (1,0), (0,0), (0,0), (-1,0) /), shape(sigma_3))
n1 = sodom2_com%n(1)
n2 = sodom2_com%n(2)
n3 = sodom2_com%n(3)

iu = lunget()
open (iu, file = sodom2%particle_output_file, recl = 300)
if (sodom2%write_as_beam_init) then
  write (iu,  '(i0)') sodom2_com%ix_ele_eval
  write (iu,  '(i0)') 1
  write (iu,  '(i0)') n1*n2*n3
  write (iu,  '(a)') 'BEGIN_BUNCH'
  write (iu,  '(2x, a)')  species_name(sodom2_com%lat%param%particle)
  write (iu,  '(2x, a)')  '0.0'
  write (iu,  '(2x, a)')  '0.0'
  write (iu,  '(2x, a)')  '0.0'
else
  write (iu,  '(3a)')     	'# lat_file                            	= ', sodom2%lat_file !quote(sodom2_com%lat)
  write (iu, '(a, a)')   		'# ele_eval                   		= ', sodom2%ele_eval
  write (iu,  '(a, 3i8)')  	'# n_samples                         	= ', sodom2%n_samples
  write (iu,  '(a, 8es16.8)')  	'# J                          		= ', sodom2%J
  write (iu,  '(a)') 	'# -----------------------------------------------------------'
  write (iu,  '(a, 8es16.8)') 	'# ADST                        		= ', sodom2_com%ADST
  write (iu,  '(a)') 	'# -----------------------------------------------------------'
  write (iu, '(a)')  '## phi_1          phi_2          phi_3          |          x              px               y              py               z              pz        |         nx              ny              nz'
endif

do i1 = 1, n1
  phi1 = (2.0_rp * PI * (i1 - 1)) / n1
  J1_a = sqrt(2*sodom2%J(1))*cos(phi1)
  J1_b = -sqrt(2*sodom2%J(1))*sin(phi1)
  jvec(1) = J1_a
  jvec(2) = J1_b

  do i2 =1, n2
    phi2 = (2.0_rp * PI * (i2 - 1)) / n2
    J2_a = sqrt(2*sodom2%J(2))*cos(phi2)
    J2_b = -sqrt(2*sodom2%J(2))*sin(phi2)
    jvec(3) = J2_a
    jvec(4) = J2_b

    do i3 =1, n3
      phi3 = (2.0_rp * PI * (i3 - 1)) / n3
      J3_a = sqrt(2*sodom2%J(3))*cos(phi3)
      J3_b = -sqrt(2*sodom2%J(3))*sin(phi3)
      jvec(5) = J3_a
      jvec(6) = J3_b

      vec = matmul(sodom2_com%n_mat, jvec)
      if (sodom2%add_closed_orbit_to_particle_output) then
        vec = vec + sodom2_com%closed_orb(sodom2_com%ix_ele_eval)%vec
      endif

      psi_n = (/ (0,0), (0,0) /)
      do j1 = -(n1-1)/2, (n1-1)/2
        do j2 = -(n2-1)/2, (n2-1)/2
          do j3 = -(n3-1)/2, (n3-1)/2
	    jdotphi = phi1*j1+phi2*j2+phi3*j3
	    row_idx = (j3 + (n3-1)/2 + 1) + (j2 + (n2-1)/2)*n3 + (j1 + (n1-1)/2)*n3*n2 
            psi_n = psi_n + sodom2_com%eig_vec(2*row_idx-1:2*row_idx, sodom2_com%idx_max_0)*exp((0,1)*jdotphi);
          enddo
        enddo
      enddo
  
      n_axis(1) = dot_product(psi_n, matmul(sigma_1, psi_n))
      n_axis(2) = dot_product(psi_n, matmul(sigma_2, psi_n))
      n_axis(3) = dot_product(psi_n, matmul(sigma_3, psi_n))
      n_axis(1:3) = n_axis(1:3)/norm2(n_axis(1:3))

      if (sodom2%write_as_beam_init) then
        write (iu, '(6es16.8, 3x, i0, 3x, i0, 3x, 3es16.8)') vec, 0, 1, n_axis(1), n_axis(2), n_axis(3)
      else
        write (iu, '(3es16.8, 3x, 6es16.8, 3x, 3es16.8)') phi1, phi2, phi3, vec, n_axis(1), n_axis(2), n_axis(3)
      endif
    enddo
  enddo
enddo
if (sodom2%write_as_beam_init) then
  write (iu, '(a)') 'END_BUNCH'
endif
close(iu)

end subroutine sodom2_write_particles

end module
