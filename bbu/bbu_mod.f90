module bbu_mod

use precision_def
use bmad

type bbu_info
  real(rp), allocatable :: time_c(:)  ! Arrival time of each cavity.
  real(rp), allocatable :: mat_c(:)   ! Transfer matrices between cavities.
  real(rp), allocatable :: freq_c(:)  ! frequency of each hom.
  real(rp), allocatable :: Q_c(:)     ! Q value of each hom.
  real(rp), allocatable :: RoQ_c(:)
  real(rp), allocatable :: angle_c(:)
  real(rp), allocatable :: power_c(:)
  integer,  allocatable :: hom(:)     ! Number of homs in each cavity
  real(rp)  current                   ! Beam current
  real(rp)  b_freq
  integer   b_time
  integer   p_time
  real(rp)  n_amp
  real(rp)  p_amp
  real(rp)  n_on
  integer   cavity                    ! Number of cavities
  integer   hom_num
end type bbu_info

type pair_num
  real(rp), allocatable ::  var(:)
  real(rp), allocatable ::  threshold(:)
end type pair_num

type hom_power
  real(rp)  power
  integer   c_ind
  integer   h_ind
end type hom_power

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine get_info(lattice, cur, bbu, Hpower, output)
!
! Extract information from the lattice for bbu simulation
!
! Modules Needed:
!   use bbu_mod
!
! Input:
!   lattice -- lat_struct: lattice file of the accelerator
!   cur     -- real: initial trial threshold current
!   bbu     -- bbu_info: will be used to store all the information needed for the bbu simulation
!   Hpower  -- HOM_power: containing the power of each HOM mode at the end of tracking
!   Output  -- logical: if true write the lattice information into files.
!
! Output
!   bbu -- bbu_info: storing all the information extracted from lattice
!-  

subroutine get_info(lattice, cur, bbu, Hpower, output)

use random_mod

implicit none

type (bbu_info) :: bbu
type (lat_struct), target  :: lattice
type (hom_power), allocatable :: Hpower(:)
type (ele_struct), pointer :: ele

real(rp), intent(in) :: cur
logical, intent(in) :: output

! Arrays used to store lattice information

real(rp) :: mat(6,6), oldmat(6,6), testmat(6,6), testmat1(6,6)
real(rp), allocatable :: erlmat(:,:,:)
real(rp), allocatable :: erltime(:)
integer,  allocatable :: cavityind(:)

real(rp)  time, rr, Vz, P0i, P0f, gamma, cavity_i
logical judge

integer i, j, l, n, u, k, ic1, ic2, homtotal
integer c_num, h_num
integer :: m_size = 4

real(rp), parameter :: noise_amp = 0.0
real(rp), parameter :: pos_amp = 0.0
real(rp), parameter :: noise_on = 0.0
real(rp), parameter :: pw = 0.1

! Compute the 6 by 6 transfer matrices 

call twiss_propagate_all(lattice)
call lat_make_mat6(lattice, -1)

! Print out the beta function

if (output) then
  open(20, file='mat', status='unknown', action='write')
endif

! Count the number of lcavities with homs in the lattice

c_num=0
if (allocated(bbu%hom)) deallocate (bbu%mat_c, bbu%time_c, bbu%Q_c, &
          bbu%hom, bbu%RoQ_c, bbu%freq_c, bbu%angle_c, bbu%power_c, Hpower)
allocate(cavityind(lattice%n_ele_track))

do i=0, lattice%n_ele_track

  ele => lattice%ele(i) 

  if (ele%key /= lcavity$) cycle
  if (.not. associated(ele%wake)) cycle
  if (all (ele%wake%lr(:)%R_over_Q < 1E-10)) cycle
  c_num = c_num + 1
  cavityind(c_num) = i

enddo

! Count the number of HOMs in each cavity

allocate(bbu%hom(c_num/2))

judge=.false.    !  logical -- true if the rf cavity has wake fields
homtotal = 0     !  Total number of HOMs
k = 1

do i= 1, c_num/2

  ele => lattice%ele(cavityind(i))

  h_num=0
  do j=1, size(ele%wake%lr)
    if (ele%wake%lr(j)%R_over_Q < 1E-10) cycle
    if (ele%wake%lr(j)%polarized) then
      h_num = h_num + 1
      homtotal = homtotal + 1
    else 
      h_num = h_num + 2
      homtotal = homtotal + 2
    endif
  enddo

  bbu%hom(k) = h_num

enddo

! Allocate the dimesion of the transport matrix arrays and time arrays

allocate(erlmat(c_num, m_size, m_size))
allocate(erltime(c_num))

allocate(bbu%mat_c((c_num-1)*m_size*m_size))
allocate(bbu%time_c(c_num))
allocate(bbu%Q_c(homtotal))
allocate(bbu%RoQ_c(homtotal))
allocate(bbu%freq_c(homtotal))
allocate(bbu%angle_c(homtotal))
allocate(bbu%power_c(homtotal))
allocate(Hpower(homtotal))

! Find the time between each cavity in the low energy lattice
! The time between two cavities is stored in erltime(k)

time = 0
k = 1     

do i = 1, lattice%n_ele_track
   
  ele => lattice%ele(i) 
  Vz = c_light * ele%value(p0c$) / ele%value(e_tot$)

  if (ele%key == LCAVITY$ ) then
    time=time+(ele%value(l$))/c_light*(1+0.5/(gamma*lattice%ele(i-1)%value(E_TOT$)/m_electron))
  else
    time = time + ele%value(l$) / Vz
  endif

  if (i == cavityind(k)) then
    erltime(k)=time
    time = 0
    k = k + 1
  endif
   
enddo

! Calculate Transport Matrices for the low energy lattice
! Matrix elements are stored in erlmat(i,j,k)

call mat_make_unit (oldmat)
call mat_make_unit (testmat)
P0i = lattice%ele(0)%value(p0c$)     ! Get the first longitudinal reference momentum
k = 1

do i = 0, lattice%n_ele_track
   
  ele => lattice%ele(i) 
  mat = matmul(ele%mat6, oldmat) 
  
  if (i /= cavityind(k)) cycle

  P0f=ele%value(p0c$)

  mat(1,2)=mat(1,2)/P0i
  mat(1,4)=mat(1,4)/P0i
  mat(2,1)=mat(2,1)*P0f
  mat(2,2)=mat(2,2)*P0f/P0i
  mat(2,3)=mat(2,3)*P0f
  mat(2,4)=mat(2,4)*P0f/P0i
  mat(3,2)=mat(3,2)/P0i
  mat(3,4)=mat(3,4)/P0i
  mat(4,1)=mat(4,1)*P0f
  mat(4,2)=mat(4,2)*P0f/P0i
  mat(4,3)=mat(4,3)*P0f
  mat(4,4)=mat(4,4)*P0f/P0i

  erlmat(k, 1:m_size, 1:m_size) = mat(1:m_size, 1:m_size)

  call mat_make_unit (mat)
  k = k + 1
  P0i = P0f 
  
enddo

! Read in HOM parameters f, Q, R/Q, theta
! Giving random initial power to each HOM

u = 1
do k = 1, c_num/2
  i = cavityind(k)
  ic1 = ele%ic1_lord
  ic2 = lattice%ic(ic1)
  h_num = 1

  do l=1, size(ele%wake%lr)
    if (ele%wake%lr(l)%R_over_Q < 1E-10) cycle
      
    call ran_gauss(rr)
    if (ele%wake%lr(l)%polarized) then
      bbu%Q_c(u)=ele%wake%lr(l)%Q
      bbu%RoQ_c(u)=2*(c_light/(2*pi*ele%wake%lr(l)%freq))**2*ele%wake%lr(l)%R_over_Q
      bbu%freq_c(u)=ele%wake%lr(l)%freq
      bbu%angle_c(u)=ele%wake%lr(l)%angle*2*pi
      bbu%power_c(u)=pw*abs(rr)
      Hpower(u)%c_ind = lattice%control(ic2)%ix_lord
      Hpower(u)%h_ind = h_num
      u = u + 1
      h_num = h_num + 1
    else
      bbu%Q_c(u)=ele%wake%lr(l)%Q
      bbu%RoQ_c(u)=2*(c_light/(2*pi*ele%wake%lr(l)%freq))**2*ele%wake%lr(l)%R_over_Q
      bbu%freq_c(u)=ele%wake%lr(l)%freq
      bbu%angle_c(u)=0
      bbu%power_c(u)=pw*abs(rr)
      Hpower(u)%c_ind = lattice%control(ic2)%ix_lord
      Hpower(u)%h_ind = h_num

      bbu%Q_c(u+1)=ele%wake%lr(l)%Q
      bbu%RoQ_c(u+1)=2*(c_light/(2*pi*ele%wake%lr(l)%freq))**2*ele%wake%lr(l)%R_over_Q
      bbu%freq_c(u+1)=ele%wake%lr(l)%freq
      bbu%angle_c(u+1)=0.5*pi
      bbu%power_c(u+1)=pw*abs(rr)
      Hpower(u+1)%c_ind = lattice%control(ic2)%ix_lord
      Hpower(u+1)%h_ind = h_num+1 

      u=u+2
      h_num = h_num + 2
    endif
  enddo
enddo

bbu%cavity = c_num/2 ! For double pass ERL

write(*,*)'The total number of rf cavities with HOMs is ', bbu%cavity

! Write bbu information into a file

bbu%time_c(1) = 0
do i=2, c_num
  bbu%time_c(i) = erltime(i)+bbu%time_c(i-1)
  write(20,'(ES20.12)') bbu%time_c(i)
enddo

do i=1, c_num-1
   do j=1, m_size
      do k=1, m_size
         bbu%mat_c((i-1)*m_size*m_size+(j-1)*m_size+k) = erlmat(i+1,j,k)
      enddo
      if (output) then
         write(20,*) i
         write(20,'(4F24.12)')bbu%mat_c((i-1)*m_size*m_size+(j-1)*m_size+1), bbu%mat_c((i-1)*m_size*m_size+(j-1)*m_size+2), &
                           bbu%mat_c((i-1)*m_size*m_size+(j-1)*m_size+3),  bbu%mat_c((i-1)*m_size*m_size+(j-1)*m_size+4)
      endif
   enddo
enddo

if (output) then
  do i=1, homtotal
     write (20,*) i, Hpower(i)%h_ind, bbu%Q_c(i), bbu%RoQ_c(i), bbu%freq_c(i), bbu%angle_c(i)
  enddo
endif

if (output) close(20)

bbu%current = cur
bbu%n_amp = noise_amp
bbu%p_amp = pos_amp 
bbu%n_on = noise_on

deallocate (erlmat, erltime, cavityind)

end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine Tracking(bbu,power_t,coor_t,Hpower,output)
!
! Particle tracking at a given current
! 
! Modules Needed
!   use bbu_mod
!
! Input
!   bbu     -- bbu_info: containing all the information needed for bbu simulation
!   power_t -- an array containing the HOM power at different time 
!   coor_t  -- an array containing the particle orbit at different time
!   Hpower  -- HOM_power: to be used for storing the power of each HOM at the end of the tracking
!   output  -- logical: if true the HOM power at different time will be written in a file.
!
! Output to the screen
!-

subroutine Tracking(bbu,power_t,coor_t,Hpower,output)

implicit none

type (bbu_info), INTENT(INOUT) :: bbu
real(rp), dimension(:), allocatable,INTENT(INOUT) :: power_t
type (coord_struct), dimension(:), allocatable,INTENT(INOUT) :: coor_t
type (HOM_power), dimension(:), INTENT(INOUT) :: Hpower
logical, INTENT(IN) :: output

real(rp), dimension(:), allocatable :: coor_t_c
integer i,j,ind
integer :: matrixsize = 4

if (allocated(power_t)) deallocate(power_t)
if (allocated(coor_t))  deallocate(coor_t)
allocate(power_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t_c(size(coor_t)*matrixsize))

! Call the c++ bbu tracking subroutine with option 0

call bbu_track(bbu%current, bbu%cavity, bbu%hom_num, bbu%mat_c, bbu%time_c, bbu%hom, bbu%Q_c, bbu%RoQ_c, bbu%freq_c, bbu%angle_c, bbu%power_c, power_t,coor_t_c, bbu%b_freq, bbu%b_time, bbu%p_time, bbu%n_amp, bbu%p_amp, bbu%n_on, 0)

! Store the HOM power at the end of the tracking

do i=1, bbu%hom_num
   Hpower(i)%power = bbu%power_c(i)
enddo

! Store the particle coordinates at the end of the tracking

do i=1, size(coor_t)
  do j=1, matrixsize
     coor_t(i)%vec(j)=coor_t_c((i-1)*4+j)
  enddo
enddo

if (output) then
   
   OPEN(10, FILE='power_c', STATUS='UNKNOWN', ACTION='WRITE')
   OPEN(20, FILE='coor_t', STATUS='UNKNOWN', ACTION='WRITE')
   ind = 1
   do i=1, bbu%cavity
      do j=1, bbu%hom(i)
         if (j==1) then
             write(10,'(I6,3X,I6,3X,E20.14)',ADVANCE = "NO")i, Hpower(ind)%c_ind, Hpower(ind)%power
         else
             write(10,'(3X,E20.14)',ADVANCE = "NO")Hpower(ind)%power
         endif
         ind = ind + 1
      enddo
      write(10,*)''
   enddo

   do i=1,size(power_t)
       write(20,'(I6, 3X, E20.14,3X,E20.14,3X,E20.14,3X,E20.14,3X,E20.14)') &
        i, power_t(i), coor_t(i)%vec(1), coor_t(i)%vec(2), coor_t(i)%vec(3), coor_t(i)%vec(4)
   enddo

   close(10)
   close(20)

endif

end subroutine



!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine Threshold(bbu,power_t,coor_t,Hpower,output)
!
! Finding the bbu threshold current for a given lattice
! 
! Modules Needed
!   use bbu_mod
!
! Input
!   bbu     -- bbu_info: containing all the information needed for bbu simulation
!   cur     -- real: initial trial threshold current
!   power_t -- an array containing the HOM power at different time 
!   coor_t  -- an array containing the particle orbit at different time
!   Hpower  -- HOM_power: to be used for storing the power of each HOM at the end of the tracking
!   output  -- logical: if true the HOM power at different time will be written in a file.
!
! Output to the screen
!-   

subroutine Threshold(bbu,cur,power_t, coor_t, Hpower, output)

implicit none

type (bbu_info), INTENT(INOUT) :: bbu
real(rp),        INTENT(OUT)  :: cur
real(rp), dimension(:), allocatable,INTENT(INOUT) :: power_t
type (coord_struct), dimension(:), allocatable,INTENT(INOUT) :: coor_t
type (HOM_power), dimension(:), INTENT(INOUT) :: Hpower
logical, INTENT(IN) :: output
real(rp), dimension(:), allocatable :: coor_t_c
integer i,j,ind
integer :: matrixsize = 4

if (allocated(power_t)) deallocate(power_t)
if (allocated(coor_t))  deallocate(coor_t)
allocate(power_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t_c(size(coor_t)*matrixsize))

! Call the c++ bbu tracking subroutine with option 1

call bbu_track(cur, bbu%cavity, bbu%hom_num,bbu%mat_c, bbu%time_c, bbu%hom, bbu%Q_c, bbu%RoQ_c, bbu%freq_c, bbu%angle_c, bbu%power_c, power_t,coor_t_c, bbu%b_freq, bbu%b_time, bbu%p_time, bbu%n_amp, bbu%p_amp, bbu%n_on, 1)

! Store the HOM power at the end of the tracking

do i=1, bbu%hom_num
   Hpower(i)%power = bbu%power_c(i)
enddo

! Store the particle coordinates at the end of the tracking

do i=1, size(coor_t)
  do j=1, matrixsize
     coor_t(i)%vec(j)=coor_t_c((i-1)*4+j)
  enddo
enddo

if (output) then

   OPEN(10, FILE='power_c', STATUS='UNKNOWN', ACTION='WRITE')
   OPEN(20, FILE='coor_t', STATUS='UNKNOWN', ACTION='WRITE')
   ind = 1
   do i=1, bbu%cavity
      do j=1, bbu%hom(i)
         if (j==1) then
             write(10,'(I6,3X,I6,3X,E20.14)',ADVANCE = "NO")i, Hpower(ind)%c_ind, Hpower(ind)%power
         else
             write(10,'(3X, E20.14)',ADVANCE = "NO")Hpower(ind)%power   
         endif
         ind = ind + 1
      enddo
      write(10,*)''
   enddo

   do i=1,size(power_t)
       write(20,'(I6, 3X, E20.14,3X,E20.14,3X,E20.14,3X,E20.14,3X,E20.14)') &
          i, power_t(i), coor_t(i)%vec(1), coor_t(i)%vec(2), coor_t(i)%vec(3), coor_t(i)%vec(4)
   enddo

   close(10)
   close(20)

endif

end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine Fdep(bbu,power_t, coor_t, f_l, f_u, N, thresh_dep, output)
!
! Threshold currents for different frequency spreads
! 
! Modules Needed
!   use bbu_mod
!
! Input
!   bbu     -- bbu_info: containing all the information needed for bbu simulation
!   power_t -- an array containing the HOM power at different time 
!   coor_t  -- an array containing the particle orbit at different time
!   f_l     -- real: lower limit of the frequency spread
!   f_u     -- real: upper limit of the frequency spread
!   N       -- integer: the number of divisions between f_l and f_u
!   thresh_dep -- pair_num: a list of threshold currents and their corresponding polarization angles
!   output  -- logical: if true the HOM power at different time will be written in a file.
!
! Output to the screen
!-

subroutine Fdep(bbu,power_t, coor_t, f_l, f_u, N, thresh_dep, output)

use random_mod

implicit none

type (bbu_info), INTENT(INOUT) :: bbu
real(rp), dimension(:), allocatable,INTENT(INOUT) :: power_t
type (coord_struct), dimension(:), allocatable,INTENT(INOUT) :: coor_t
real(rp), INTENT(IN) :: f_l
real(rp), INTENT(IN) :: f_u
integer,  INTENT(IN) :: N
type (pair_num), INTENT(INOUT) :: thresh_dep
logical, INTENT(IN) :: output
real(rp), dimension(:), allocatable :: coor_t_c
integer :: matrixsize = 4
real(rp), DIMENSION(:), ALLOCATABLE :: Newf_c
real(rp), DIMENSION(:), ALLOCATABLE :: hompower
real(rp)  r, rr, spread
integer   i,j,k

if (allocated(power_t)) deallocate(power_t)
if (allocated(coor_t))  deallocate(coor_t)
allocate(power_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t_c(size(coor_t)*matrixsize))
allocate(Newf_c(bbu%hom_num))
allocate(hompower(size(bbu%power_c)))

if (output) then
   OPEN(10, FILE='FDEP', STATUS='UNKNOWN', ACTION='WRITE')
endif

! Assign values to the HOM frequencies

r = 1.0/(N-1)
do i=0, N-1

   do k=1, size(bbu%power_c)
      hompower(k) = bbu%power_c(k)
   enddo

   spread = f_l + i*(f_u-f_l)*r
   do k=1, bbu%hom_num
      call ran_gauss(rr)
      Newf_c(k) = bbu%freq_c(k) * (1 + spread * rr)
   enddo

! Call the c++ bbu tracking subroutine with option 1

   call bbu_track(bbu%current,bbu%cavity,bbu%hom_num,bbu%mat_c,bbu%time_c,bbu%hom,bbu%Q_c,bbu%RoQ_c,Newf_c,bbu%angle_c,hompower,power_t,coor_t_c,bbu%b_freq,bbu%b_time,bbu%p_time,bbu%n_amp,bbu%p_amp,bbu%n_on,1)

! Record the threshold current and its corresponding frequency spread
 
   if (output) then
       write(10,'(E20.14, 3X, E20.14)')f_l+i*(f_u-f_l)*r, bbu%current
   endif
   thresh_dep%var(i+1) = f_l+i*(f_u-f_l)*r
   thresh_dep%threshold(i+1) = bbu%current
enddo

if (output) then
   close(10)
endif

do i=1, size(coor_t)
  do j=1, matrixsize
     coor_t(i)%vec(j)=coor_t_c((i-1)*4+j)
  enddo
enddo

end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!-
! Subroutine Pdep(bbu,power_t, coor_t, a_l, a_u, check, N, thresh_dep, output)
!
! Threshold currents for different polarization angles
! 
! Modules Needed
!   use bbu_mod
!
! Input
!   bbu     -- bbu_info: containing all the information needed for bbu simulation
!   power_t -- an array containing the HOM power at different time 
!   coor_t  -- an array containing the particle orbit at different time
!   a_l     -- real: lower limit of the polarization angle
!   a_u     -- real: upper limit of the polarization angle
!   check   -- integer: if 0, then change the polarization angle of every HOM
!   N       -- integer: the number of divisions between a_l and a_u
!   thresh_dep -- pair_num: a list of threshold currents and their corresponding polarization angles
!   output  -- logical: if true the HOM power at different time will be written in a file.
!
! Output to the screen
!+

subroutine Pdep(bbu,power_t, coor_t, a_l, a_u, check, N, thresh_dep, output)

implicit none


type (bbu_info), INTENT(INOUT) :: bbu
real(rp), dimension(:), allocatable, INTENT(INOUT) :: power_t
type (coord_struct), dimension(:), allocatable, INTENT(INOUT) :: coor_t
real(rp), INTENT(IN) :: a_l
real(rp), INTENT(IN) :: a_u
integer, INTENT(IN) :: check
integer, INTENT(IN) :: N
type (pair_num), INTENT(INOUT) :: thresh_dep
logical, INTENT(IN) :: output
real(rp), dimension(:), allocatable :: coor_t_c
integer :: matrixsize = 4
real(rp), DIMENSION(:), ALLOCATABLE :: Newa_c
real(rp), DIMENSION(:), ALLOCATABLE :: hompower
integer, DIMENSION(:), ALLOCATABLE :: start
real(rp)  r
integer   i,j,k

if (allocated(power_t)) deallocate(power_t)
if (allocated(coor_t))  deallocate(coor_t)
allocate(power_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t_c(size(coor_t)*matrixsize))
allocate(Newa_c(bbu%hom_num))
allocate(start(bbu%cavity))
allocate(hompower(size(bbu%power_c)))

if (output) then
    OPEN(10, FILE='PDEP', STATUS='UNKNOWN', ACTION='WRITE')
endif

! Find the index of the first HOM in each cavity.

start(1) = 1
do i=2, bbu%cavity
   start(i) = start(i-1) + bbu%hom(i-1)
enddo

! Assign values to the HOM polarization angle

r = 1.0/(N-1)
do i=0, N-1
  do k=1, size(bbu%power_c)
      hompower(k) = bbu%power_c(k)
  enddo
  if (check==0) then                ! check = 0: change the polariation angles of all HOMs
   do k=1, bbu%hom_num
      Newa_c(k) = (a_l + i*(a_u-a_l)*r)/180*pi
   enddo
  else
   do k=1, bbu%hom_num
      Newa_c(k) = bbu%angle_c(k)    
   enddo
   do k=1, bbu%cavity
      Newa_c(start(k)+check-1) = (a_l + i*(a_u-a_l)*r)/180*pi     ! check = n: change the polarization angle of the n-th HOM
   enddo
  endif

! Call the c++ bbu tracking subroutine with option 1
   
  call bbu_track(bbu%current,bbu%cavity,bbu%hom_num,bbu%mat_c,bbu%time_c,bbu%hom,bbu%Q_c,bbu%RoQ_c,bbu%freq_c,Newa_c,hompower,power_t,coor_t_c,bbu%b_freq,bbu%b_time,bbu%p_time,bbu%n_amp,bbu%p_amp,bbu%n_on,1)

! Record the threshold current and its corresponding polarization angle
  
  if (output) then
     write(10,'(E20.14, 3X, E20.14)')(a_l+i*(a_u-a_l)*r)/180*pi, bbu%current
  endif
  thresh_dep%var(i+1) = (a_l+i*(a_u-a_l)*r)/180*pi
  thresh_dep%threshold(i+1) = bbu%current  
enddo

if (output) then
    close(10)
endif

do i=1, size(coor_t)
  do j=1, matrixsize
     coor_t(i)%vec(j)=coor_t_c((i-1)*4+j)
  enddo
enddo

end subroutine

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine Qdep(bbu,power_t, coor_t, Q_l, Q_u, ck, N, thresh_dep, output)
!
! Threshold currents for different Qs
! 
! Modules Needed
!   use bbu_mod
!
! Input
!   bbu     -- bbu_info: containing all the information needed for bbu simulation
!   power_t -- an array containing the HOM power at different time 
!   coor_t  -- an array containing the particle orbit at different time
!   Q_l     -- real: lower limit of Q
!   Q_u     -- real: upper limit of Q
!   ck      -- integer: if 0, then change the polarization angle of every HOM
!   N       -- integer: the number of divisions between a_l and a_u
!   thresh_dep -- pair_num: a list of threshold currents and their corresponding polarization angles
!   output  -- logical: if true the HOM power at different time will be written in a file.
!
! Output to the screen
!-

subroutine Qdep(bbu,power_t,coor_t,Q_l,Q_u,ck,N,thresh_dep, output)

implicit none

type (bbu_info), INTENT(INOUT) :: bbu
real(rp), dimension(:), allocatable, INTENT(INOUT) :: power_t
type (coord_struct), dimension(:), allocatable, INTENT(INOUT) :: coor_t
real(rp), INTENT(IN) :: Q_l
real(rp), INTENT(IN) :: Q_u
integer, INTENT(IN) :: ck
integer, INTENT(IN) :: N
type (pair_num), INTENT(INOUT) :: thresh_dep
logical, INTENT(IN) :: output
real(rp), dimension(:), allocatable :: coor_t_c
integer :: matrixsize = 4
real(rp), DIMENSION(:), ALLOCATABLE :: hompower
integer, DIMENSION(:), ALLOCATABLE :: start
real(rp), DIMENSION(:), ALLOCATABLE :: NewQ_c
real(rp)  r
integer   i,j,k

if (allocated(power_t)) deallocate(power_t)
if (allocated(coor_t))  deallocate(coor_t)
allocate(power_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t(bbu%b_time/bbu%p_time+1))
allocate(coor_t_c(size(coor_t)*matrixsize))
allocate(NewQ_c(bbu%hom_num))
allocate(start(bbu%cavity))
allocate(hompower(size(bbu%power_c)))

if (output) then
    OPEN(10, FILE='QDEP', STATUS='UNKNOWN', ACTION='WRITE')
endif

! Find the index of the first HOM in each cavity.

start(1) = 1
do i=2, bbu%cavity
   start(i) = start(i-1) + bbu%hom(i-1)
enddo

! Assign values to the HOM Q

r = 1.0/(N-1)
do i=0, N-1
   do k=1, size(bbu%power_c)
      hompower(k) = bbu%power_c(k)
   enddo
  if (ck==0) then      ! check = 0: change the Qs of all HOMs
   do k=1, bbu%hom_num
      NewQ_c(k) = Q_l + i*(Q_u-Q_l)*r
   enddo
  else                 ! check = n: change the Q of the n-th HOM
   do k=1, bbu%hom_num
      NewQ_c(k) = bbu%Q_c(k)
   enddo
   do k=1, bbu%cavity
      NewQ_c(start(k)+ck-1) = Q_l + i*(Q_u-Q_l)*r
   enddo
  endif

! Call the c++ bbu tracking subroutine with option 1

  call bbu_track(bbu%current,bbu%cavity,bbu%hom_num,bbu%mat_c,bbu%time_c,bbu%hom,NewQ_c,bbu%RoQ_c,bbu%freq_c,bbu%angle_c,hompower,power_t,coor_t_c,bbu%b_freq,bbu%b_time,bbu%p_time,bbu%n_amp,bbu%p_amp,bbu%n_on,1)

! Record the threshold current and its corresponding Q value.
  
  if (output) then
      write(10,'(E20.14, 3X, E20.14)')Q_l+i*(Q_u-Q_l)*r, bbu%current
  endif
  thresh_dep%var(i+1) = Q_l+i*(Q_u-Q_l)*r
  thresh_dep%threshold(i+1) = bbu%current     
enddo

if (output) then
   close(10)
endif

do i=1, size(coor_t)
  do j=1, matrixsize
     coor_t(i)%vec(j)=coor_t_c((i-1)*4+j)
  enddo
enddo

end subroutine

end module
