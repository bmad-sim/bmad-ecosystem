! This is the main program for transverse bbu simulation
! 
! 1. Particle tracking at a given current.
! 2. Finding the threshold current for a given lattice.
! 3. Finding the dependence of the threshold current on the HOM Q.
! 4. Finding the dependence of the threshold current on the HOM frequency spread.
! 5. Finding the dependence of the threshold current on the HOM polarization angle.
!
! Module Needed:
! use bmad
! use bbu_mod
!
! Input from the screen
!   response1 -- Type of simulation
!   cur       -- Initial trial threshold current
!   response2 -- Output choice
!
! Output to the screen
! 

program csbbu


use bbu_mod

implicit none

type (bbu_info) bbu
type (lat_struct) lattice
type (pair_num) thresh_dep
type (HOM_power), dimension(:), allocatable :: Hpower
integer response1, response2, ck,Np,i
logical plot,output
real(rp) cur,Q_l,Q_u,f_l,f_u,a_l,a_u
real(rp), dimension(:), allocatable ::  power_t
type (coord_struct), dimension(:), allocatable ::  coor_t
integer Nt, n_prt
real, parameter :: bunch_freq = 1.3E9
real, parameter :: beam_time = 1.0E-4
real, parameter :: print_time = 1.0E-8

! Read in the lattices

call bmad_parser ("erl.lat", lattice)

! Calculate the total number of bunches

Nt = NINT(beam_time*bunch_freq)
n_prt = NINT(print_time*bunch_freq)

bbu%b_time = Nt
bbu%p_time = n_prt
bbu%b_freq = bunch_freq

! Read in simulation choices

print *, 'What would you like to do with this lattice'
print *, '1. Do the tracking at a given current'
print *, '2. Find the threshold current of a given lattice'
print *, '3. Find the dependence of the threshold current on the HOM Q.'
print *, '4. Find the dependence of the threshold current on the HOM frequency spread.'
print *, '5. Find the dependence of the threshold current on the HOM polarization angle.'
read *, response1

print *, 'Please input the trial current in mA'
read *, cur


print *, 'Do you want to check the plot of Power v.s. t at the end?'
print *, 'If you do, enter 1, otherwise enter 2'
read *, response2 



if (response2 == 1) then 
    plot = .true.
    output = .true.
else 
    plot = .false.
    output = .false.
endif

! Extract lattice information

call Get_info(lattice,cur,bbu,Hpower,output)

! Particle tracking

select case (response1)

case(1)  !  Do the tracking at a given current
 
   call Tracking(bbu,power_t,coor_t,Hpower,output)
   if (plot ==.true.) then
      call check_result1(power_t,size(power_t))
   endif

case(2)  !  Find the threshold current of a given lattice
  
   call Threshold(bbu,cur,power_t,coor_t,Hpower,output)
   if (plot ==.true.) call check_result1(power_t,size(power_t))

case(3)  !  Find the dependence of the threshold current on the HOM Q
 
   Q_l = 10000    ! lower limit of Q values 
   Q_u = 40000      ! upper limit of Q values
   ck  = 1          ! ck = 0 -- change the Qs of all HOMs, ck = n -- change the Q of the n-th HOM 
   Np  = 10         ! Np points are sampled within [Q_l, Q_u]
   allocate(thresh_dep%var(Np))
   allocate(thresh_dep%threshold(Np))
   call Qdep(bbu,power_t,coor_t,Q_l,Q_u,ck,Np,thresh_dep,output)
   do i=1,Np
      write(*,*) thresh_dep%var(i), thresh_dep%threshold(i)
   enddo  
   if (plot ==.true.) call check_result2(thresh_dep%var,thresh_dep%threshold,Np)

case(4)  !  Find the dependence of the threshold current on the HOM frequency spread
 
!   f_l = 5.33635015E-3         ! lower limit of the frequency spread
!   f_u = 5.33635015E-3         ! upper limit of the frequency spread
   f_l = 0.         ! lower limit of the frequency spread
   f_u = 0.         ! upper limit of the frequency spread
!   Np  = 500                   ! Np points are sampled within [f_l, f_u]
   Np  = 1                   ! Np points are sampled within [f_l, f_u]
   allocate(thresh_dep%var(Np))
   allocate(thresh_dep%threshold(Np))
   call Fdep(bbu,power_t,coor_t,f_l,f_u,Np,thresh_dep,output)
   do i=1,Np
      write(*,*) thresh_dep%var(i), thresh_dep%threshold(i)
   enddo  
   if (plot ==.true.) call check_result2(thresh_dep%var,thresh_dep%threshold,Np)

case(5)  !  Find the dependence of the threshold current on the HOM polarization angle
   
   a_l = 0                  ! lower limit of polarization angles 
   a_u = 90                 ! upper limit of polarization angles
   ck  = 1                  ! ck = 0 -- change the polarization angles of all HOMs, ck = n -- change the polarization angle of the n-th HOM 
   Np  = 10                 ! Np points are sampled within [a_l, a_u]
   allocate(thresh_dep%var(Np))
   allocate(thresh_dep%threshold(Np))
   call Pdep(bbu,power_t,coor_t,a_l,a_u,ck,Np,thresh_dep,output)
   do i=1,Np
      write(*,*) thresh_dep%var(i), thresh_dep%threshold(i)
   enddo  
   if (plot ==.true.) call check_result2(thresh_dep%var,thresh_dep%threshold,Np)
end select

end program


