program modulated_map
use madx_ptc_module
use pointer_lattice
use c_TPSA
implicit none

    interface
       subroutine build_lattice_als(ALS,MIS,error,exact,sl,thin,onecell)
         use madx_ptc_module
         use pointer_lattice
         implicit none
         type(layout), target :: ALS
         logical(lp) mis
         real(dp),optional :: error(6)
         logical, optional :: exact,sl,thin,onecell
       end subroutine build_lattice_als
    end interface
type(layout), pointer:: ALS
real(dp) prec,closed_orbit(6) 
real(dp) L,Kq,k0,mu_mod,beta,dmu,mu_x,circ
 
type(internal_state),target :: state 
logical(lp) :: mis=.false. 
type(c_damap)  one_turn_map, drift_map, quad_map,id
type(c_normal_form) normal_form
type(c_taylor) q(4)
integer :: pos =1, nind(11)
integer i,map_order,mf,mf1,mfmap
type(probe) ray_closed
type(probe_8) ray
!!!!!!!!!!!!!!!!!!!!!


c_verbose=.false.
prec=1.d-10 ! for printing
longprint=.false. 
 
state=only_2d0+modulation0

call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start

! sl= true produces the small 1-d lattice
!  Qf1 = QUADRUPOLE(" QF1 ",0.d0, K1= 0.01d0 ); L1  = drift("L1 ",0.1d0);
! ALS=L1+QF1;

!call build_lattice_als(ALS,mis,exact=.false.,sl=.true.) 

call kanalnummer(mf,"result_with_modulation.txt")
call kanalnummer(mfmap,"maps.txt") 

!!!! Let us do the baby example with maps
map_order=3
call init_all(state,map_order,0)

call alloc(one_turn_map, drift_map, quad_map,id)
call alloc(q,4)
call alloc(normal_form)
call alloc(ray)

!!!!!!! Simple example !!!!!!! 
!  without using PTC
do i=1,4
 q(i)=1.d0.cmono.i   ! q_1, q_2, q_3 and q_4 are created as TPSA variables
enddo
l=0.1d0; Kq=.1d0; k0=.2d0; mu_mod=twopi*0.12345d0;  ! (1a)
! Analytic calculation 
beta=1.d0/sqrt(Kq)/sqrt(1.d0-Kq*L**2/4)             ! (1b)

mu_x=acos(1.d0-Kq*L**2/2)                           !  (2a)
dmu=-2*(sin(2*mu_x+mu_mod)/(1-cos(2*mu_x+mu_mod)) & !  (2b)
+sin(2*mu_x-mu_mod)/(1-cos(2*mu_x-mu_mod)))*4*(-beta*k0*L/16)**2

write(mf,*);
write(mf,*) "  Analytical tune in radians =     ",mu_x
write(mf,*) "  Analytical tune shift in radians = ",dmu
write(mf,*);

drift_map=1
quad_map=1

drift_map%v(1)=q(1)+L*q(2)    ! (3) drift
drift_map%v(2)=q(2)  
drift_map%v(3)=cos(mu_mod)*q(3) + sin(mu_mod)*q(4)
drift_map%v(4)=cos(mu_mod)*q(4) - sin(mu_mod)*q(3) 

quad_map%v(1)=q(1) 
quad_map%v(2)=q(2)-L*(Kq+k0*q(3))*q(1) ! (4) quadrupole

! Map of system is made
one_turn_map=quad_map*drift_map   ! total map   (5)

 write(mfmap,*); write(mfmap,*) " Map hardwired in main program " ; write(mfmap,*);
 call print(one_turn_map,mfmap)

! Map is normalised
call  c_normal(one_turn_map,normal_form)      !  (6)
write(mf,*);
write(mf,*) " Result from the normal form algorithm for hardwired map ";
write(mf,*);

write(mf,*) " Normal form result for tune in radians =       ", &  !  (7)
 -aimag(normal_form%ker%f(1)%v(1).sub.'1000')   
write(mf,*) " Normal form result for tune shift in radians = ", &  !  (8)
 -aimag(normal_form%ker%f(3)%v(1).sub.'1011')
 
! Same calculation from with PTC
! using the following little lattice
!  Qf1 = QUADRUPOLE(" QF1 ",L=0.d0, K1= 0.01d0 ); L1  = drift("L1 ",L0=.1d0);
!  ALS=L1+QF1;
call build_lattice_als(ALS,mis,exact=.false.,sl=.true.) 

!!!! circ is the circumference of the ring !!!! 
call get_length(als,circ)
!!!! AC_modulate.txt sets the magnet QF1 as a modulated magnet !!!! 
 call kanalnummer(mf1,file="AC_modulation.txt")
  write(mf1,*) "select layout"                  
  write(mf1,*) 1
  write(mf1,*) " MODULATE"               
  write(mf1,*) " QF1  1"       ! name and number of clocks           
  write(mf1,*) "1.d0 0 0       !DC_ac,A_ac,theta_ac"
  write(mf1,*) "1.d0   2       ! D_ac,n_ac  "
  write(mf1,*) "2 0.02d0 0      ! n d_bn(n) d_an(n)  "  ! (A)
  write(mf1,*) "0  0 0 "
  write(mf1,*) " return "
 close(mf1)
 call read_ptc_command77("AC_modulation.txt")

!!!! set a modulation clock !!!!!!
ray_closed%ac(1)%om=mu_mod/circ ! (B1) ! differs from the first edition   
ray_closed%ac(1)%x=0.d0 ;       ! (B2) ! differs from the first edition   
write(6,*) " Modulation tune in radians =",circ*ray_closed%ac%om

closed_orbit=0.d0;                                                   ! (C)

call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=pos)   ! (D)
 
ray_closed=closed_orbit     ! (E)
id=1;   
! ray= closed orbit + identity map  
ray=id+ray_closed;          ! (F1)
write(mf,*); write(mf,*) " Initial value of the clock ( in type probe_8) "
write(mf,*);call print(ray%ac(1),mf ) ! (F2) ! differs from the first edition   
                          
call propagate(als,RAY,state,fibre1=pos)  ! (G)

! Six polymorphs and the fluctuationsare E_ij 
! are promoted to Taylor maps  
 
one_turn_map=ray                         ! (H)
 write(mfmap,*); write(mfmap,*) " Map produced by code " ; write(mfmap,*);  ! (I)
 call print(one_turn_map,mfmap)
 
call  c_normal(one_turn_map,normal_form)    ! (J)

write(mf,*);write(mf,*);
write(mf,*) " Result from the normal form algorithm for the code ";
write(mf,*) 


write(mf,*) " Normal form result for tune in radians =       ", &    ! (K)
-aimag(normal_form%ker%f(1)%v(1).sub.'1000')
write(mf,*) " Normal form result for tune shift in radians = ", &    ! (L)
-aimag(normal_form%ker%f(3)%v(1).sub.'1011')
write(mf,*)





call kill(one_turn_map, drift_map, quad_map,id)
call kill(normal_form)
 
call kill(ray)
 


close(mfmap)
close(mf)

call ptc_end 

end program modulated_map

