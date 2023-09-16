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
type(c_damap)  one_turn_map, quasi_diagonal, diagonal,a_ac, id
type(c_normal_form) normal_form
integer :: pos =1 
integer i,map_order,mf1,mfmap
type(probe) ray_closed
type(probe_8) ray
type(real_8) y(6)
complex(dp) g1,g2,al1,al2
!!!!!!!!!!!!!!!!!!!!!


c_verbose=.false.
prec=1.d-10 ! for printing
longprint=.false. 

call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start


 
call kanalnummer(mfmap,"maps.txt") 

state= nocavity0 + modulation0 !

map_order=1
call init_all(state,map_order,0)


call alloc(y)
call alloc(one_turn_map,quasi_diagonal, diagonal,id,a_ac)
call alloc(normal_form)
call alloc(ray)


call build_lattice_als(ALS,mis,exact=.false.) 

!!!! circ is the circumference of the ring !!!! 
call get_length(als,circ)
!!!! AC_modulate.txt sets the magnet QF1 as a modulated magnet !!!! 
 call kanalnummer(mf1,file="AC_modulation.txt")
  write(mf1,*) "select layout"                  
  write(mf1,*) 1
  write(mf1,*) " MODULATE"               
  write(mf1,*) " BEND1 1" ! name and number of frequencies              
  write(mf1,*) "1.d0 0 0       !DC_ac,A_ac,theta_ac"
  write(mf1,*) "1.d0   1       ! D_ac,n_ac  "
  write(mf1,*) "1 0.001d0 0      ! n d_bn(n) d_an(n)  "  ! (A)
  write(mf1,*) "0  0 0 "
  write(mf1,*) " return "
 close(mf1)
 call read_ptc_command77("AC_modulation.txt")

!!!! set a modulation clock !!!!!!
mu_mod=twopi*0.12345d0; 
ray_closed%ac(1)%om=mu_mod/circ ! (B1) differs from the first edition 
ray_closed%ac(1)%x=0.d0 ;       ! (B2) differs from the first edition 
write(6,*) " Modulation tune in radians =",circ*ray_closed%ac(1)%om

closed_orbit=0.d0;                                                   ! (C)

call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=pos)   ! (D)
 
ray_closed=closed_orbit     ! (E)

id=1;    
! ray= closed orbit + identity map  

ray=id+ray_closed;          ! (F)
                  
call propagate(als,RAY,state,fibre1=pos)  ! (G)
 
 one_turn_map=ray                         ! (H)
 write(mfmap,*); write(mfmap,*) " Map produced by code " ; write(mfmap,*);  
 call print(one_turn_map,mfmap)

do_linear_ac_longitudinal=.false.           ! (I)
call  c_normal(one_turn_map,normal_form)    ! (J1)
id=normal_form%a_t*from_phasor()
do_linear_ac_longitudinal=.true. 
call  c_normal(one_turn_map,normal_form)    ! (J2)
normal_form%a_t=normal_form%a_t*from_phasor()

write(mfmap,*);write(mfmap,*);
write(mfmap,*) " Correct A (M=ARA^-1) from the algorithm for the code ";
write(mfmap,*);write(mfmap,*);

call  print(normal_form%a_t,mfmap,prec)  
 
diagonal=normal_form%a_t**(-1)*one_turn_map*normal_form%a_t   ! (K1)

quasi_diagonal=id**(-1)*one_turn_map*id                       ! (K2)

g1=quasi_diagonal%v(6).sub.'00000010'                         ! (L1)
g2=quasi_diagonal%v(6).sub.'00000001'                         ! (L2)

al1=-g1/(1.d0-(quasi_diagonal%v(7).sub.'00000010'))           ! (L3)
al2=-g2/(1.d0-(quasi_diagonal%v(8).sub.'00000001'))           ! (L4)

a_ac=1

a_ac%v(6)=a_ac%v(6)+(al1.cmono.'00000010')+(al2.cmono.'00000001') ! (L5)

!!!!!!         Print all the resulting maps         !!!!!!

write(mfmap,*);write(mfmap,*);
write(mfmap,*) " Correct R (M=ARA^-1) from the algorithm for the code ";
write(mfmap,*);write(mfmap,*);
 

call  print(diagonal,mfmap,prec)  

write(mfmap,*);write(mfmap,*);
write(mfmap,*) " Quasi-diagonal R (M=ARA^-1) from the algorithm for the code ";
write(mfmap,*);write(mfmap,*);
 
call  print(quasi_diagonal,mfmap,prec)  

quasi_diagonal=a_ac**(-1)*quasi_diagonal*a_ac

write(mfmap,*);write(mfmap,*);
write(mfmap,*) " Quasi-diagonal should now be diagonal ";
write(mfmap,*);write(mfmap,*);
 

call  print(quasi_diagonal,mfmap,prec)  

call kill(y)
call kill(one_turn_map,quasi_diagonal, diagonal,id,a_ac)
call kill(normal_form)
call kill(ray)
 
 close(mfmap)
call ptc_end(graphics_maybe=1,flat_file=.false.)

end program modulated_map

