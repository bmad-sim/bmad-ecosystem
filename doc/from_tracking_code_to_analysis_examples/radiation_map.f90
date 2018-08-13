program radiation_map
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
real(dp) prec,closed_orbit(6),mat(6,6),a(6,6),ai(6,6),del,error(6),emi(3),sij(3)
real(dp) Ia(6,6,3),Sa(6,6,3),S(6,6),Ka(6,6,3),Ba(6,6,3),Ha(6,6,3),Ea(6,6,3),xij 
complex(dp) mc(6,6)
type(internal_state),target :: state,sta(2)
logical(lp) :: mis=.false. 
type(c_damap)  one_turn_map, Id,a0,a_cs
type(real_8) y(6)
type(c_normal_form) normal_form
type(c_taylor) e(3),z_ij,ave_FLOQUET
integer expo(6),pos
character*48 command_gino,fmd,fmd1
integer i,map_order,j,mf,k,t,cas
type(probe) ray_closed
type(probe_8) ray
!!!!!!!!!!!!!!!!!!!!!
type(fibre), pointer :: p1
type(integration_node), pointer :: ti
!!!!!!!!!!!!!!!!!!!!!
c_verbose=.false.
prec=1.d-10 ! for printing
longprint=.false. 
 del=0.d0
use_info=.true.
fmd= '(a12,1X,a3,I1,a3,i1,a4,D18.11,1x,D18.11)'
fmd1='(1X,a3,I1,a3,i1,a4,2(D18.11,1x),(f10.3,1x),a2)'
sta(1)=default0
sta(2)=default0 +radiation0+envelope0

Ia=0.d0; Sa=0.d0; S=0.d0;Ba=0.d0;Ha=0.d0;Ea=0.d0

do i=1,3
 S(2*i-1,2*i) = 1.d0; S(2*i,2*i-1) = -1.d0;
 Ia(2*i-1,2*i-1,i) =1.d0;  Ia(2*i,2*i,i)   = 1.d0;
 Sa(2*i-1,2*i,i)   =1.d0;  Sa(2*i,2*i-1,i) =-1.d0;
enddo

call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start

write(6,*) " misalignments type t, otherwise f"
read(5,*) mis

If(mis) then
 write(6,'(a56,/)') " The lattice produced will have errors and thus coupling"
else
 write(6,'(a56,/)') " The lattice produced is ideal: mid-plane symmetric     "
endif

pos=15
error=0.d0
error(6)=1.d-5 ! tilt along z-axis

call build_lattice_als(ALS,mis,error,exact=.false.) 

call kanalnummer(mf,"result_with_stochastic_radiation.txt")

do cas=1,2

 state=sta(cas)

map_order=1
call init_all(state,map_order,0)

call alloc(one_turn_map,id,a0,a_cs)
call alloc(y) 
call alloc(normal_form)
call alloc(e)
call alloc(z_ij,ave_FLOQUET)
call alloc(ray)

closed_orbit=0.d0;                                                   ! (1)

call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=pos)   ! (2)
if(cas==2) then
 write(mf,*) " ";write(mf,*) "  Constant part of the map";
 write(mf,'(a16,6(1x,g12.5))') " closed orbit = ",closed_orbit(1:6)    
endif 

ray_closed=closed_orbit     ! (3)
id=1;   
! ray= closed orbit + identity map  
ray=ray_closed+id;          ! (4)
                                  

call propagate(als,RAY,state,fibre1=pos)  ! (5)

! Six polymorphs and the fluctuationsare E_ij 
! are promoted to Taylor maps  

one_turn_map=ray                         ! (6a)
one_turn_map=one_turn_map.sub.1          ! (6b)    

if(cas==2) then
 write(mf,*) " "; write(mf,*) "  The linear map";
 call print(one_turn_map,mf,prec)                    ! (7)                                  
endif

call  c_normal(one_turn_map,normal_form)             ! (8)

if(cas==2) write(mf,'(a8,3(1x,g20.13))') " tune = ",normal_form%tune(1:3)     ! (9)


a=normal_form%a_t                                                    ! (10a)
ai=normal_form%a_t**(-1)                                             ! (10b)

mc=from_phasor(-1)*normal_form%a_t**(-1)*one_turn_map*normal_form%a_t*from_phasor()  ! (11)

if(cas==2) then
write(mf,'(/,a37,15x,a41,/)')" i  j    Real(mc(i,j))    Im(mc(i,j))", &
                             "damping  + i *   phase  ( =log(mc(i,j)) )"
do i=1,c_%nd2
do j=1,c_%nd2
if(abs(mc(i,j))>1.e-10_dp) write(mf,'(i2,1x,i2,2(5x,g12.5),5x,2(5x,g12.5))')&
                                      i,j,mc(i,j),log(mc(i,j))                       ! (12)
enddo 
enddo
endif

if(cas==1) then
!!!!!!! Lattice functions !!!!!!!
!! coefficient of invariant  !!

do i=1,c_%nd
 Ba(1:6,1:6,i)=matmul(matmul(a,Sa(1:6,1:6,i)),ai)  ! (13a)
 Ha(1:6,1:6,i)=matmul(matmul(a,Ia(1:6,1:6,i)),ai)  ! (13b)
 Ka(1:6,1:6,i)=-matmul(S,Ba(1:6,1:6,i))            ! (13c)
 Ea(1:6,1:6,i)=-matmul(Ba(1:6,1:6,i),S)            ! (13d)
enddo

endif
 

if(cas==2) then 
write(mf,*)
write(mf,'(16X,a50)') "   Equilibrium moments in Phasors Basis           "
 do i=1,6
 do j=i,6 
  if(abs(normal_form%s_ijr(i,j))>1.d-20) then
   write(mf,fmd) " Phasors -> ","<x_",i," x_",j,"> = ",  &   ! (14)
                    c_clean(normal_form%s_ijr(i,j),1.d-20)  
  endif
 enddo
 enddo 

emi(1)=real(normal_form%s_ijr(1,2))/2        ! (15a)
emi(2)=real(normal_form%s_ijr(3,4))/2        ! (15b)
emi(3)=real(normal_form%s_ijr(5,6))/2        ! (15c)

write(mf,*) 

write(mf,'(16X,a54)') "    Exact              Chao         (Exact-Chao)/Exact"
 do i=1,6
 do j=1,6
   xij=0.d0
   do k=1,3
     xij= Ea(i,j,k)*emi(k) + xij                              ! (16)
   enddo
 if(abs(normal_form%s_ij0(i,j))>1.e-15_dp) then
   write(mf,fmd1) "<x_",i," x_",j,"> = ", real(normal_form%s_ij0(i,j)), &
                 xij ,abs(100*(real(normal_form%s_ij0(i,j))-xij)/       &
   real(normal_form%s_ij0(i,j)))," %"
 endif
 enddo
 enddo 

 sij(1)=emi(1)*Ea(1,6,1)+emi(2)*Ea(1,6,2)+emi(3)*Ea(1,6,3)    ! (17a)
 sij(2)=emi(1)*Ea(6,6,1)+emi(2)*Ea(6,6,2)+emi(3)*Ea(6,6,3)    ! (17b)
 sij(3)=emi(1)*Ea(1,1,1)+emi(2)*Ea(1,1,2)+emi(3)*Ea(1,1,3)    ! (17c)

 write(mf,*) 
 write(mf,'(a56)') " Exact Crab Angle with beam envelope                    "
 write(mf,'((5x,D18.11))') real(normal_form%s_ij0(1,6)) &
/(real(normal_form%s_ij0(6,6))-real(normal_form%s_ij0(1,1)))                  ! (18)
 write(mf,'(a25)') " Chao Exact Crab Angle " 
 write(mf,'((5x,D18.11))') sij(1)/(sij(2)-sij(3))                             ! (19)
 write(mf,'(a56)') " Almost Exact Crab Angle ignoring transverse emittances "    
 write(mf,'((5x,D18.11))') -Ka(5,2,3)/(Ka(5,5,3)-Ka(2,2,3))                   ! (20)
 write(mf,'(a56)') " Approximate Crab Angle ignoring transverse emittances  "
 write(mf,'((5x,D18.11))')  Ha(1,6,3)-Ka(5,6,3)*Ha(1,5,3)/Ka(5,5,3)           ! (21)






write(mf,*); 
write(mf,*) "!!!  Raising the maps with moments to the power 2**100  !!!"

do i=1,100
  one_turn_map=one_turn_map*one_turn_map  ! (22a)
enddo

call print(one_turn_map,mf)               ! (22b)

endif


call kill(one_turn_map,id,a0,a_cs)
call kill(y) 
call kill(normal_form)
call kill(e)
call kill(z_ij,ave_FLOQUET)
call kill(ray)
enddo



close(mf)

call ptc_end(graphics_maybe=1,flat_file=.false.)

end program radiation_map

