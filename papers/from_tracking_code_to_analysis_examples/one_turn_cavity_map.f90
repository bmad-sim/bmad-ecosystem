program one_turn_cavity_map
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
real(dp) prec,closed_orbit(6),mat(6,6),a(6,6),ai(6,6),del,error(6),mom(6,6),ray(5,6),cav_tune
real(dp) Ia(6,6,3),Sa(6,6,3),S(6,6),Ka(6,6,3),Ba(6,6,3),Ha(6,6,3),Ea(6,6,3),inv,inv1,tune_s
REAL(DP) F(6,6),alpha_s,m(6,6), alphas(2),c
complex(dp) mc(6,6)
type(internal_state),target :: state
logical(lp) :: mis=.false. 
type(c_damap)  one_turn_map, Id,a0,a_cs
type(real_8) y(6)
type(c_normal_form) normal_form
type(c_taylor) e(3),z_ij,ave_FLOQUET
integer expo(6),pos
character*48 command_gino,fmd
integer i,map_order,j,cas,mf,nd2a,nda,k,t
!!!!!!!!!!!!!!!!!!!!!
type(fibre), pointer :: p1
type(integration_node), pointer :: ti
!!!!!!!!!!!!!!!!!!!!!
c_verbose=.false.
prec=1.d-10 ! for printing
longprint=.false. 
 del=0.d0
use_info=.true.
fmd='(1x,g17.10)'

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

 
 do cas=2,0,-1

if(cas==0) then
 state=nocavity0
 nd2a=4
 nda=2
call kanalnummer(mf,"result_no_cavity.txt")
elseif(cas==1) then
 state=default0
 nd2a=6
 nda=3
call kanalnummer(mf,"result_with_cavity.txt")
elseif(cas==2) then
 state=radiation0
 nd2a=6
 nda=3
call kanalnummer(mf,"result_with_cavity_and_radiation.txt")
endif

map_order=2
call init_all(state,map_order,0)

call alloc(one_turn_map,id,a0,a_cs)
call alloc(y) 
call alloc(normal_form)
call alloc(e)
call alloc(z_ij,ave_FLOQUET)

closed_orbit=0.d0; closed_orbit(5)=del;                              ! (0)

call find_orbit_x(als,closed_orbit(1:6),STATE,1.e-8_dp,fibre1=pos)     ! (1)
del=closed_orbit(5)

id=1   ! map is set to identity                                      ! (2)

! map is added to closed orbit and put into the 6 polymorphs
y(1:6)=closed_orbit(1:6)+id                                          ! (3)

write(mf,*) " ";write(mf,*) "  Constant part of the map";
write(mf,'(a16,6(1x,g12.5))') " closed orbit = ",closed_orbit(1:6)   ! (4)
!write(mf,'(a16,6'//fmd//')') " closed orbit = ",closed_orbit(1:6)   ! (4)

call propagate(als,y(1:6),state,fibre1=pos)                            ! (5)

one_turn_map=y(1:6) ! Six polymorphs are promoted to Taylor maps     ! (6)
closed_orbit=y                                                       ! (7)

one_turn_map=one_turn_map.sub.1                                      ! (8)

write(mf,*) " "; write(mf,*) "  The linear map";
                                   
call print(one_turn_map,mf,prec)                                       


call  c_normal(one_turn_map,normal_form)                             ! (9)


write(mf,'(a8,3(1x,g20.13))') " tune = ",normal_form%tune(1:nda)     ! (10)
if(cas==1) tune_s=normal_form%tune(3)
write(mf,*) " ";write(mf,*) "  Constant part of the map";
 write(mf,'(a16,6(1x,g12.5))') " closed orbit = ",closed_orbit(1:6)  ! (11)

a=normal_form%a_t                                                    ! (12a)
ai=normal_form%a_t**(-1)                                             ! (12b)

mc=from_phasor(-1)*normal_form%a_t**(-1)*one_turn_map*normal_form%a_t*from_phasor()  ! (13)

write(mf,'(/,a37,15x,a41,/)')" i  j    Real(mc(i,j))    Im(mc(i,j))", &
                             "damping  + i *   phase  ( =log(mc(i,j)) )"
do i=1,c_%nd2
do j=1,c_%nd2
if(abs(mc(i,j))>1.e-10_dp) write(mf,'(i2,1x,i2,2(5x,g12.5),5x,2(5x,g12.5))')&
                                      i,j,mc(i,j),log(mc(i,j)) ! (14)
enddo 
enddo
 
!!!!!!! Lattice functions !!!!!!!
!! coefficient of invariant  !!

do i=1,c_%nd
 Ba(1:6,1:6,i)=matmul(matmul(a,Sa(1:6,1:6,i)),ai)  ! (15a)
 Ha(1:6,1:6,i)=matmul(matmul(a,Ia(1:6,1:6,i)),ai)  ! (15b)
 Ka(1:6,1:6,i)=-matmul(S,Ba(1:6,1:6,i))            ! (15c)
 Ea(1:6,1:6,i)=-matmul(Ba(1:6,1:6,i),S)            ! (15d)
enddo
call clean_mat(Ba,prec);call clean_mat(Ha,prec);   ! (16a)
call clean_mat(Ka,prec);call clean_mat(Ea,prec);   ! (16b)

if(cas/=2) then
  do i=1,nda
      expo=0
      expo(2*i-1)=2
     e(i)=1.d0.cmono.expo
      expo=0
      expo(2*i)=2
     e(i)=e(i)+(1.d0.cmono.expo)
     e(i)=e(i)*normal_form%a_t**(-1)               ! (17)
  enddo

  do i=1,nda
write(mf,*)" "; write(mf,*) " Invariant (q^2+p^2) o A^(-1)  in plane ",i; write(mf,*)" "; 
   call print(e(i),mf,prec)
   do j=1,2*nda
    write(mf,'(6'//fmd//')')Ka(j,1:6,i)          ! (18)
   enddo
  enddo
 
  do i=1,2*nda
   do j=i,2*nda
      z_ij=(1.d0.cmono.i)*(1.d0.cmono.j)
      call C_AVERAGE(z_ij,normal_form%a_t,ave_FLOQUET)          ! (19a)
      write(mf,*);write(mf,'(a3,i1,a3,i1,a1)')  "<z_",i," z_",j,">"
      write(mf,'(3(5x,g12.5))') Ea(i,j,1)/2,Ea(i,j,2)/2,Ea(i,j,3)/2
      call print(ave_FLOQUET,mf,prec)
   enddo
  enddo
endif

  do i=1,3
   write(mf,*)" "; write(mf,*) " Matrix H ",i; write(mf,*)" "; 
   do j=1,6
    write(mf,'(6'//fmd//')') Ha(j,1:6,i)
   enddo
  enddo

  if(cas/=0) then
   alphas(cas)=Ba(6,5,3)*sin((normal_form%tune(3))*twopi)
   write(mf,'(a33)') " Approximate time slip  "
   write(mf,'((5x,g12.5))') alphas(cas)
     write(mf,'(a42)') " H based dispersion as in Chao-Sands paper"
   do j=1,6
    write(mf,*) i, Ha(j,5,3)/Ha(5,5,3)
   enddo

  endif

   if(cas==0) then
   write(mf,'(a34)') " Standard time slip without cavity"
   write(mf,'((5x,g12.5))') real(mc(6,5))
   write(mf,'(a34)') " Time slip estimated with a cavity"
   write(mf,'((5x,g12.5))') alphas(1)
   write(mf,'(a48)') " Time slip estimated with a cavity and radiation"
   write(mf,'((5x,g12.5))') alphas(2)


!!!!  Computation by "hand" of the time slip alpha_s
call c_full_canonise(normal_form%A_t,a_cs,a0=a0)   ! (0)
    mat=one_turn_map
    f=a0                                           ! (1)

! not trusting f(6,1:4) from normal form
    f(6,1:4)=0.0_dp 
    do i=1,4
    do j=1,4
     f(6,j)=f(6,j)-f(i,5)*S(i,j)                               ! (2)
    enddo
    enddo

    m(6,1:4) = mat(6,1:4) !(mat(6,1:4)-matmul(f(6,1:4),mat(1:4,1:4)))     ! (3a)
    alpha_s=mat(6,5)
    do i=1,4
     alpha_s =  m(6,i)* f(i,5)+ alpha_s ! - f(6,i)*mat(i,5) + alpha_s    !(3b)
    enddo
    write(mf,*) " Time slip alpha_s from computation and normal form"
    write(mf,*)   alpha_s,real(mc(6,5))
    call move_to(als,p1,"CAV",pos,reset=.true.)
    cav_tune=(alpha_s/2)*(p1%mag%volt*1e-3_dp/p1%mag%P%P0C)*twopi*p1%mag%freq/CLIGHT
    write(mf,*) "tune using alpha_s ", sign(acos(1.0_dp+cav_tune)/2/pi,p1%mag%volt)  
    write(mf,*) "Exact temporal tune", tune_s ; write(6,*) " ";
     write(mf,'(a21)') " Standard dispersion "
    do i=1,6
     write(mf,*) i, a(i,5)      
    enddo
   endif
   if(cas/=0) then
    write(mf,'(a56)') " Almost Exact Crab Angle ignoring transverse emittances "
    write(mf,'((5x,g12.5))') -Ka(5,2,3)/Ka(5,5,3)
    write(mf,'(a56)') " Approximate Crab Angle ignoring transverse emittances  "
    write(mf,'((5x,g12.5))')  Ha(1,6,3)-Ka(5,6,3)*Ha(1,5,3)/Ka(5,5,3)
  endif


!!!!!           Evaluation of the Kinetic Invariants          !!!!!!     
ray(1,1:6)= (/3,4,3,6,1,3/)
ray(2,1:6)= (/3,5,2,6,1,0/)
ray(3,1:6)= (/2,4,3,6,1,1/)
ray(4,1:6)= (/3,5,3,6,1,2/)
ray(5,1:6)= (/2,3,5,8,1,0/)

 mat=one_turn_map

do t=1,10


 do k=1,5
   ray(k,1:6)=matmul(mat,ray(k,1:6))
 enddo

do i=1,6
do j=i,6
 mom(i,j)=0.d0
do k=1,5
 mom(i,j)=ray(k,i)*ray(k,j)/5+mom(i,j)
enddo
enddo
enddo

inv=0.D0
  
do i=1,3
inv=4*(mom(2*i-1, 2*i-1)*mom(2*i, 2*i)-mom(2*i-1, 2*i)**2)+inv
enddo
inv1=inv

inv=inv+8*(mom(1, 3)*mom(2, 4)-mom(1, 4)*mom(2, 3)) &
       +8*(mom(1, 5)*mom(2, 6)-mom(1, 6)*mom(2, 5)) &
       +8*(mom(3, 5)*mom(4, 6)-mom(3, 6)*mom(4, 5))  

write(mf,'(a19,i2,a2,1x,g12.5,a18,1x,g12.5)') &
" Invariant at turn ",t,  " =",inv,' 1-d-f invariant =',inv1
enddo



call kill(one_turn_map,id,a0,a_cs)
call kill(y) 
call kill(normal_form)
call kill(e)
call kill(z_ij,ave_FLOQUET)

close(mf)
 enddo ! cas

call ptc_end(graphics_maybe=1,flat_file=.false.)

end program one_turn_cavity_map

subroutine clean_mat(H,prec)
use precision_constants
implicit none
real(dp) H(6,6,3),prec
integer i,j,k

do i=1,6
do j=1,6
do k=1,3
 if(abs(H(i,j,k))<prec) H(i,j,k)=0
enddo
enddo
enddo

end subroutine clean_mat