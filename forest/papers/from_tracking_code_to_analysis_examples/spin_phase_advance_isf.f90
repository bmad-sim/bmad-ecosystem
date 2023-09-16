program spin_phase_advance_isf
use madx_ptc_module
use pointer_lattice
implicit none

type(probe) xs0,xs1,XST
type(probe_8) xs
type(layout), pointer :: als
INTEGER MF,mft,mfisf,I,N,k,pos,no,kp,nturn,mfa
TYPE(FIBRE),POINTER:: P
type(internal_state) state
real(dp)  prec,cut,n_isf(3), closed(6), x(6), theta0 
logical first
logical :: mis=.false.,thin=.false.
type(c_damap) c_map,c_spin0,U,U_c,D,f,A,b,R ,id_s,D_tilde
type(c_taylor) phase(3),nu_spin, fonction,fonction_FLOQUET,a12,a34
type(c_ray) cray
type(c_normal_form) c_n
TYPE(c_spinor) ISF,S_ISF,ISFoM,dISF,O
type(spinor) ISF_strobo   ! real spinor 
type(c_vector_field) h_vector_field
type(c_taylor)  h_poisson_bracket
integer expo(4)
real(dp) DX_AVERAGE_DCS,betax_1,betax_2
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

!-----------------------------------

first=.true.;Lmax = 10.d0;use_info = .true.;prec=1.d-16;thin=.false.

call ptc_ini_no_append
call append_empty_layout(m_u)
ALS=>m_u%start
 
call build_lattice_als(ALS,mis,exact=.false.,thin=thin) 
 
p=>als%start

!!!!!  Make the lattice linear by zeroing all sextupoles !!!!!! 
!!!
courant_snyder_teng_edwards=.true.
time_lie_choice=.true.
p=>als%start
do i=1,als%n
 IF(P%MAG%P%NMUL>=3) THEN
  if(first) then
   CALL ADD(P,2,0,0.1D0)
   call make_it_knob(p%magp%bn(2),1) ! (0)
   first=.false.
  endif
 ENDIF

 p=>p%next
enddo

x(1:3)=1.d-5; x(6)=1.d-6; x(5)=1.d-3; x(4)=0.d0!
cut=4.d0;
call MESS_UP_ALIGNMENT(als,x,cut)  ! (1)
x=0.d0;DX_AVERAGE_DCS=0.d0;
closed=0.d0

state=nocavity0    ! (2)  

CALL FIND_ORBIT(ALS,CLOSED,1,STATE,c_1d_5)  ! (3)
 
state=state+SPIN0

write(6,*) " If order>=4, the phase advance will jump by 100 positions"
write(6,*) " so the fractional tunes will be off"
write(6,*) " but canonisation results will be examined at position 100 "
write(6,*)
write(6,*) " with no=1,2,3 a full Twiss will be done " 
write(6,*) " with the computation <x^2> thrown in if no>1 "
write(6,*);write(6,*);
write(6,*) " Give order no : 1,2,3,4,..."
read(5,*) no


call init_all(STATE,no,1)

call alloc(c_map)
call alloc(c_n)
call alloc(c_spin0)
CALL ALLOC(ISF); call alloc(S_ISF); call alloc(ISFoM);call alloc(dISF);call alloc(O);
call alloc(U,U_c,D,f,A,b,R,id_s,D_tilde)
call alloc(phase);call alloc(xs);
call alloc(nu_spin,fonction,fonction_FLOQUET,a12,a34,h_poisson_bracket)
call alloc(h_vector_field)

call kanalnummer(mft,"spin_twiss.txt")
if(thin) call kanalnummer(mfa,"analytical_x_average.txt")

!!!! Polymorphic probe is created in the usual manner 
   XS0=CLOSED    
   ID_S=1        
   XS=XS0+ID_S 

!!!! get spin polymorphic probe after one turn   
CALL propagate(ALS,XS,+STATE,FIBRE1=1)  ! (4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! Copy probe_8 into a complex damap 
c_map=XS ! (5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
call c_normal(c_map,c_n,dospin=my_true)  ! (6)

 write(6,'(4(1x,g21.14))') c_n%tune(1:3), c_n%spin_tune

U=c_n%As*c_n%A_t    ! att=c_n%A_t*c_n%As
 
! id_s is a rotation
id_s=U**(-1)*c_map*U ! (7a)

! a trick to compute the fractional tunes and time slip
call c_full_canonise(id_s,U_c,D,f,A,b,R,phase,nu_spin) ! (7b)

write(mft,*);write(mft,*) " Fractional tune x"; write(mft,*);

call c_clean_taylor(phase(1),phase(1),prec)
call print(phase(1),mft)

write(mft,*);write(mft,*) " Fractional tune y"; write(mft,*);
call c_clean_taylor(phase(2),phase(2),prec)
call print(phase(2),mft)

if(state%nocavity) then
 write(mft,*);write(mft,*) "  Time "; write(mft,*);
 call c_clean_taylor(phase(3),phase(3),prec)   ! inconsequential error in book
 call print(phase(3),mft)
endif

write(mft,*);write(mft,*) " Fractional tune spin"; write(mft,*);
call c_clean_taylor(nu_spin,nu_spin,prec)
call print(nu_spin,mft)

U=c_n%As*c_n%A_t ! Non-descript U exiting normal form ! (8a)

call c_full_canonise(U,U_c,D,F,A,b,R,phase,nu_spin) ! (8b)

phase(1)=0.d0 ;phase(2)=0.d0 ; ;phase(3)=0.d0; nu_spin=0.d0;

   XS=XS0+U_c   ! (8c)

p => als%start
 
do i=1,als%n

 CALL propagate(ALS,XS,+STATE,FIBRE1=i,fibre2=i+1)  ! (9a)

if((mod(i,100) == 0.or.i==als%n.or.i==1).or.no<=3) then
  if(mod(i,100)==0) write(6,*) " Position ",i 
  
  xs0=xs ! Saving orbit  ! (9b)
  U=XS ! copying in map  ! (9c)

  ! U = U_c o  R = D o f o A o b o R  
  call c_full_canonise(U,U_c,D,F,A,b,R,phase,nu_spin) ! (10)

 if(no>=4) then   ! (A0)  special interlude
   call kanalnummer(mf,"check_canonisation.txt")

   write(mf,*);Write(mf,*) " Time slip factor "; write(mf,*);
   h_vector_field=log(F)
   h_poisson_bracket=getpb(h_vector_field) ! (A1)
   write(mf,*);Write(mf,*) "Fixed point map"; write(mf,*);
   call print(F,mf,prec)
   write(mf,*);Write(mf,*) " Lie exponent of the fixed point map"; write(mf,*);
   call print(h_poisson_bracket,mf,prec)
  
  
   do k=1,4
     id_s%v(k)=F%v(k)-(1.d0.cmono.k) ! (A2)
   enddo
   h_poisson_bracket=h_poisson_bracket &  ! (A3)
  +(id_s%v(1)*(1.d0.cmono.2)-id_s%v(2)*(1.d0.cmono.1)) &
  +(id_s%v(3)*(1.d0.cmono.4)-id_s%v(4)*(1.d0.cmono.3))
   write(mf,*);Write(mf,*) "  Comparing with the canonical form "; write(mf,*);
   call print(h_poisson_bracket,mf,prec)

  expo=0;expo(2)=1;
  a12=A%v(1).par.expo   ! (B1)
  expo=0;expo(4)=1;
  a34=A%v(3).par.expo   ! (B2)
   write(mf,*); Write(mf,*) " Checking Courant-Snyder-Teng-Edwards "
   write(mf,*);Write(mf,*) " A_12 should be zero "; write(mf,*);
   call print(a12,mf,prec)
   write(mf,*);Write(mf,*) " A_34 should be zero "; write(mf,*);
   call print(a34,mf,prec)

   h_vector_field=log(b)    ! (C1)
   h_poisson_bracket=getpb(h_vector_field) ! (C2)
   h_poisson_bracket=h_poisson_bracket*from_phasor()  ! (C3)

   write(mf,*);Write(mf,*) " Lie exponent of the nonlinear part "; write(mf,*);
   call print(h_poisson_bracket,mf,prec)

  !  The original code, in my book, uses U_c=f*A*b 
  !  and D_tilde=to_phasor()*U_c**(-1)*D*U_c*from_phasor().
  !  This messes up the phase advance loop if no>=4
  U=f*A*b     
  D_tilde=to_phasor()*U**(-1)*D*U*from_phasor()  ! (D1)
  O=log(D_tilde%s) ! (D2)

   write(mf,*);Write(mf,*) " Vertical spinor O_y of the canonised D~ "; write(mf,*);
   call print(O%v(2),mf,prec)
 close(mf)
  endif

  !!!!!!!!!!!!!!!   doing something  !!!!!!!!!!!!!!!!!

if(no>1) then
  fonction =2*(1.d0.cmono.1)**2   ! 2*x**2  (Ea) 
  call C_AVERAGE(fonction,U_c,fonction_FLOQUET) ! (Eb)
endif
 write(mft,*) "position, Element ", i, p%mag%name

 betax_1=(U_c%v(1).sub.'1000')**2+(U_c%v(1).sub.'0100')**2 ! (Ec)
 betax_2=(U_c%v(1).sub.'0010')**2+(U_c%v(1).sub.'0001')**2 ! (Ed)
 write(mft,*) " Ripken Beta_x_1 Beta_x_2 ",betax_1,betax_2
 write(mft,*) " 2< x^2 > "
if(no>1) then
  call print(fonction_FLOQUET,mft)
 else
  write(mft,*) "Unfortunately no information if no = 1"
 endif


ISF=2   ! (Fa)
ISF=D%s*ISF ! (Fb)

Write(mft,*) " ISF vector n "
call print(ISF,mft)

write(mft,*) " phase x"
call c_clean_taylor(phase(1),phase(1),prec)
call print(phase(1),mft)

write(mft,*) " phase y"
call c_clean_taylor(phase(2),phase(2),prec)
call print(phase(2),mft)


if(state%nocavity) then
 write(mft,*) "  Time "
 call c_clean_taylor(phase(3),phase(3),prec)
 call print(phase(3),mft)
endif

write(mft,*) " phase spin"
call c_clean_taylor(nu_spin,nu_spin,prec)
call print(nu_spin,mft)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

If(thin) then
 if(p%mag%name(1:2)=="SF".or.p%mag%name(1:2)=="SD") then
  DX_AVERAGE_DCS=(betax_1)**1.5_DP*p%mag%BN(3)/4.0_DP &   ! (Fa) 
    *(-SIN(PHASE(1)*TWOPI)+SIN((PHASE(1)-c_n%TUNE(1))*TWOPI)) &
    /(1.0_DP-COS(c_n%TUNE(1)*TWOPI)) + DX_AVERAGE_DCS
 endif
endif
write(mft,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 

XS=XS0+U_c ! (11)
endif

p=>p%next
enddo

if(thin) then ; 

 DX_AVERAGE_DCS=DX_AVERAGE_DCS*SQRT(betax_1) ! (Fb) 
 write(mfa,'(a11,F20.13,a20)')'  d<x>/dCS ', dx_average_dCS, " < ---- analytical  "

 if(no==2.or.no==3) then

 fonction =(1.d0.cmono.1)  ! x  
  write(mfa,*);write(mfa,*) "Full <x> "
  call C_AVERAGE(fonction,U_c,fonction_FLOQUET) ! (Fc)
  call print(fonction_FLOQUET,mfa)
  write(mfa,*);write(mfa,*) "Full x-dispersion  "
  call print(F%v(1),mfa)
 endif
endif

call kanalnummer(mfisf,"checking_isf.txt")

write(mfisf,*) "!!!!  Exploring the ISF  at the end of the lattice !!!!"

Write(mfisf,*) " Testing Barber's Equation S ISF = ISF o m "

S_ISF = c_map%s*ISF   !  (12a)
ISFoM = ISF*c_map     !  (12b)

Write(mfisf,*) "  |S ISF-  ISF o m|/ |S ISF| "

do i=1,3
dISF%v(i)=S_ISF%v(i)-ISFoM%v(i)
write(mfisf,*) i,full_abs(dISF%v(i)),full_abs(dISF%v(i))/full_abs(S_ISF%v(i)) ! (12c)
enddo

x=0.d0
x(1)=0.001d0 ;x(3)=0.001d0 ;    ! (13a)
xs1=closed+x
 
xst=0
cray%x=0.d0
cray%x(1:6)=x
do i=1,3
 n_isf(i) = ISF%v(i).o.cray  ! (13b)
enddo

Write(6,*) "  Stroboscopic Average 5000 turns : patience "
nturn=5000
kp=1000
call stroboscopic_average(als,xs1,xst,1,STATE,nturn,kp,ISF_strobo,mfisf) ! (14c)

Write(mfisf,*) "  Stroboscopic Average "
 
write(mfisf,*); 
write(mfisf,'(a19,4(1x,g20.13),a19,i4)') " ISF  for x(1:4) = " &
,x(1:4), " number of turns = ", nturn
write(mfisf,'(a24,3(1x,g20.13))') " Stroboscopic average ",ISF_strobo
write(mfisf,'(a24,3(1x,g20.13))') " From the normal form ",n_isf
write(mfisf,'(a4,20x,3(1x,g20.13))')" n0 ",  real(ISF%v(1).sub.'0'), &
real(ISF%v(2).sub.'0'),real(ISF%v(3).sub.'0')
 
close(mfisf)
close(mft)
if(thin) close(mfa)

 call ptc_end(graphics_maybe=1,flat_file=.false.)

end program spin_phase_advance_isf


