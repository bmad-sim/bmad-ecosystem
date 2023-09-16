!+
! Sections refer to the Book:
!   From Tracking Code to Analysis -- Generalised Courant-Snyder theory for any accelerator models
!   Etienne Forest
!
! A copy of this book is included in the directory forest/doc
!-

program spin_phase_advance_isf

use madx_ptc_module
use pointer_lattice

implicit none

type(probe) xs0,xs1,XST
type(probe_8) xs
type(layout), pointer :: als
INTEGER MF,mfisf,I,N,k,pos,no,kp,nturn,mfa
TYPE(FIBRE),POINTER:: P
type(internal_state) state
real(dp)  prec,sig(6),cut,n_isf(3), closed(6), x(6), theta0 
logical break_symmetry
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
logical canonize
!-----------------------------------

prec=1.d-16

! True => Transform lattice elements (sextupoles) into thin lenses (drift-kick-drift) so can do
! Averaging formula of 7.41 exactly.
thin=.false.   

call ptc_ini_no_append        ! Init PTC with empty mad_universes
call append_empty_layout(m_u) ! Append in m_u mad_universe an initial layout
 
als => m_u%start                ! Point to the initial layout
call build_lattice_als0 (ALS,mis,exact=.false.,thin=thin) ! Creates the lattice. This subroutine below. 

p => als%start                  ! Point to initial fibre

! Make the lattice linear by zeroing all sextupoles

courant_snyder_teng_edwards = .true.  ! Ensure linear phase advance agrees with teng-edwards formalism
                                      ! See Sections 7.4 and 7.7.3.2 
time_lie_choice = .true.              ! See Sections 7.3.2 and 7.7.3.1. Only relevant to Jordan normal form. 

! Break symmetry of 1-turn map to make the 1-turn tune more relevant (break the super-symmetry of the lattice.)

do i=1, als%n
  if (p%mag%p%nmul >= 3) then   ! If element has sextupole (or higher) component
    call add (P,2,0,0.1D0)    ! Adds quadrupole component
    call make_it_knob (p%magp%bn(2), 1) ! Make quadrupole component into a knob (DA variable).
    exit  ! Only do this with first sextupole
  endif
  p=>p%next    ! Next element
enddo

! Misalign lattice.
! If we don't misalign and run thin = True, the theoretical formula of 7.41 will agree with simulation.
! We misalign to illustrate the limitations of this formula.

x = [1.0d-5, 1.0d-5, 1.0d-5, 0.0d0, 1.0d-3, 1.0d-6]   ! Misalignment vector
cut = 4.d0                                            ! Gaussian sigma cutoff
call MESS_UP_ALIGNMENT(als, x, cut)                   ! Misalign with Gaussian distribution

!

x = 0.d0

closed = 0.d0

state = nocavity0    ! (2) No cavity but 6D phase space (For Jordan normal form.)

CALL FIND_ORBIT (ALS, CLOSED, 1, STATE, c_1d_5)  ! (3) Find closed orbit
 
state = state + SPIN0    ! Track the spin


write(6, *) " Give Taylor order (1, 2, 3, 4, ...) "
read(5, *) no
write(6, *) " Canonized: true or false "
read(5, *) canonize


call init_all(STATE, no, 1)  ! Increase number of DA variables by 1 to handle quadrupole knob defined above.

call alloc(c_map)
call alloc(c_n)
call alloc(c_spin0)
CALL ALLOC(ISF); call alloc(S_ISF); call alloc(ISFoM);call alloc(dISF);call alloc(O);
call alloc(U, U_c, D, f, A, b, R, id_s, D_tilde)
call alloc(phase);call alloc(xs);
call alloc(nu_spin, fonction, fonction_FLOQUET, a12, a34, h_poisson_bracket)
call alloc(h_vector_field)

if (thin) call kanalnummer(mfa, "analytical_x_average.txt")

!!!! Polymorphic probe is created in the usual manner 
   XS0 = CLOSED    ! Init probe xs0
   ID_S = 1        ! Init c_damap id_s to identity
   XS = XS0+ID_S   ! init probe_8 xs

! get spin polymorphic probe after one turn (1-turn map)
! Propagate 1-turn. The "+" activates tracking with the quad knob. 

CALL propagate(ALS, XS, +STATE, FIBRE1 = 1)  ! (4)

! Copy probe_8 into a complex damap 
c_map = XS ! (5)

! Calculate Normal form c_n including spin.
 
call c_normal(c_map, c_n, dospin = my_true)  ! (6)

! Print:
!   transverse orbital tunes c_n%tune(1:2) 
!   path length slip c_n%tune(3) meters / delta_p/p / turn
!   spin tune. Note: spin tune is negative due to FPP convention.

 write(6, '(4(1x, g21.14))') c_n%tune(1:3),  c_n%spin_tune

if(canonize) then
! spin-orbital similarity transformation

U = c_n%As * c_n%A_t    ! att = c_n%A_t*c_n%As
 
! id_s is a normal form rotation

 call c_full_canonise(U, U_c, D, f, A, b, R, phase, nu_spin) ! (7b)
else
 D=c_n%As
endif

if(canonize) then 
 call kanalnummer(mfisf, "checking_isf_canonized.txt")
else 
 call kanalnummer(mfisf, "checking_isf_noncanonized.txt")
endif

ISF = 2         ! (Fa)   (0,1,0) Normalizing with respect to the y-direction
ISF = D%s * ISF ! (Fb)    

write (mfisf, *) 'Printing ISF as function of DA variables: phase space + quadrupole knob'
write (mfisf, *) 'Analysis is complex but notice that complex part is essentially zero.'

call print(ISF, mfisf)  

write(mfisf, *) "!!!!  Exploring the ISF (invarient spin field) at the end of the lattice !!!!"
Write(mfisf, *) " Testing Barber's Equation S ISF  =  ISF o m "

S_ISF  =  c_map%s * ISF   !  (12a)  Propagating spin over 1-turn using spin matrix.
ISFoM  =  ISF * c_map     !  (12b)  Propagating ISF over 1-turn using phase space map.

Write (mfisf, *) 'Difference between Propagating spin over 1-turn using spin matrix and'
write (mfisf, *) 'Propagating ISF over 1-turn using phase space map:  |S ISF-  ISF o m|/ |S ISF| '
write (mfisf, *) 'Differnece should be zero.'

do i = 1, 3
  dISF%v(i) = S_ISF%v(i)-ISFoM%v(i)
  write(mfisf, *) i, full_abs(dISF%v(i)), full_abs(dISF%v(i))/full_abs(S_ISF%v(i)) ! (12c)
enddo

!

x = 0.d0
x(1) = 0.0002d0 ;x(3) = 0.0002d0 ;    ! (13a)
xs1 = closed+x
 
xst = 0
cray%x = 0.d0
cray%x(1:6) = x
do i = 1, 3
 n_isf(i)  =  ISF%v(i) .o. cray  ! (13b)
enddo

Write(6, *) "  Stroboscopic Average 5000 turns : patience "
nturn = 50000
kp = 5000
call stroboscopic_average(als, xs1, xst, 1, STATE, nturn, kp, ISF_strobo, 6) ! (14c)

Write(mfisf, *) "  Stroboscopic Average "
 
write(mfisf, *); 
write(mfisf, '(a19, 4(1x, g20.13), a19, i8)') " ISF  for x(1:4)  =  " , x(1:4),  " number of turns  =  ",  nturn
write(mfisf, '(a24, 3(1x, g20.13))') " Stroboscopic average ", ISF_strobo
write(mfisf, '(a24, 3(1x, g20.13))') " From the normal form ", n_isf
write(mfisf, '(a4, 20x, 3(1x, g20.13))')" n0 ",   real(ISF%v(1).sub.'0'),  &
real(ISF%v(2).sub.'0'), real(ISF%v(3).sub.'0')
 
close(mfisf)


call ptc_end(graphics_maybe = 1, flat_file = .false.)

contains 

!-----------------------------------------------------------------------
! This subroutine creates a lattice without using Bmad.
! This has several disadvantages so this method of constructing lattices it is strongly discouraged.

subroutine  build_lattice_als0(ALS, mis, error, exact, sl, thin, onecell)
use madx_ptc_module
use pointer_lattice
implicit none

type(layout),  target :: ALS
real(dp), optional :: error(6)
logical,  optional :: exact, sl, thin, onecell
real(dp) :: alpha, lbend,  cut,  ksd,  ksf 
type(fibre)  L1, L2, L3, L4, L5, L6, L7, L8, L9, L10 
type(fibre)  L11, L12, L13, L14, L15, L16, L17, L18, L19, L20, CAVM
type(fibre)  L21, L22, L23, L24, L25, L26, L27, L27A, L27B, L27C, L27D, DS
type(fibre)  QF1, QF2, QD1, QD2, QFA1, QFA2, sf, sd, cav, bend, vc5, bend1 
type(layout) :: sfline, sdline, sup1, supb
logical(lp) :: mis, thi = .false., oneperiod
!-----------------------------------
if (present(thin)) thi = thin

call make_states(.true.)
exact_model  =  .false.;oneperiod  =  .false.
if (present(exact)) exact_model = exact
if (present(onecell)) oneperiod = onecell
call update_states
madlength  =  .false.


call set_mad(energy  =  1.5d0,  method  =  2,  step  =  1)

madkind2  =  matrix_kick_matrix


L1   =  drift("L1 ",   2.832695d0);L2   =  drift("L2 ",   0.45698d0);
L3   =  drift("L3 ",   0.08902d0);L4   =  drift("L4 ",   0.2155d0);
L5   =  drift("L5 ",   0.219d0);L6   =  drift("L6 ",   0.107078d0);
L7   =  drift("L7 ",   0.105716d0);L8   =  drift("L8 ",   0.135904d0);
L9   =  drift("L9 ",   0.2156993d0);L10  =  drift("L10",   0.089084d0);
L11 =  drift("L11",   0.235416d0);L12 =  drift("L12",   0.1245d0);
L13 =  drift("L13",   0.511844d0);L14 =  drift("L14",   0.1788541d0);
L15 =  drift("L15",   0.1788483d0);L16 =  drift("L16",   0.511849d0);
L17 =  drift("L17",   0.1245d0);L18 =  drift("L18",   0.235405d0);
L19 =  drift("L19",   0.089095d0);L20 =  drift("L20",   0.2157007d0);
L21 =  drift("L21",   0.177716d0);L22 =  drift("L22",   0.170981d0);
L23 =  drift("L23",   0.218997d0);L24  =  drift ("L24",   0.215503d0);
L25  =  drift ("L25",   0.0890187d0);L26  =  drift ("L26",   0.45698d0);
L27  =  drift ("L27",   2.832696d0);L27a   =  drift (" L27a",   0.8596d0);
L27b   =  drift (" L27b",   0.1524d0);L27c   =  drift (" L27c",   0.04445d0);
L27d   =  drift (" L27d",   1.776246d0);ds   =  drift (" DS  ",  0.1015d0);

QF1  =  QUADRUPOLE(" QF1 ", 0.344D0,  K1 =  2.2474D0+6.447435260914397D-03)
QF2  =  QUADRUPOLE(" QF2 ", 0.344D0,  K1 =  2.2474D0)
QD1  =  QUADRUPOLE(" QD1 ", 0.187D0,  K1 =  -2.3368D0-2.593018157427161D-02); 
QD2  =  QUADRUPOLE(" QD2 ", 0.187D0,  K1 =  -2.3368D0);  
QFA1 =  QUADRUPOLE(" QFA1", 0.448D0,  K1 =  2.8856D0);  
QFA2 =  QUADRUPOLE(" QFA2", 0.448D0,  K1 =  2.8856D0);  

!!! 1/2 mad-x value
ksf = -41.3355516397069748d0;
ksd = 56.2564709584745489d0;

sf = sextupole ("sf", 2.d0*0.1015d0,  K2 =  ksf);
sd =  sextupole("sd",  2.d0*0.1015d0,  K2 =  ksd);

 VC5 = marker("vc5");
ALPHA = 0.17453292519943295769236907684886d0;
 
LBEND = 0.86621d0;
 
BEND  =  RBEND("BEND",  LBEND,  ANGLE = ALPHA).q.(-0.778741d0)
BEND1  =  RBEND("BEND1",  LBEND,  ANGLE = ALPHA).q.(-0.778741d0)
 
CAVM = MARK("CAVM");
CAV = RFCAVITY("CAV", L = 0.0000d0, VOLT = -1.0d0, REV_FREQ = 500.0d6)

if (thi) then
 sf = sextupole ("sf", 0.d0,  K2 =  ksf*0.203d0);
 sd =  sextupole("sd",  0.d0,  K2 =  ksd*0.203d0);
  sfline = (ds+sf+ds);
  sdline = (ds+sd+ds);
else
 sfline = 1*sf;
 sdline = 1*sd;
endif

SUP1 = L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27+cavm;

SUPb = L1+L2+L3+QF1+VC5+L4+L5+QD1+L6+L7+L8+VC5+BEND+VC5+L9+sfline+L10+&
           L11+QFA1+L12+sdline+L13+ &
           L14+BEND+L15+L16+sdline+L17+ &
           QFA2+L18+L19+sfline+L20+BEND1+L21+&
           L22+QD2+L23+L24+QF2+L25+ &
           L26+VC5+L27+cav;

if (oneperiod) then
 ALS  =  sup1;  !11*sup1+supb;
else
 ALS  =  11*sup1+supb;
endif
if (present(sl)) then
L1   =  drift("L1 ",   2.832695d0);
 if ( sl ) then
  Qf1  =  QUADRUPOLE(" QF1 ", L = 0.d0,  K1 =  0.01d0 ); L1   =  drift("L1 ", L = 0.1d0);
  ALS = L1+QF1;
 endif 
endif

ALS  =  .ring.ALS

call survey(ALS)

if (mis) then
 sig = 1.d-5; cut = 4.d0; 
 if (present(error)) sig = error
 call MESS_UP_ALIGNMENT(ALS, SIG, cut);
endif

end subroutine build_lattice_als0

end program spin_phase_advance_isf


