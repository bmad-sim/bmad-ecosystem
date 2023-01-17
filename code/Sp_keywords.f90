!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module madx_keywords
  use S_fitting
  implicit none
  public
  logical(lp)::mad8=my_false
  integer :: ifield_name=0
  logical(lp),private :: do_survey =my_false
  logical(lp) :: print_marker =my_true
  type(tree_element), private, allocatable :: t_e(:),t_ax(:),t_ay(:)
  real(dp), private :: a_(3),ent_(3,3), b_(3),exi_(3,3),angcsp,xcsp,dcsp,hcsp,vcsp
  integer, private :: metwig,nstwig 
  logical :: tpsa_pancake=.true.
  logical :: old_name_vorname = .false.
  logical :: readingmaps = .true.
  type keywords
     character*20 magnet
     character*20 model
     logical(lp) FIBRE_flip
     INTEGER FIBRE_DIR
     integer method
     integer nstep
     logical(lp) exact
     logical(lp) madLENGTH
     logical(lp) mad8
     real(dp) tiltd
     type(el_list) LIST
  end type keywords

  type MADX_SURVEY
     REAL(DP) ALPHA,TILT,LD
     REAL(DP) PHI,THETA,PSI
     TYPE(CHART) CHART
  END type MADX_SURVEY

  include "a_namelists.inc"


 ! INTERFACE read_lattice_append
!     MODULE PROCEDURE read_universe_database 
 ! END  INTERFACE


contains


  subroutine create_fibre_append(append,mylat,key,EXCEPTION,magnet_only,br,bri)  
    implicit none

!    type(mad_universe), target, intent(inout)  :: m_u
    type(layout), target, intent(inout)  :: mylat
    logical(lp), optional :: magnet_only
    type(keywords) key
    INTEGER EXCEPTION  !,NSTD0,METD0
    logical(lp) doneit,append
    type(fibre), pointer :: current
    type (taylor),optional, INTENT(INout):: br(:,:)
    integer,optional, INTENT(INout):: bri(:,:)
    call set_metc_for_pancake(key%method)
    if(append) then
     call append_empty(mylat)
    else
     if(associated(mylat%end)) then
      IF(ASSOCIATED(mylat%T)) THEN
         CALL kill_NODE_LAYOUT(mylat%T)  !  KILLING THIN LAYOUT
         nullify(mylat%T)
        if(lielib_print(12)==1) WRITE(6,*) " NODE LAYOUT HAS BEEN KILLED "
       ENDIF      
        mylat%end=-1
       else
        call append_empty(mylat)
     endif
    endif
     call  create_fibre(mylat%end,key,EXCEPTION,magnet_only,br,bri)
     
    if(.not.append) then
     mylat%closed=my_true

     doneit=my_true
     call ring_l(mylat,doneit)
     call line_l(mylat,doneit)
     mylat%closed=.false.
   !  call survey(mylat)
     call MAKE_NODE_LAYOUT( mylat)     
    endif
  end subroutine create_fibre_append
   
  subroutine change_method_in_create_fibre(ptc_key,nterm,change)
   implicit none
   type(keywords) ptc_key
   integer kind00,met,nst,nterm
   logical change

   kind00=0
   if(ptc_key%magnet=='wiggler') kind00=kindwiggler
   if(ptc_key%magnet=='INTERNALPANCAKE') kind00=kindpa
   if(ptc_key%magnet=='PANCAKEBMAD    ') kind00=kindpa
   if(ptc_key%magnet=='PANCAKEBMADZERO') kind00=kindpa
   if(ptc_key%magnet=='wiggler') then 
     limit_int0_new=limit_int0_new*nterm
       call against_the_method(ptc_key%method,ptc_key%nstep,met,nst,kind00,change)
     limit_int0_new=limit_int0_new/nterm
   else
    call against_the_method(ptc_key%method,ptc_key%nstep,met,nst,kind00,change)
    if(lielib_print(17)==1.and.change) then
     write(6,*) ptc_key%LIST%NAME, "recut ",met,nst," to ",ptc_key%method,ptc_key%nstep
    endif
    if(switch_to_drift_kick.and.change) then
      ptc_key%model = 'DRIFT_KICK'
      if(lielib_print(17)==1)write(6,*) "also changed to drift-kick-drift "
    endif 

  endif
  end subroutine change_method_in_create_fibre

  subroutine create_fibre(el,key,EXCEPTION,magnet_only,br,bri)
    implicit none
    integer ipause, mypause,i
    type(fibre), target, intent(inout)::el
    logical(lp), optional :: magnet_only
    type(keywords) key
    type(el_list) blank
    character*255 magnet
    character*17 MODEL
    INTEGER EXCEPTION  !,NSTD0,METD0
    LOGICAL(LP) EXACT0,magnet0
    logical(lp) FIBRE_flip0,MAD0
    logical(lp) :: t=my_true,f=my_false
    INTEGER FIBRE_DIR0,IL
    real(dp) e1_true,norm
    type (taylor),optional, INTENT(INout):: br(:,:)
    integer,optional, INTENT(INout):: bri(:,:)

    IL=15

    if(present(magnet_only)) then
       magnet0=magnet_only
    else
       magnet0=my_false
    endif

    blank=0
    magnet=key%magnet
    call context(magnet)
    model=key%model
    call context(model)

    CALL SET_MADX_(t,magnet0)


    select case(MODEL)
    CASE("DRIFT_KICK       ")
       MADTHICK=drift_kick_drift
    CASE("MATRIX_KICK      ")
       MADTHICK=matrix_kick_matrix
    CASE("DELTA_MATRIX_KICK")
       MADTHICK=kick_sixtrack_kick
    CASE DEFAULT
 
       EXCEPTION=1
       ipause=mypause(444)
       RETURN
    END SELECT
    !   MADTHICK=drift_kick_drift
    !    NSTD0=NSTD
    !    METD0=METD
    EXACT0=EXACT_MODEL
    FIBRE_FLIP0= FIBRE_FLIP
    FIBRE_DIR0=FIBRE_DIR
    MAD0=MAD

    KEY%LIST%nst=KEY%NSTEP
    KEY%LIST%method=KEY%METHOD
    EXACT_MODEL=KEY%EXACT
    FIBRE_FLIP = KEY%FIBRE_FLIP
    FIBRE_DIR  = KEY%FIBRE_DIR
    MADLENGTH=KEY%MADLENGTH

    !     real(dp) L,LD,LC,K(NMAX),KS(NMAX)
    !     real(dp) ang(3),t(3)
    !     real(dp) angi(3),ti(3)
    !     integer patchg
    !     real(dp) T1,T2,B0
    !     real(dp) volt,freq0,harmon,lag,DELTA_E,BSOL
    !     real(dp) tilt
    !     real(dp) FINT,hgap,h1,h2,X_COL,Y_COL
    !     real(dp) thin_h_foc,thin_v_foc,thin_h_angle,thin_v_angle  ! highly illegal additions by frs
    !     CHARACTER(120) file
    !     CHARACTER(120) file_rev
    !    CHARACTER(nlp) NAME
    !     CHARACTER(vp) VORNAME
    !     INTEGER KIND,nmul,nst,method
    !     LOGICAL(LP) APERTURE_ON
    !     INTEGER APERTURE_KIND
    !     REAL(DP) APERTURE_R(2),APERTURE_X,APERTURE_Y
    !     LOGICAL(LP) KILL_ENT_FRINGE,KILL_EXI_FRINGE,BEND_FRINGE,PERMFRINGE
    !     REAL(DP) DPHAS,PSI,dvds
    !     INTEGER N_BESSEL

    if(sixtrack_compatible) then
       EXACT_MODEL=my_false
       KEY%LIST%method=2
       MADTHICK=drift_kick_drift
    endif


    SELECT CASE(magnet(1:IL))
    CASE("DRIFT          ")
       BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("SUPERDRIFT     ")
       BLANK=SUPERDRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("SOLENOID       ")
       if(sixtrack_compatible) stop 1
       if(KEY%LIST%L/=0.0_dp) then
          BLANK=SOLENOID(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
          !  BLANK%bend_fringe=key%list%bend_fringe
       else
          write(6,*) "switch solenoid to dubious thin multipole "
          BLANK=MULTIPOLE_BLOCK(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       endif

    CASE("THICKMULTIPOLE ")
       BLANK=multipoleTILT(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe

    CASE("QUADRUPOLE     ")
       BLANK=QUADRUPOLE(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe
    CASE("SEXTUPOLE     ")
       BLANK=SEXTUPOLE(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe
    CASE("OCTUPOLE      ")
       BLANK=OCTUPOLE(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe
    CASE("SBEND         ")
       BLANK=SBEND(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("TRUERBEND     ")
       if(sixtrack_compatible) stop 2

  !     e1_true= KEY%LIST%b0/2.0_dp+ KEY%LIST%t1
       BLANK=rbend(KEY%LIST%NAME,l=KEY%LIST%l,angle=KEY%LIST%b0,list=KEY%LIST)
 !      BLANK=rbend(KEY%LIST%NAME,l=KEY%LIST%l,angle=KEY%LIST%b0,e1=e1_true,list=KEY%LIST)

    CASE("WEDGRBEND     ")
       if(sixtrack_compatible) stop 3

       BLANK=rbend(KEY%LIST%NAME,l=KEY%LIST%l,angle=KEY%LIST%b0,e1=KEY%LIST%t1,e2=KEY%LIST%t2,list=KEY%LIST)

    CASE("RBEND         ")
       if(sixtrack_compatible) stop 4
       KEY%LIST%T1=KEY%LIST%T1+KEY%LIST%B0/2.0_dp
       KEY%LIST%T2=KEY%LIST%T2+KEY%LIST%B0/2.0_dp
       BLANK=SBEND(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("KICKER         ","VKICKER        ","HKICKER        ")
       BLANK=KICKER(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       BLANK%bend_fringe=key%list%bend_fringe
    CASE("MONITOR        ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       endif
    CASE("HMONITOR        ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST) ;BLANK%KIND=KIND12;
       endif
    CASE("VMONITOR       ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST) ;BLANK%KIND=KIND13;
       endif
    CASE("INSTRUMENT     ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=MONITOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST) ;BLANK%KIND=KIND14;
       endif
    CASE("MARKER         ")
       BLANK=MARKER(KEY%LIST%NAME,list=KEY%LIST)
    CASE("CHANGEREF      ")
       if(sixtrack_compatible) stop 5
       BLANK=CHANGEREF(KEY%LIST%NAME,KEY%LIST%ANG,KEY%LIST%T,KEY%LIST%PATCHG)
    CASE("RFCAVITY       ")
       if(sixtrack_compatible) then
          If(KEY%LIST%L/=0.0_dp) stop 60
          If(KEY%LIST%N_BESSEL/=0.0_dp) stop 61
          norm=0.0_dp
          do i=1,nmax
             norm=norm+abs(KEY%LIST%k(i))+abs(KEY%LIST%ks(i))
          enddo
          norm=norm-abs(KEY%LIST%k(2))
          if(norm/=0.0_dp) then
             write(6,*) norm
             stop 62
          endif
       endif
       BLANK=RFCAVITY(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("TWCAVITY       ")
       if(sixtrack_compatible) stop 7
       BLANK=TWCAVITY(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("ELSEPARATOR    ")
       if(sixtrack_compatible) stop 8
       BLANK=ELSEPARATOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("MULTIPOLE_BLOCK","MULTIPOLE      ")
       BLANK=MULTIPOLE_BLOCK(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("SMI            ","SINGLE_LENS    ")
       if(sixtrack_compatible) stop 9
       BLANK=SINGLE_LENS(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("RCOLLIMATOR    ")
       if(sixtrack_compatible) stop 10
       BLANK=RCOLLIMATOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("ECOLLIMATOR    ")
       if(sixtrack_compatible) then
          BLANK=DRIFT(KEY%LIST%NAME,LIST=KEY%LIST)
       else
          BLANK=ECOLLIMATOR(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
       endif
    CASE("WIGGLER        ")
       if(sixtrack_compatible) stop 12
       BLANK=WIGGLER(KEY%LIST%NAME,t=tilt.is.KEY%tiltd,LIST=KEY%LIST)
    CASE("HELICALDIPOLE  ")
       if(sixtrack_compatible) stop 13
       BLANK=HELICAL(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("PANCAKE        ")
       if(sixtrack_compatible) stop 13
       BLANK=pancake(KEY%LIST%NAME,KEY%LIST%file)
    CASE("ABELL_DRAGT    ")
       if(sixtrack_compatible) stop 13
       BLANK=abell_dragt(KEY%LIST%NAME,LIST=KEY%LIST)
    CASE("INTERNALPANCAKE")
       if(sixtrack_compatible) stop 13
       BLANK=pancake(KEY%LIST%NAME,br=br)
    CASE("PANCAKEBMAD")
       if(sixtrack_compatible) stop 13
       BLANK=pancake_bmad(KEY%LIST%NAME,br=bri)
    CASE("PANCAKEBMADZERO")
       if(sixtrack_compatible) stop 13
       BLANK=pancake_bmad_empty(KEY%LIST%NAME)
    CASE DEFAULT
       WRITE(6,*) " "
       WRITE(6,*) " THE MAGNET"
       WRITE(6,*) " "
       WRITE(6,*) "  --->   ",MAGNET(1:IL)
       WRITE(6,*) " "
       WRITE(6,*)  " IS NOT PERMITTED "
       STOP 666
    END SELECT

    BLANK%VORNAME = KEY%LIST%VORNAME
    CALL EL_Q_FOR_MADX(EL,BLANK)
    !  added 2007.07.09
    el%mag%parent_fibre =>el
    el%magp%parent_fibre=>el
    !  end of added 2007.07.09


    CALL SET_MADX_(f,f)

    !    NSTD=NSTD0
    !    METD=METD0
    EXACT_MODEL=EXACT0
    FIBRE_FLIP= FIBRE_FLIP0
    FIBRE_DIR=FIBRE_DIR0
    MAD=MAD0

    !    IF(ASSOCIATED(EL%PREVIOUS)) THEN
    !     if(.not.associated(EL%POS))allocate(EL%POS)
    !     EL%POS=EL%PREVIOUS%POS+1
    !    ELSE
    !     if(.not.associated(EL%POS))allocate(EL%POS)
    !     EL%POS=1
    !    ENDIF
    if(key%list%BEND_FRINGE) then
       el%mag%p%bend_fringe=my_true
       el%magp%p%bend_fringe=my_true
    endif

    if(el%mag%kind==kind4) then
       el%mag%c4%CAVITY_TOTALPATH=key%list%CAVITY_TOTALPATH
       el%magp%c4%CAVITY_TOTALPATH=key%list%CAVITY_TOTALPATH
    endif

  end subroutine create_fibre

  subroutine zero_key(key)
    implicit none

    type(keywords) , intent(out):: key
    key%magnet="CROTTE"
    select case(MADTHICK)
    CASE(drift_kick_drift)
       key%model="DRIFT_KICK       "
    CASE(matrix_kick_matrix)
       key%model="MATRIX_KICK      "
    CASE(kick_sixtrack_kick)
       key%model="DELTA_MATRIX_KICK"
    END SELECT


    key%FIBRE_flip=FIBRE_flip
    key%FIBRE_DIR=FIBRE_DIR
    key%method=METD
    key%nstep=NSTD
    key%exact=EXACT_MODEL
    key%madLENGTH=madLENGTH
    key%LIST%NMUL = 1
    key%mad8 = mad8
    key%tiltd=0.0_dp
    key%LIST=0

  end subroutine zero_key

  !  PRINTING FIBRES FOR FLAT FILES
 

 
 

  subroutine print_pancake_field(el,filename)
    implicit none
    type(pancake), pointer :: el
    integer mf,nst,i,j
    character(*) filename
    real(dp) brho,cl
    type(real_8) b(8)
    type(taylor) bt(8)

    call kanalnummer(mf)
    open(unit=mf,file=filename) !,recl=200)
    if(el%p%method==4.or.el%p%method==2) then
     nst=2*el%p%nst+1
    else
     nst=7*el%p%nst+1
    endif
    cl=(clight/1e8_dp)
    BRHO=el%p%p0c*10.0_dp/cl


    call init(EL%B(1)%no,2)
    CALL ALLOC(B)
    CALL ALLOC(Bt)

!    write(mf,*) nst,el%p%ld,el%p%b0,EL%B(1)%no,my_false
    write(mf,*)  el%p%ld,el%p%b0
    write(mf,*) nst,el%p%method,EL%B(1)%no
    write(mf,*) el%p%lc,EL%hc
    write(mf,*) el%dc,el%vc,el%xc  
    write(mf,*) el%angc

!  type  tree_element   !@1  USED FOR FAST TRACKING IN O_TREE_ELEMENT.F90
!     real(dp) ,  DIMENSION(:), POINTER :: CC
!     real(dp) ,  DIMENSION(:), POINTER :: fixr,fix,fix0
!     integer,  DIMENSION(:), POINTER :: JL,JV
!     INTEGER,POINTER :: N,NP,no
!     real(dp), pointer :: e_ij(:,:)
!     real(dp), pointer :: rad(:,:)
!     real(dp), pointer :: ds,beta0,eps
!     logical, pointer :: symptrack,usenonsymp,factored
!  end  type tree_element

 
    do i=1,nst
       B(1)=morph(1.0_dp.mono.1)
       B(2)=morph(1.0_dp.mono.2)
       B(3)=0.0_dp;
       B(4)=0.0_dp;
       B(5)=0.0_dp;
       B(6)=0.0_dp;
       B(7)=0.0_dp;
       B(8)=0.0_dp;
       CALL trackg(EL%B(i),B)
       write(mf,*) i
       do j=1,3
          bt(j)=b(j)*brho
          call print(bt(j),mf)
       enddo

    enddo
    CALL kill(B)
    CALL kill(Bt)
   close(mf)

  end subroutine print_pancake_field

 
  subroutine print_pancake_field_integer(el,filename)
    use dabnew_pancake
    implicit none
    type(pancake), pointer :: el
    integer mf,nst,i,j
    character(*) filename
    real(dp) brho,cl
    integer b(8)

    call kanalnummer(mf)
    open(unit=mf,file=filename) !,recl=200)
 
    if(el%p%method==4.or.el%p%method==2) then
     nst=2*el%p%nst+1
    else
     nst=7*el%p%nst+1
    endif

    cl=(clight/1e8_dp)
    BRHO=el%p%p0c*10.0_dp/cl


!    call init(EL%B(1)%no,2)
    call init_pancake(EL%B(1)%no,2)

  !  CALL ALLOC(B)
    CALL alloc_pancake(B)

!    write(mf,*) nst,el%p%ld,el%p%b0,EL%B(1)%no,my_false
    write(mf,*)  el%p%ld,el%p%b0
    write(mf,*) nst,el%p%method,EL%B(1)%no
    write(mf,*) el%p%lc,EL%hc
    write(mf,*) el%dc,el%vc,el%xc  
    write(mf,*) el%angc

!  type  tree_element   !@1  USED FOR FAST TRACKING IN O_TREE_ELEMENT.F90
!     real(dp) ,  DIMENSION(:), POINTER :: CC
!     real(dp) ,  DIMENSION(:), POINTER :: fixr,fix,fix0
!     integer,  DIMENSION(:), POINTER :: JL,JV
!     INTEGER,POINTER :: N,NP,no
!     real(dp), pointer :: e_ij(:,:)
!     real(dp), pointer :: rad(:,:)
!     real(dp), pointer :: ds,beta0,eps
!     logical, pointer :: symptrack,usenonsymp,factored
!  end  type tree_element

 
    do i=1,nst
  do j=1,8
    call dacon_pancake(b(j),0.0_dp)
   enddo
 !      B(1)=morph(1.0_dp.mono.1)
 !     B(2)=morph(1.0_dp.mono.2)
 !      B(3)=0.0_dp;
 !      B(4)=0.0_dp;
 !!      B(5)=0.0_dp;
 !      B(6)=0.0_dp;
 !      B(7)=0.0_dp;
 !      B(8)=0.0_dp;
 call davar_pancake(B(1),0.0_dp,1)
 call davar_pancake(B(2),0.0_dp,2)
 
       CALL track_g_pancake(EL%B(i),B)
       write(mf,*) i

       do j=1,3
!          b(j)=b(j)*brho
          call dacmu_pancake(b(j),brho,b(j))
          call dapri77_pancake(b(j),mf)
       enddo

    enddo
    CALL kill_pancake(B)
   close(mf)

  end subroutine print_pancake_field_integer
 




  subroutine read_pancake_new(el,filename,t_em,met)
    implicit none
    type(ELEMENT), pointer :: el
    integer mf,nstc,ORDER,I,ii,met
    real(dp)  L,hc,cl,BRHO,ld,xc,dc,vc,hd,lc,angc,ANGLE,ds,a
    logical(lp) REPEAT
    character(*) filename
    type(tree_element), allocatable :: t_em(:)
    TYPE(TAYLOR) B(nbe),ba(nbe),bf(nbe),bn(nbe),it  !,ax(2),ay(2)

    cl=(clight/1e8_dp)
    BRHO=el%p%p0c*10.0_dp/cl




    a=0.0_dp
   ! file_fitted=file


    call kanalnummer(mf)
    open(unit=mf,file=filename) !,recl=200)


 
    read(mf,*) LD,hD  !,REPEAT   ! L and Hc are geometric
    read(mf,*) nstc,met, ORDER 
    read(mf,*) LC,hc
    read(mf,*) dc,vc,xc
    read(mf,*) angc

       angcsp=angc
       xcsp=xc
       dcsp=dc
       hcsp=hc
       vcsp=vc

    ds=LC/nstc
    ii=0
!    if(present(no)) order=no
 


    order=order+1   
 CALL INIT(ORDER,1,0,0)

    CALL ALLOC(B)
    CALL ALLOC(Bf)
    CALL ALLOC(Ba)
    CALL ALLOC(Bn)
    call alloc(it) 
bf(1)=0.0_dp;bf(2)=0.0_dp;bf(3)=0.0_dp;
ba(1)=0.0_dp;ba(2)=0.0_dp;ba(3)=0.0_dp;


!    IF(REPEAT.AND.NST==0) NST=NSTD

    ALLOCATE(t_em(NSTc))  

    read(mf,*) ii 
          CALL READ(Bf(1),mf);CALL READ(Bf(2),mf);CALL READ(Bf(3),mf);

          Bf(1)=Bf(1)/BRHO
          Bf(2)=Bf(2)/BRHO
          Bf(3)=Bf(3)/BRHO
    read(mf,*) ii 
          CALL READ(Ba(1),mf);CALL READ(Ba(2),mf);CALL READ(Ba(3),mf);


          Ba(1)=Ba(1)/BRHO
          Ba(2)=Ba(2)/BRHO
          Ba(3)=Ba(3)/BRHO

    DO I=3,NSTc 



    read(mf,*) ii 
          CALL READ(B(1),mf);CALL READ(B(2),mf);CALL READ(B(3),mf);

 
          B(1)=B(1)/BRHO
          B(2)=B(2)/BRHO
          B(3)=B(3)/BRHO

         if(i==3) then
          Bn(1)=Bf(1)
          Bn(2)=Bf(2)
          Bn(3)=Bf(3)
          bn(4)=-(bn(3).i.2)  ! ax
          it=1.0_dp+hc*(1.0_dp.mono.1)
          bn(5)=(4*b(4)-ba(4)-3*bf(4))/ds/2-it*bn(2)   !  d/dx (1+hx)A_s 
          bn(6)= it*bn(1)   !  d/dy (1+hx)A_s  
          bn(7)=bn(4).d.1   !  d/dx Ax
          bn(8)=bn(4).d.2   !  d/dy Ax   
          CALL SET_TREE_g(t_em(1),Bn)
         elseif(i==nstc) then
          Bn(1)=B(1)
          Bn(2)=B(2)
          Bn(3)=B(3)
          bn(4)=-(bn(3).i.2)  ! ax
          it=1.0_dp+hc*(1.0_dp.mono.1)
          bn(5)=(3*b(4)+bf(4)-4*ba(4))/ds/2-it*bn(2)   !  d/dx (1+hx)A_s 
          bn(6)= it*bn(1)   !  d/dy (1+hx)A_s  
          bn(7)=bn(4).d.1   !  d/dx Ax
          bn(8)=bn(4).d.2   !  d/dy Ax   
          CALL SET_TREE_g(t_em(i),Bn)
         endif

          Bn(1)=Ba(1)
          Bn(2)=Ba(2)
          Bn(3)=Ba(3)
          bn(4)=-(bn(3).i.2)  ! ax
          it=1.0_dp+hc*(1.0_dp.mono.1)
          bn(5)=(b(4)-bf(4))/ds/2-it*bn(2)   !  d/dx (1+hx)A_s 
          bn(6)= it*bn(1)   !  d/dy (1+hx)A_s  
          bn(7)=bn(4).d.1   !  d/dx Ax
          bn(8)=bn(4).d.2   !  d/dy Ax   
          CALL SET_TREE_g(t_em(i-1),Bn)
 
          Bf(1)=Ba(1)
          Bf(2)=Ba(2)
          Bf(3)=Ba(3)  
          Bf(4)=Ba(4)
          Bf(5)=Ba(5)
          Bf(6)=Ba(6)
          Bf(7)=Ba(7)
          Bf(8)=Ba(8)    
  
          Ba(1)=B(1)
          Ba(2)=B(2)
          Ba(3)=B(3)
          Ba(4)=B(4)
          Ba(5)=B(5)
          Ba(6)=B(6)
          Ba(7)=B(7)
          Ba(8)=B(8)

    enddo
    call KILL(B)
    call KILL(Bf)
    call KILL(Ba)
    call KILL(Bn)
    call KILL(it) 


 
 !  else
 !    NST=size(t_em)
 !  endif
    ANGLE=LD*HD


    !    IF(ANG/=zero.AND.R/=zero) THEN
    if(hc/=0.0_dp) then
       el%p%LC=2.0_dp*SIN(ANGLE/2.0_dp)/hD
    else
       el%p%LC=LD
    endif
!    el%p%B0=hD                     wrong in fibre
!    el%p%LD=LD
!    el%p%L=lc


    close(mf)

  end subroutine read_pancake_new



  subroutine read_pancake_new_integer(el,filename,t_em,met)
    use dabnew_pancake
    implicit none
    type(ELEMENT), pointer :: el
    integer mf,nstc,ORDER,I,ii,k,met
    real(dp)  L,hc,cl,BRHO,ld,xc,dc,vc,hd,lc,angc,ANGLE,ds,a,brhoi
    logical(lp) REPEAT
    character(*) filename
    type(tree_element), allocatable :: t_em(:)
    integer B(nbe),ba(nbe),bf(nbe),bn(nbe)   !,it   

    cl=(clight/1e8_dp)
    BRHO=el%p%p0c*10.0_dp/cl




    a=0.0_dp
   ! file_fitted=file


    call kanalnummer(mf)
    open(unit=mf,file=filename) !,recl=200)


 
    read(mf,*) LD,hD  !,REPEAT   ! L and Hc are geometric
    read(mf,*) nstc,met, ORDER 
    read(mf,*) LC,hc
    read(mf,*) dc,vc,xc
    read(mf,*) angc

       angcsp=angc
       xcsp=xc
       dcsp=dc
       hcsp=hc
       vcsp=vc

    ds=LC/nstc
    ii=0
!    if(present(no)) order=no
 
brhoi=1.0_dp/brho

    order=order+1   
 !CALL INIT(ORDER,1,0,0)
 CALL INIT_pancake(ORDER,2)

    CALL ALLOC_pancake(B)
    CALL ALLOC_pancake(Bf)
    CALL ALLOC_pancake(Ba)
    CALL ALLOC_pancake(Bn)
 !   call daall0_pancake(it) 
!bf(1)=0.0_dp;bf(2)=0.0_dp;bf(3)=0.0_dp;
!ba(1)=0.0_dp;ba(2)=0.0_dp;ba(3)=0.0_dp;
  do i=1,nbe
    call dacon_pancake(b(i),0.0_dp)
    call dacon_pancake(bf(i),0.0_dp)
    call dacon_pancake(ba(i),0.0_dp)
    call dacon_pancake(bn(i),0.0_dp)
   enddo
 !   call dacon_pancake(it,0.0_dp)
!    IF(REPEAT.AND.NST==0) NST=NSTD

    ALLOCATE(t_em(NSTc))  

    read(mf,*) ii 

!          CALL READ(Bf(1),mf);CALL READ(Bf(2),mf);CALL READ(Bf(3),mf);
         CALL darea77_pancake(Bf(1),mf);CALL darea77_pancake(Bf(2),mf);CALL darea77_pancake(Bf(3),mf);

    do k=1,3
     call dacmu_pancake(Bf(k),brhoi,bf(k))

    enddo

 

 !         Bf(1)=Bf(1)/BRHO
 !         Bf(2)=Bf(2)/BRHO
 !         Bf(3)=Bf(3)/BRHO
    read(mf,*) ii 
 !         CALL READ(Ba(1),mf);CALL READ(Ba(2),mf);CALL READ(Ba(3),mf);
          CALL darea77_pancake(Ba(1),mf);CALL darea77_pancake(Ba(2),mf);CALL darea77_pancake(Ba(3),mf);

    do k=1,3
     call dacmu_pancake(Ba(k),brhoi,ba(k))
    enddo

     !     Ba(1)=Ba(1)/BRHO
     !     Ba(2)=Ba(2)/BRHO
     !     Ba(3)=Ba(3)/BRHO

    DO I=3,NSTc 



    read(mf,*) ii 
!          CALL READ(B(1),mf);CALL READ(B(2),mf);CALL READ(B(3),mf);
          CALL darea77_pancake(B(1),mf);CALL darea77_pancake(B(2),mf);CALL darea77_pancake(B(3),mf);
     do k=1,3
     call dacmu_pancake(B(k),brhoi,b(k))
    enddo

    !      B(1)=B(1)/BRHO
    !     B(2)=B(2)/BRHO
    !      B(3)=B(3)/BRHO

         if(i==3) then
     do k=1,3
     call dacop_pancake(Bf(k),bn(k))
    enddo

    !      Bn(1)=Bf(1)
    !      Bn(2)=Bf(2)
    !      Bn(3)=Bf(3)
    
! not implemented
   !      bn(4)=-(bn(3).i.2)  ! ax
   !       it=1.0_dp+hc*(1.0_dp.mono.1)
   !       bn(5)=(4*b(4)-ba(4)-3*bf(4))/ds/2-it*bn(2)   !  d/dx (1+hx)A_s 
   !       bn(6)= it*bn(1)   !  d/dy (1+hx)A_s  
   !       bn(7)=bn(4).d.1   !  d/dx Ax
   !       bn(8)=bn(4).d.2   !  d/dy Ax   
          CALL SET_TREE_G_pancake(t_em(1),Bn)
         elseif(i==nstc) then
     do k=1,3
     call dacop_pancake(B(k),bn(k))
    enddo
      !    Bn(1)=B(1)
      !    Bn(2)=B(2)
      !    Bn(3)=B(3)
! not implemented

     !     bn(4)=-(bn(3).i.2)  ! ax
     !     it=1.0_dp+hc*(1.0_dp.mono.1)
     !     bn(5)=(3*b(4)+bf(4)-4*ba(4))/ds/2-it*bn(2)   !  d/dx (1+hx)A_s 
     !     bn(6)= it*bn(1)   !  d/dy (1+hx)A_s  
     !     bn(7)=bn(4).d.1   !  d/dx Ax
     !     bn(8)=bn(4).d.2   !  d/dy Ax   
          CALL SET_TREE_G_pancake(t_em(i),Bn)
         endif
     do k=1,3
     call dacop_pancake(Ba(k),bn(k))
    enddo

   !       Bn(1)=Ba(1)
   !       Bn(2)=Ba(2)
   !       Bn(3)=Ba(3)

! not implemented

       !  bn(4)=-(bn(3).i.2)  ! ax
       !   it=1.0_dp+hc*(1.0_dp.mono.1)
       !   bn(5)=(b(4)-bf(4))/ds/2-it*bn(2)   !  d/dx (1+hx)A_s 
       !   bn(6)= it*bn(1)   !  d/dy (1+hx)A_s  
       !   bn(7)=bn(4).d.1   !  d/dx Ax
       !   bn(8)=bn(4).d.2   !  d/dy Ax   
          CALL SET_TREE_G_pancake(t_em(i-1),Bn)
 
     do k=1,8
     call dacop_pancake(Ba(k),bf(k))
    enddo
     do k=1,8
     call dacop_pancake(B(k),ba(k))
    enddo
    !      Bf(1)=Ba(1)
    !      Bf(2)=Ba(2)
    !      Bf(3)=Ba(3)  
    !      Bf(4)=Ba(4)
    !      Bf(5)=Ba(5)
    !      Bf(6)=Ba(6)
    !      Bf(7)=Ba(7)
    !      Bf(8)=Ba(8)    
  
    !      Ba(1)=B(1)
    !      Ba(2)=B(2)
    !      Ba(3)=B(3)
    !      Ba(4)=B(4)
    !      Ba(5)=B(5)
    !      Ba(6)=B(6)
    !      Ba(7)=B(7)
    !      Ba(8)=B(8)

    enddo
    call kill_pancake(B)
    call kill_pancake(Bf)
    call kill_pancake(Ba)
    call kill_pancake(Bn)
 !   call dadal1_pancake(it) 


 
 !  else
 !    NST=size(t_em)
 !  endif
    ANGLE=LD*HD


    !    IF(ANG/=zero.AND.R/=zero) THEN
    if(hc/=0.0_dp) then
       el%p%LC=2.0_dp*SIN(ANGLE/2.0_dp)/hD
    else
       el%p%LC=LD
    endif
!    el%p%B0=hD                     wrong in fibre
!    el%p%LD=LD
!    el%p%L=lc


    close(mf)

  end subroutine read_pancake_new_integer

  SUBROUTINE track_g_pancake(T,XI)
    use dabnew_pancake
    IMPLICIT NONE
    TYPE(TREE_ELEMENT), INTENT(IN) :: T
    integer, INTENT(INOUT) :: XI(:)
    integer XT(lno),XF(lnv),XM(lno+1),XX
    INTEGER JC,I,IV,itemp

    CALL ALLOC_pancake(XT)
    CALL ALLOC_pancake(XF)
    CALL ALLOC_pancake(XM)
    CALL daall0_pancake(XX)
    CALL daall0_pancake(itemp)
 


    do i=1,T%np
 !      xt(i)=xi(i)
       call dacop_pancake(xi(i),xt(i))
    enddo
    do i=1,T%np
  !     xf(i) = T%cc(i)
       call dacon_pancake(xf(i),T%cc(i))
    enddo
    call dacon_pancake(XM(1),1.0_dp)
 !   XM(1) = 1.0_dp
    JC=T%np

    do i=1,(T%N-T%np)/T%np
       !
 !      xx = xm(T%jl(JC+1))*xt(T%jV(JC+1))
       call damul_pancake(xm(T%jl(JC+1)),xt(T%jV(JC+1)),xx)

       call dacop_pancake(xx,xm(T%jl(JC+1)+1))
 !      xm(T%jl(JC+1)+1) = xx
       !
       do iv=1,T%np
          jc=jc+1

          call dacmu_pancake(xx,t%cc(jc),itemp)
 
          call daadd_pancake(xf(iv),itemp,xf(iv))
   !       xf(iv) = xf(iv) + t%cc(jc) * xx
       enddo
    enddo

    do i=1,size(xi)
       call dacop_pancake(xF(i),xI(i))
!       xI(i)=xF(i)
    enddo

    CALL KILL_pancake(XT)
    CALL KILL_pancake(XF)
    CALL KILL_pancake(XM)
    CALL dadal1_pancake(XX)
    CALL dadal1_pancake(itemp)

  END SUBROUTINE track_g_pancake


 

 

 
!!!!!!!!!!  pointed at layout !!!!!!!!!!!!!!
 

!!!!!!!!!!!!  switching routines !!!!!!!!!!!!!
  SUBROUTINE switch_layout_to_cavity( L,name,sext,a,r,freq,t )  ! switch to kind7
    implicit none
    TYPE (layout), target :: L
    TYPE (fibre), pointer :: p
    character(*) name
    real(dp),OPTIONAL:: a,r,freq,t
    INTEGER I
    logical(lp) sext


    p=>L%start
    do i=1,L%n

       call  switch_to_cavity( p,name,sext,a,r,freq,t)

       p=>p%next
    enddo



  end SUBROUTINE switch_layout_to_cavity



  SUBROUTINE switch_to_cavity( el,name,sext,a,r,freq,t )  ! switch to kind7
    implicit none
    TYPE (fibre), target :: el
    character(*) name
    integer i,nm,EXCEPTION
    real(dp),OPTIONAL:: a,r,freq,t
    real(dp), allocatable :: an(:),bn(:)
    type(keywords) key
    logical(lp) sext
    ! This routines switches to cavity
    nm=len_trim(name)
    select case(el%mag%kind)
    case(kind10,kind16,kind2,kind7,kind6,KIND20)
       if(el%mag%name(1:nm)==name(1:nm)) then
          if(sext.and.el%mag%p%nmul>2)then
             write(6,*) el%mag%name
             call add(el,3,0,0.0_dp)
             call add(el,-3,0,0.0_dp)
          else
             write(6,*) el%mag%name, " not changed "
          endif
       endif
       if(el%mag%p%b0/=0.0_dp.or.el%mag%name(1:nm)==name(1:nm)) return
       write(6,*) el%mag%name
       if(el%mag%p%nmul/=size(el%mag%an)) then
          write(6,*) "error in switch_to_cavity "
          stop 666
       endif
       allocate(an(el%mag%p%nmul),bn(el%mag%p%nmul))
       an=el%mag%an*el%mag%p%p0c
       bn=el%mag%bn*el%mag%p%p0c
       call zero_key(key)
       key%magnet="rfcavity"
       key%list%volt=0.0_dp
       key%list%lag=0.0_dp
       key%list%freq0=0.0_dp
       IF(PRESENT(FREQ)) THEN
          key%list%freq0=FREQ
       ENDIF
       key%list%n_bessel=0
       key%list%harmon=1.0_dp
       key%list%l=el%mag%L
       key%list%name=el%mag%name
       key%list%vorname=el%mag%vorname
       EXACT_MODEL=el%mag%p%exact
       key%nstep=el%mag%p%nst
       key%method=el%mag%p%method


       el%mag=-1
       el%magp=-1
       el%mag=0
       el%magp=0
       call create_fibre(el,key,EXCEPTION,my_true)
       do i=size(an),1,-1
          call add(el,-i,0,an(i))
          call add(el,i,0,bn(i))
       enddo
       el%mag%c4%a=1.0_dp
       el%magp%c4%a=1.0_dp
       IF(PRESENT(a)) THEN
          el%mag%c4%a=a
          el%magp%c4%a=a
       ENDIF
       el%mag%c4%r=0.0_dp
       el%magp%c4%r=0.0_dp
       IF(PRESENT(r)) THEN
          el%mag%c4%r=r
          el%magp%c4%r=r
       ENDIF
       el%mag%c4%PHASE0=0.0_dp
       el%mag%c4%PHASE0=0.0_dp
       el%mag%c4%always_on=my_true
       el%magp%c4%always_on=my_true
       el%mag%c4%xprime=my_false
       el%magp%c4%xprime=my_false
 
       IF(PRESENT(FREQ)) THEN
          el%mag%FREQ=FREQ
          el%magP%FREQ=FREQ
       ENDIF
       IF(PRESENT(T).and.PRESENT(FREQ)) THEN
          el%mag%C4%T=T/(el%mag%C4%freq/CLIGHT)
          el%magP%C4%T=el%mag%C4%T
       ENDIF
       deallocate(an,bn)
    CASE(KIND4)
       IF(el%mag%c4%always_on) THEN
          IF(PRESENT(FREQ)) THEN
             el%mag%FREQ=FREQ
             el%magP%FREQ=FREQ
          ENDIF
          IF(PRESENT(T)) THEN
             el%mag%C4%T=T/(el%mag%C4%freq/CLIGHT)
             el%magP%C4%T=el%mag%C4%T
          ENDIF
          IF(PRESENT(r)) THEN
             el%mag%c4%r=r
             el%magp%c4%r=r
          ENDIF
          IF(PRESENT(a)) THEN
             el%mag%c4%a=a
             el%magp%c4%a=a
          ENDIF

       ENDIF

    case default
       return
    end select

  END SUBROUTINE switch_to_cavity

  SUBROUTINE switch_to_kind7( el )  ! switch to kind7
    implicit none
    TYPE (fibre), target :: el
    ! This routines switches to kind7 (not exact) from kind2,10,16
    select case(el%mag%kind)
    case(kind10,kind16,kind2,KIND20)
       el%magp=-1
       el%mag%L=el%mag%p%ld
       el%mag%p%lc=el%mag%p%ld
       el%mag%p%exact=my_false
       el%magp=0
    end select

    select case(el%mag%kind)
    case(kind10)
       EL%MAG%TP10=-1
       deallocate(EL%MAG%TP10)
       el%mag%kind=KIND7
       CALL SETFAMILY(EL%MAG)
       CALL COPY(EL%MAG,EL%MAGP)
    case(kind16,KIND20)
       EL%MAG%k16=-1
       deallocate(EL%MAG%k16)
       el%mag%kind=KIND7
       CALL SETFAMILY(EL%MAG)
       CALL COPY(EL%MAG,EL%MAGP)
    case(KIND2)
       el%mag%kind=KIND7
       CALL SETFAMILY(EL%MAG)
       CALL COPY(EL%MAG,EL%MAGP)
    end select

  END SUBROUTINE switch_to_kind7

!!!!!! New flat file for data base universe  M_T !!!!!!!!!!!!!!!!! 
  
subroutine  print_new_flat(ring,filename,last,com)

implicit none
type(layout), target :: ring
type(fibre), pointer :: f
logical(lp), optional ::last
character(6), optional ::com
character(*) filename
integer i,mf
character(120) line
logical(lp) fin
character (6) comt
comt='REWIND'
fin=my_true
if(present(last)) fin=last
!goto 1
if(present(com)) comt=com
call kanalnummer(mf)
open(unit=mf,file=filename,position=comt) !,recl=200) !comt could be append for a complex universe 

   write(mf,'(a120)') ring%name                        ! Sagan depedent line
   write(mf,*) highest_fringe  , " highest fringe "    !  Sagan depedent line DO NOT CHANGE
   write(mf,*) lmax, " Maximum Length for Orbit "
   write(MF,*) ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING  , "ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING "
   write(mf,*) SECTOR_NMUL,SECTOR_NMUL_MAX , " SECTOR_NMUL,SECTOR_NMUL_MAX "
    
 write(mf,*) " $$$$$$$$$$$$$$$$$ START OF LAYOUT $$$$$$$$$$$$$$$$$"


  
f=>ring%start

call Print_initial_chart(f,mf)

do i=1,ring%n
  call el_el0(f%mag,my_true,mf)
  call fib_fib0(f,my_true,mf)
  CALL MC_MC0(f%MAG%P,my_true,mf)
  CALL print_ElementLIST(f%mag,MY_TRUE,mf)
  if(f%patch%patch/=0.or.f%patch%time/=0.or.f%patch%energy/=0) call patch_patch0(f%patch,my_true,mf)
  if(f%mag%mis) call CHART_CHART0(f%chart,my_true,mf)
 write(mf,*) " $$$$$$$$$$$$$$$$$ END OF FIBRE $$$$$$$$$$$$$$$$$"
 f=>f%next    
enddo

 write(mf,*) "&ELENAME"
if(fin) then
 write(mf,*)  "ELE0%NAME_VORNAME='alldone','alldone',"
else
 write(mf,*)  "ELE0%NAME_VORNAME='endhere','endhere',"
endif 
write(mf,*) "/"


close(mf)


end subroutine print_new_flat

subroutine  print_universe(un,filename)
!call print_universe(m_u,'junk2.txt')
!call print_universe_pointed(m_u,m_t,'junk3.txt')
implicit none
type(mad_universe),target :: un
type(layout),pointer :: r
character(*) filename
integer i



r=>un%start


if(associated(r,un%end)) then
 call  print_new_flat(r,filename,last=my_true)
else
          do while(.not.associated(r,un%end))
            if(associated(r,un%start)) then
              call  print_new_flat(r,filename,last=my_false)  
            else
              call  print_new_flat(r,filename,last=my_false,com="APPEND")  
            endif
            r=>r%next
          enddo
              call  print_new_flat(r,filename,last=my_true,com="APPEND")  
 
endif
call print_universe_siamese(un,filename)
call print_universe_girders(un,filename)
end subroutine  print_universe

subroutine  print_universe_girders(un,filename)
implicit none
type(mad_universe),target :: un
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
character(*) filename
integer i,k,i1,i2,j1,j2,MF



call TIE_MAD_UNIVERSE(un)

p=>un%start%start
p0=>p
p=>p%next


k=0
do while(.not.associated(p0,p))


if(associated(p%mag%girders)) then
 
 if(associated(p%mag%GIRDER_FRAME)) then
 k=k+1

endif
endif
 
 p=>p%n
enddo

call kanalnummer(mf)
open(unit=mf,file=filename,position='APPEND')
write(MF,*) k, " girders in the universe "




p=>un%start%start
p0=>p
p=>p%next


i=0
do while(.not.associated(p0,p))
i=i+1

if(associated(p%mag%girders)) then
 
 
 if(associated(p%mag%GIRDER_FRAME)) then
  j1=0
  j2=0
  ps=>p
  call locate_in_universe(ps,i1,i2)
     write(MF,*) p%mag%girder_frame%a    
     write(MF,*) p%mag%girder_frame%ent(1,1:3)
     write(MF,*) p%mag%girder_frame%ent(2,1:3)
     write(MF,*) p%mag%girder_frame%ent(3,1:3)
     write(MF,*) p%mag%girder_frame%b    
     write(MF,*) p%mag%girder_frame%exi(1,1:3)
     write(MF,*) p%mag%girder_frame%exi(2,1:3)
     write(MF,*) p%mag%girder_frame%exi(3,1:3)



     write(MF,*) i1,i2,ps%loc,PS%MAG%NAME
  do while(j1/=i1.or.j2/=i2)
      ps=>ps%mag%girders%parent_fibre    
      call locate_in_universe(ps,j1,j2)
      if (j1/=i1.or.j2/=i2) then
        write(MF,*) j1,j2 ,ps%loc
      else
        write(MF,*) 0,0,0
      endif
  enddo
 endif

endif
 
 p=>p%n
enddo


CLOSE(MF)

end subroutine  print_universe_girders

subroutine  read_universe_girders(un,mf,ns)
implicit none
type(mad_universe),target :: un
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
integer i,k,j,i1,i2,j1,MF,ns


call TIE_MAD_UNIVERSE(un)




do i=1,ns

call read_initial_chart(mf)

read(mf,*) i1,i2,j1

 p0=>un%start%start 
 do j=2,j1
  p0=>p0%n
 enddo
 ps=>p0

 do k=1,1000000 
  read(mf,*) i1,i2,j1
  if(i1==0) exit
  p=>un%start%start
  do j=2,j1
   p=>p%n
  enddo
  ps%mag%girders=>p%mag
  ps=>p
 enddo
  ps%mag%girders=>p0%mag

  allocate(p0%mag%girder_frame)
  call NULL_af(p0%mag%girder_frame)
  allocate(p0%mag%girder_frame%a(3))
  allocate(p0%mag%girder_frame%b(3))
  allocate(p0%mag%girder_frame%ent(3,3))
  allocate(p0%mag%girder_frame%exi(3,3))
  p0%mag%girder_frame%a=a_
  p0%mag%girder_frame%ent(1,1:3)=ent_(1,1:3)
  p0%mag%girder_frame%ent(2,1:3)=ent_(2,1:3)
  p0%mag%girder_frame%ent(3,1:3)=ent_(3,1:3)
  p0%mag%girder_frame%b=b_
  p0%mag%girder_frame%exi(1,1:3)=exi_(1,1:3)
  p0%mag%girder_frame%exi(2,1:3)=exi_(2,1:3)
  p0%mag%girder_frame%exi(3,1:3)=exi_(3,1:3)

enddo


end subroutine  read_universe_girders

subroutine  print_universe_siamese(un,filename)
implicit none
type(mad_universe),target :: un
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
character(*) filename
integer i,k,i1,i2,j1,j2,MF



call TIE_MAD_UNIVERSE(un)

p=>un%start%start
p0=>p
p=>p%next


k=0
do while(.not.associated(p0,p))


if(associated(p%mag%siamese)) then
 
 if(associated(p%mag%SIAMESE_FRAME)) then
 k=k+1

endif
endif
 
 p=>p%n
enddo

call kanalnummer(mf)
open(unit=mf,file=filename,position='APPEND') !,recl=200)
write(MF,*) k, " siamese in the universe "




p=>un%start%start
p0=>p
p=>p%next


i=0
do while(.not.associated(p0,p))
i=i+1

if(associated(p%mag%siamese)) then
 
 
 if(associated(p%mag%SIAMESE_FRAME)) then
  j1=0
  j2=0
  ps=>p
  call locate_in_universe(ps,i1,i2)
     write(MF,*) p%mag%SIAMESE_FRAME%ANGLE
     write(MF,*) p%mag%SIAMESE_FRAME%d
     write(MF,*) i1,i2,ps%loc,PS%MAG%NAME
  do while(j1/=i1.or.j2/=i2)
      ps=>ps%mag%siamese%parent_fibre    
      call locate_in_universe(ps,j1,j2)
      if (j1/=i1.or.j2/=i2) then
        write(MF,*) j1,j2 ,ps%loc
      else
        write(MF,*) 0,0,0
      endif
  enddo
 endif

endif
 
 p=>p%n
enddo


CLOSE(MF)

end subroutine  print_universe_siamese

subroutine  read_universe_siamese(un,mf,ns)
implicit none
type(mad_universe),target :: un
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
integer i,k,j,i1,i2,j1,MF,ns
real(dp) a(3),d(3)


call TIE_MAD_UNIVERSE(un)




do i=1,ns

read(mf,*) a
read(mf,*) d

read(mf,*) i1,i2,j1

 p0=>un%start%start 
 do j=2,j1
  p0=>p0%n
 enddo
 ps=>p0

 do k=1,1000000 
  read(mf,*) i1,i2,j1
  if(i1==0) exit
  p=>un%start%start
  do j=2,j1
   p=>p%n
  enddo
  ps%mag%siamese=>p%mag
  ps=>p
 enddo
  ps%mag%siamese=>p0%mag

  allocate(p0%mag%siamese_frame)
  call NULL_af(p0%mag%siamese_frame)
  allocate(p0%mag%siamese_frame%angle(3))
  allocate(p0%mag%siamese_frame%d(3))
  p0%mag%siamese_frame%angle=a
  p0%mag%siamese_frame%d=d
enddo


end subroutine  read_universe_siamese

subroutine  read_universe_database(un,filename,arpent)
! the universes should be empty
!call read_universe_database(m_u,'junk2.txt',arpent=my_false)
!call read_universe_pointed(M_u,M_t,'junk3.txt')
!call create_dna(M_u,m_t)
!arpent= false => the database should not be surveyed.
! DNA is automatically created in create_dna
implicit none
type(mad_universe),target :: un
logical(lp), optional :: arpent
character(*) filename
integer mf,ns
ELE0%NAME_VORNAME(1)=' '
ELE0%NAME_VORNAME(2)=' '

 call kanalnummer(mf,filename(1:len_trim(filename)))
          do while(ELE0%NAME_VORNAME(1)/="alldone")
           call append_empty_layout(un)  
           call set_up(un%end)
           call read_lattice(un%end,filename,mf,arpent)

          enddo
read(mf,*) ns   ! number of siamese
 call read_universe_siamese(un,mf,ns)
read(mf,*) ns   ! number of girders
call read_universe_girders(un,mf,ns)
close(mf)
end subroutine  read_universe_database

subroutine  read_lattice_append(un,filename,arpent)
! the universes should be empty
!call read_universe_database(m_u,'junk2.txt',arpent=my_false)
!call read_universe_pointed(M_u,M_t,'junk3.txt')
!call create_dna(M_u,m_t)
!arpent= false => the databaseshould not be surveyed.
! DNA is automatically created in create_dna
implicit none
type(mad_universe),target :: un
logical(lp), optional :: arpent
character(*) filename
integer mf
ELE0%NAME_VORNAME(1)=' '
ELE0%NAME_VORNAME(2)=' '
 call kanalnummer(mf,filename(1:len_trim(filename)),old=.true.)

           call append_empty_layout(un)  
           call set_up(un%end)
           call read_lattice(un%end,filename,mf,arpent)
close(mf)
end subroutine  read_lattice_append



subroutine  read_lattice(r,filename,mfile,arpent)
implicit none
type(layout),target :: r
character(*) filename
logical(lp), optional :: arpent
integer , optional :: mfile
logical(lp) doneit,surv,do_extrasurvey,change
character(120) line
type(fibre),pointer :: s22
type(element),pointer :: s2
type(elementp), pointer :: s2p
integer se2,se1
logical :: excess,switch
integer mf,n,met,nst  !,met0,nst0
character(120) fff
excess=.true.
switch=.true.

if(present(mfile)) then
 mf=mfile
else
 call kanalnummer(mf,filename(1:len_trim(filename)),old=.true.)
endif
surv=my_true





!-----------------------------------


   read(mf,'(a120)') r%name
   read(mf,*) highest_fringe  
   read(mf,*) lmax  
   read(MF,*) ALWAYS_EXACTMIS,ALWAYS_EXACT_PATCHING  
   read(mf,*) se2,se1
     call input_sector(se2,se1)
call make_states(my_false)
call set_mad(energy=2.0e0_dp)


   
 read(mf,'(a120)') line


 
 read(mf,'(a120)') line
call read_initial_chart(mf)
 read(mf,'(a120)') line



n=0
do while(.true.) 
    call zero_ele0
   read(mf,NML=ELEname,end=999) !!  basic stuff : an,bn,L, b_sol, etc...
!write(6,*) ELE0%name_vorname


   if(ELE0%NAME_VORNAME(1)== "endhere".or.ELE0%NAME_VORNAME(1)=="alldone") goto 99
 !write(6,NML=ELEname)
   call zero_fib0
   read(mf,NML=FIBRENAME,end=999) !
!     real(dp) GAMMA0I_GAMBET_MASS_AG(4)  !GAMMA0I,GAMBET,MASS ,AG  BETA0 is computed
 !    real(dp)  CHARGE
 !    integer  DIR !DIR,CHARGE
 !    integer patch 
   call zero_MAGL0
   read(mf,NML=MAGLNAME,end=999) ! reads magnet frame: kill_fringe, nst, method, etc... 
 !write(6,NML=MAGLNAME)
 call read_ElementLIST(ELE0%kind,MF)  ! read individual elments and aperture at the end

                   
 if(fib0%patch/=0) then
  call zero_patch0
  read(mf,NML=patchname,end=999)    ! patch read if present
 endif

if(ele0%recut_even_electric_MIS(4)) then
 call zero_CHART0
 read(mf,NML=CHARTname)  ! reading misalignment
endif


 read(mf,'(a120)') line

 call append_empty(r)  ! appends a fibre
 s22=>r%end
  call  nullify_for_madx(s22)


  

    s2=>s22%mag;
    s2p=>s22%magp;



    call fib_fib0(s22,my_false)
IF(ELE0%kind==KIND10) THEN
 IF(MAGL0%METHOD_NST_NMUL_permfringe_highest(3)>SECTOR_NMUL_MAX) MAGL0%METHOD_NST_NMUL_permfringe_highest(3)=SECTOR_NMUL_MAX
ENDIF
     S2 = MAGL0%METHOD_NST_NMUL_permfringe_highest(3)    
     


!write(6,*) associated(s2%p)
 !pause 78
    call MC_MC0(s2%p,my_false)

 

    call el_el0(s2,dir=my_false)

 
 change=.false.
metwig=s2%p%method
nstwig=s2%p%nst
 
if(check_excessive_cutting) then
 call against_the_method(s2%p%method,s2%p%nst,met,nst,ELE0%kind,change)
  if(change) then
    if(lielib_print(17)==1) then
fff="((1x,i4,1x,i4,1x,a16,1x,i4,1x,i4,1x,i4))"
      write(6,fff) s2%p%method,s2%p%nst, ELE0%name_vorname(1)(1:16),met,nst,ELE0%kind
    endif
 !   s2%p%nst=s2%p%nst/faclim

if(excess) write(6,*) "At least one magnet had its method changed due to excessive cutting ",ELE0%name_vorname(1)
   excess=.false.
  endif
 !if(switch_to_drift_kick.and.change) then
 ! ! Restoring the user input times faclim
 ! s2%p%method=met
 ! s2%p%nst=nst*faclim
 !endif
! if(lielib_print(17)==1) write(6,*) 1,s2%p%method,s2%p%nst
     if(switch_to_drift_kick.and.change) then
       if(switch) write(6,*) "all changed magnets are switched to drift-kick-drift"
     if(lielib_print(17)==1) write(6,*)s2%p%method,s2%p%nst, ELE0%name_vorname(1)
      
!     s2%p%nst=s2%p%nst/faclim
      switch=.false.
      if(lielib_print(17)==1) then
      fff="((a30,i4,1x,i4,1x))"
      write(6,fff) "changed to drift-kick-drift ",s2%p%method,s2%p%nst
      endif
     endif
! else
 elseif(switch_to_drift_kick) then
       if(switch) write(6,*) "all  magnets are switched to drift-kick-drift"
     nst=s2%p%nst
!     s2%p%nst=s2%p%nst/faclim
     change=.true.
      switch=.false.
       if(lielib_print(17)==1) then
      fff="((a30,1x,i4,1x,i4,a4,1x,i4))"
           write(6,fff)"changed to drift-kick-drift ",s2%p%method,nst," - > ",s2%p%nst
        endif
    ! endif
else
 
 call against_the_method(s2%p%method,s2%p%nst,met,nst,ELE0%kind,change)
  if(change) then
   if(excess.or.(lielib_print(17)==1)) then
     write(6,*) " Looks like excessive cutting might take place "
     write(6,*) ELE0%name_vorname(1)
     write(6,*) " met0,nst0,met,nst ", metwig,nstwig,met,nst
     excess=.false.
   endif
  endif
s2%p%method=metwig 
s2%p%nst=nstwig

endif

 if(change.and.lielib_print(17)==1) write(6,*) " $$$$$$$$$$$$$$$$$$$$$$$$$$ "
 

    if(s2%kind/=kindpa) then
       CALL SETFAMILY(S2) 
    else
  !    read(mf,'(a)') line
  !   line=line(1:len_trim(line))
       if(tpsa_pancake) then
             call read_pancake_new(s2,s2%vorname,t_e,S2%P%METHOD)
       else
             call read_pancake_new_integer(s2,s2%vorname,t_e,S2%P%METHOD)
       endif
       CALL SETFAMILY(S2,t=T_E)
   !    S2%P%METHOD=4
       s2%pa%angc=angcsp
       s2%pa%xc=xcsp
       s2%pa%dc=dcsp
       s2%pa%hc=hcsp
       s2%pa%vc=vcsp
       s2%pa%xprime=xprime_pancake
       deallocate(T_E)
!     read(mf,'(a)') line
 !    write(6,*) line(1:len_trim(line))
     endif    

 


    call print_ElementLIST(s2,dir=my_false)

    s2p=0   
 
 !pause 665

    call copy(s2,s2p)
    s22%mag%parent_fibre =>s22
    s22%magp%parent_fibre=>s22
       s22%CHART=0
       s22%PATCH=0
     if(fib0%patch/=0) then
       !s22%PATCH%patch=fib0%patch
      call patch_patch0(s22%patch,my_false)
    endif
   if(ele0%recut_even_electric_MIS(4)) call CHART_CHART0(s22%chart,my_false)





n=n+1
enddo


   100 CONTINUE


99 write(6,*) ELE0%NAME_VORNAME(1)
999 write(6,*) "Read ",n

 s22=>r%start
! if(s22%dir==1) then
  s22%chart%f%a=a_
  s22%chart%f%ent=ent_
! else
  s22%chart%f%b=b_
  s22%chart%f%exi=exi_
! endif
    r%closed=.true.

    doneit=.true.
    call ring_l(r,doneit)
 if(present(arpent)) surv=arpent

   if(surv) then
     call check_mis_presence(r,do_extrasurvey)
     if(do_extrasurvey) then
      do_mis=.false.
      call survey(r)
      do_mis=.true.
     endif
     call survey(r)
   endif

1000 continue

if(.not.present(mfile)) close(mf)

end subroutine read_lattice

subroutine against_the_method(m,n,met,nst,kind00,change)
implicit none
!type(element), target :: s2
integer, intent(inout) :: met,nst
integer  m,n,m0,n0,nt
integer kind00
logical change
change=.false.
m0=m  !s2%p%method
n0=n  !s2%p%nst
met=m  !s2%p%method
nst=n  !s2%p%nst

 if(kind00==kindpa.or.kind00==kind0) return
if(m<=2) then
 if(n>limit_int0_new(1).and.n<=limit_int0_new(2)) then
 n=n/3
 nt=n*faclim
 if(nt>0) n=nt
 m=4
 change=.true.
if(.not.check_excessive_cutting) then
 met=m
 nst=n
 m=m0  !s2%p%method
 n=n0  !s2%p%nst
endif
 return
 endif
endif

 
if(m<=4) then
 if(n>limit_int0_new(2).and.n<=limit_int0_new(3)) then
 n=n/7
 nt=n*faclim
 if(nt>0) n=nt
 m=6
 change=.true.
if(.not.check_excessive_cutting) then
 met=m
 nst=n
 m=m0  !s2%p%method
 n=n0  !s2%p%nst
endif
return
 endif
endif

if(m<=6) then
if(n>limit_int0_new(3)) then
 change=.true.
 if(kind0==kindwiggler) then
  n=n/7
  nt=n*faclim
  if(nt>0) n=nt
  m=6
else
 n=n/15
 nt=n*faclim
 if(nt>0) n=nt
 m=8
 endif
if(.not.check_excessive_cutting) then
 met=m
 nst=n
 m=m0  !s2%p%method
 n=n0  !s2%p%nst
endif
return
endif
endif

end subroutine against_the_method

subroutine check_mis_presence(r,do_extrasurvey)
implicit none
type(layout), target :: r
logical do_extrasurvey
type(fibre), pointer :: p
integer i

do_extrasurvey=.false.
p=>r%start
do i=1,r%n
if(p%mag%mis) then
 do_extrasurvey=.true.
 return
endif
p=>p%next
enddo


end subroutine check_mis_presence


  subroutine read_ElementLIST(kind,mf)
    implicit none

    integer mf,i,kind
    LOGICAL dir
    character*255 line

 
    select case(kind)
    CASE(KIND0,KIND1,kind2,kind6,kind7,kind8,kind9,KIND11:KIND15,kind17)
  case(kind3)
     read(mf,NML=thin30name)
    case(kind4)
     read(mf,NML=CAVname)
    case(kind5)
     read(mf,NML=sol50name)
    case(kind10)
      read(mf,NML=tp100name)
    case(kindabell)
      read(mf,NML=ab0name)
    case(kind16,kind20)

     read(mf,NML=k160name)

    case(kind18)

    case(kind19)


    case(kindhel)
      read(mf,NML=helname)
    case(kind21)
     read(mf,NML=tCAVname)
    case(KINDWIGGLER)
      read(mf,NML=wigname)
    case(KINDpa)

   case default
       write(MF,*) " not supported in print_specific_element",kind
 !      stop 101
    end select
    
    CALL  r_ap_aplist(mf) 


  END SUBROUTINE read_ElementLIST

subroutine  fib_fib0(f,dir,mf)
implicit none
type(fibre), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
! fib0%t(1)=f%BETA0
 fib0%GAMMA0I_GAMBET_MASS_AG(1)=f%GAMMA0I
 fib0%GAMMA0I_GAMBET_MASS_AG(2)=f%GAMBET
 fib0%GAMMA0I_GAMBET_MASS_AG(3)=f%MASS
 fib0%GAMMA0I_GAMBET_MASS_AG(4)=f%AG
 fib0%DIR=f%DIR
 fib0%CHARGE=f%CHARGE
 fib0%patch=f%patch%patch+7*f%patch%energy+49*f%patch%time
 
    if(present(mf)) then
     write(mf,NML=fibrename)
    endif   
else
    if(present(mf)) then
     read(mf,NML=fibrename)
    endif   
    ! f%BETA0=fib0%t(1)
 f%GAMMA0I=fib0%GAMMA0I_GAMBET_MASS_AG(1)
 f%GAMBET=fib0%GAMMA0I_GAMBET_MASS_AG(2)
 f%MASS=fib0%GAMMA0I_GAMBET_MASS_AG(3)
 f%AG=fib0%GAMMA0I_GAMBET_MASS_AG(4)
 f%BETA0=sqrt(1.0_dp-f%GAMMA0I**2)   
 f%DIR=fib0%DIR 
 f%CHARGE=fib0%CHARGE
 !f%patch%patch=fib0%patch     ! f%patch%patch is not yet allocated
endif
endif

end subroutine fib_fib0

subroutine  patch_patch0(f,dir,mf)
implicit none
type(PATCH), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
! fib0%t(1)=f%BETA0
 
 patch0%A_X1=f%A_X1
 patch0%A_X2=f%A_X2
 patch0%B_X1=f%B_X1
 patch0%B_X2=f%B_X2
 patch0%A_D=f%A_D
 patch0%B_D=f%B_D
 patch0%A_ANG=f%A_ANG
 patch0%B_ANG=f%B_ANG
 patch0%A_T=f%A_T
 patch0%B_T=f%B_T
 patch0%A_L=f%A_L
 patch0%B_L=f%B_L
 patch0%ENERGY=f%ENERGY
 patch0%TIME=f%TIME
 patch0%geometry=f%patch
 patch0%track=f%track

    if(present(mf)) then
     write(mf,NML=patchname)
    endif   
else
    if(present(mf)) then
     read(mf,NML=patchname)
    endif   

f%A_X1= patch0%A_X1
f%A_X2= patch0%A_X2
f%B_X1= patch0%B_X1
f%B_X2= patch0%B_X2

 f%A_D=patch0%A_D
 f%B_D=patch0%B_D
 f%A_ANG=patch0%A_ANG
 f%B_ANG=patch0%B_ANG
 f%track=patch0%track
 f%A_T=patch0%A_T
 f%B_T=patch0%B_T
 f%A_L=patch0%A_L
 f%B_L=patch0%B_L
 f%ENERGY=patch0%ENERGY
 f%TIME=patch0%TIME
 f%patch=patch0%geometry

endif
endif
end subroutine patch_patch0


subroutine  CHART_CHART0(f,dir,mf)
implicit none
type(chart), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then    

 CHART0%D_IN=f%D_IN
 CHART0%D_OUT=f%D_OUT
 CHART0%ANG_IN=f%ANG_IN
 CHART0%ANG_OUT=f%ANG_OUT
 
    if(present(mf)) then
     write(mf,NML=CHARTname)
    endif   
else
    if(present(mf)) then
     read(mf,NML=CHARTname)
    endif   

 f%D_IN=CHART0%D_IN
 f%D_OUT=CHART0%D_OUT
 f%ANG_IN=CHART0%ANG_IN
 f%ANG_OUT=CHART0%ANG_OUT
 

endif
endif
end subroutine CHART_CHART0



subroutine  MC_MC0(f,dir,mf)
implicit none
type(MAGNET_CHART), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
 MAGL0%LC_LD_B0_P0(1)=f%LC
 MAGL0%LC_LD_B0_P0(2)=f%LD
 MAGL0%LC_LD_B0_P0(3)=f%B0
 MAGL0%LC_LD_B0_P0(4)=f%P0C
 
 MAGL0%TILTD_EDGE(1)=f%TILTD
 MAGL0%TILTD_EDGE(2)=f%EDGE(1)
 MAGL0%TILTD_EDGE(3)=f%EDGE(2)

 MAGL0%KILL_SPIN(1)=f%KILL_ENT_SPIN
 MAGL0%KILL_SPIN(2)=f%KILL_EXI_SPIN

 MAGL0%KIN_KEX_BENDFRINGE_EXACT(1)=f%KILL_ENT_FRINGE
 MAGL0%KIN_KEX_BENDFRINGE_EXACT(2)=f%KILL_EXI_FRINGE
 MAGL0%KIN_KEX_BENDFRINGE_EXACT(3)=f%bend_fringe
 MAGL0%KIN_KEX_BENDFRINGE_EXACT(4)=f%EXACT

 MAGL0%METHOD_NST_NMUL_permfringe_highest(1)=f%METHOD
 MAGL0%METHOD_NST_NMUL_permfringe_highest(2)=f%NST
 MAGL0%METHOD_NST_NMUL_permfringe_highest(3)=f%NMUL
 MAGL0%METHOD_NST_NMUL_permfringe_highest(4)=f%permfringe
 MAGL0%METHOD_NST_NMUL_permfringe_highest(5)=f%highest_fringe

 if(present(mf)) then
     write(mf,NML=MAGLname)
    endif   
else
    if(present(mf)) then
     read(mf,NML=MAGLname)
    endif   
 f%LC=MAGL0%LC_LD_B0_P0(1)
 f%LD=MAGL0%LC_LD_B0_P0(2)
 f%B0=MAGL0%LC_LD_B0_P0(3)
 f%P0C=MAGL0%LC_LD_B0_P0(4)
 
 f%TILTD=MAGL0%TILTD_EDGE(1)
 f%EDGE(1)=MAGL0%TILTD_EDGE(2)
 f%EDGE(2)=MAGL0%TILTD_EDGE(3)

 f%KILL_ENT_SPIN=MAGL0%KILL_SPIN(1)
 f%KILL_EXI_SPIN=MAGL0%KILL_SPIN(2)

 f%KILL_ENT_FRINGE=MAGL0%KIN_KEX_BENDFRINGE_EXACT(1)
 f%KILL_EXI_FRINGE=MAGL0%KIN_KEX_BENDFRINGE_EXACT(2)
 f%bend_fringe=MAGL0%KIN_KEX_BENDFRINGE_EXACT(3)
 f%EXACT=MAGL0%KIN_KEX_BENDFRINGE_EXACT(4)

 f%METHOD=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(1)
 f%NST=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(2)
 f%NMUL=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(3)
 f%permFRINGE=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(4)
 f%highest_fringe=MAGL0%METHOD_NST_NMUL_permFRINGE_highest(5)

endif
endif
end subroutine MC_MC0

subroutine  el_el0(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf
character(nlp+3) nname
integer n,np,no,inf,i
!logical, optional :: ch
!logical change

!change=.false.
!if(present(ch)) change=ch

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
 ELE0%KIND=F%KIND
! ELE0%name_vorname(1)="'"//f%name//"' "
! ELE0%name_vorname(2)="'"//f%vorname//"' "

 ELE0%name_vorname(1)= f%name 
 ELE0%name_vorname(2)= f%vorname 
  call context(ELE0%name_vorname(1))
  call context(ELE0%name_vorname(2))
 if(.not.old_name_vorname) then
  call context(ELE0%name_vorname(1),dollar=my_true)
  call context(ELE0%name_vorname(2),dollar=my_true)
 endif
 ele0%an=0.0_dp
 ele0%bn=0.0_dp
if(f%p%nmul>0) then
  ele0%an(1:f%p%nmul)=f%an(1:f%p%nmul)
  ele0%bn(1:f%p%nmul)=f%bn(1:f%p%nmul)
endif
 ele0%VOLT_FREQ_PHAS=0.0_dp
 ele0%B_SOL=0.0_dp
 
   ele0%fint_hgap_h1_h2_va_vs(1)=f%fint(1)
   ele0%fint_hgap_h1_h2_va_vs(2)=f%fint(2)
   ele0%fint_hgap_h1_h2_va_vs(3)=f%hgap(1)
   ele0%fint_hgap_h1_h2_va_vs(4)=f%hgap(2)
   ele0%fint_hgap_h1_h2_va_vs(5)=f%h1
   ele0%fint_hgap_h1_h2_va_vs(6)=f%h2
   ele0%fint_hgap_h1_h2_va_vs(7)=f%va
   ele0%fint_hgap_h1_h2_va_vs(8)=f%vs


   ele0%L=f%L
   IF(ASSOCIATED(f%B_SOL)) ele0%B_SOL=f%B_SOL
 
 if(associated(f%volt)) ele0%VOLT_FREQ_PHAS(1)=f%VOLT
 if(associated(f%FREQ)) ele0%VOLT_FREQ_PHAS(2)=f%FREQ
 if(associated(f%PHAS)) ele0%VOLT_FREQ_PHAS(3)=f%PHAS
 if(associated(f%THIN)) ele0%THIN=f%THIN


ele0%slow_ac= f%slow_ac
ele0%recut_even_electric_MIS(1) = f%recut
ele0%recut_even_electric_MIS(2) = f%even
ele0%recut_even_electric_MIS(3) = f%electric
ele0%recut_even_electric_MIS(4) = f%MIS
 ele0%usebf_do1bf(1)=f%useb
 ele0%usebf_do1bf(2)=f%usef 
 ele0%skipptcbf(1)=f%skip_ptc_b 
 ele0%skipptcbf(2)=f%skip_ptc_f 
 ele0%usebf_do1bf(3)=f%do1mapb 
 ele0%usebf_do1bf(4)=f%do1mapf
 ele0%filef=trim(f%filef)
 ele0%fileb=TRIM(f%fileb)
 



if(associated(f%forward)) then

  if(present(mf)) then
   call kanalnummer(inf,f%filef)
    call print_tree_elements(f%forward,inf)
   close(inf)
 endif
endif
if(associated(f%backward)) then
 if(present(mf)) then
   call kanalnummer(inf,f%fileb)
    call print_tree_elements(f%backward,inf)
   close(inf)
 endif
endif


    if(present(mf)) then
     write(mf,NML=ELEname)
    endif  


else



    if(present(mf)) then
     read(mf,NML=ELEname)
    endif   

if(switch_to_drift_kick) then
 if(ELE0%KIND==kind7.or.ELE0%KIND==kind6 ) then
 IF(.not.f%p%exact) then
  ELE0%KIND = kind2
 endif
 ! if(lielib_print(17)==1)write(6,*) "element ",ELE0%name_vorname(1)," changed to drift-kick "
endif
endif

 F%KIND=ELE0%KIND  
  call context(ELE0%name_vorname(1))
  call context(ELE0%name_vorname(2))
nname=ELE0%name_vorname(1)
 f%name=nname(1:nlp)
nname=ELE0%name_vorname(2)
 f%vorname=nname(1:nlp)
if(f%p%nmul>0) then
 f%an=0.0_dp
 f%bn=0.0_dp
 f%an(1:f%p%nmul)=ele0%an(1:f%p%nmul)
 f%bn(1:f%p%nmul)=ele0%bn(1:f%p%nmul)
endif

f%fint(1)= ele0%fint_hgap_h1_h2_va_vs(1)
f%fint(2)=ele0%fint_hgap_h1_h2_va_vs(2)
f%hgap(1)= ele0%fint_hgap_h1_h2_va_vs(3)
f%hgap(2)= ele0%fint_hgap_h1_h2_va_vs(4)
f%h1  = ele0%fint_hgap_h1_h2_va_vs(5)
f%h2  = ele0%fint_hgap_h1_h2_va_vs(6)
f%va  = ele0%fint_hgap_h1_h2_va_vs(7)
f%vs  = ele0%fint_hgap_h1_h2_va_vs(8)


if(f%kind==kind4.or.f%kind==kind21) then ! cavities
 if(.not.associated(f%volt)) then
  ALLOCATE(f%VOLT,f%FREQ,f%PHAS,f%DELTA_E,f%THIN,f%lag)
 endif
  f%VOLT=ele0%VOLT_FREQ_PHAS(1)
  f%FREQ=ele0%VOLT_FREQ_PHAS(2)
  f%PHAS=ele0%VOLT_FREQ_PHAS(3)
  f%THIN=ele0%THIN
  f%delta_e=0.0_dp
endif

if(f%kind==KIND15) then   ! electrip separetor
 if(.not.associated(f%volt)) then
  ALLOCATE(f%VOLT,f%PHAS)
 endif
  f%VOLT=ele0%VOLT_FREQ_PHAS(1)
  f%PHAS=ele0%VOLT_FREQ_PHAS(3)
endif

    if(f%kind==kind22) then  ! helical dipole
       if(.not.associated(f%freq)) ALLOCATE(f%FREQ,f%PHAS)
       f%FREQ=ele0%VOLT_FREQ_PHAS(2)
       f%PHAS=ele0%VOLT_FREQ_PHAS(3)
    endif

 f%slow_ac = ele0%slow_ac
 f%recut = ele0%recut_even_electric_MIS(1)
 f%even = ele0%recut_even_electric_MIS(2)
 f%electric = ele0%recut_even_electric_MIS(3)
 f%MIS = ele0%recut_even_electric_MIS(4)
 solve_electric=f%electric
   F%L=ele0%L



!     logical(lp) usebf_do1bf(4)!
!	 integer skipptcbf(2)

 f%useb=ele0%usebf_do1bf(1)
 f%usef=ele0%usebf_do1bf(2)
 f%skip_ptc_b=ele0%skipptcbf(1)
 f%skip_ptc_f=ele0%skipptcbf(2)
 f%do1mapb=ele0%usebf_do1bf(3) 
 f%do1mapf=ele0%usebf_do1bf(4)
 f%fileb= ele0%fileb
 f%filef= ele0%filef
if(ele0%filef/=' '.and.readingmaps) then
 if(.not.associated(f%forward)) then 
  allocate(f%forward(3))
 endif
 call kanalnummer(inf,ele0%filef)

  do i=1,3
    read(inf,*) n,np,no
    CALL ALLOC_TREE(f%forward(i),N,NP)
    f%forward(i)%N=n
    f%forward(i)%NP=np
    f%forward(i)%no=no
    call read_tree_element(f%forward(i),inf)
  enddo
close(inf)
endif

if(ele0%fileb/=' '.and.readingmaps) then
 if(.not.associated(f%backward)) then 
  allocate(f%backward(3))
 endif
 call kanalnummer(inf,ele0%fileb)

  do i=1,3
    read(inf,*) n,np,no
    CALL ALLOC_TREE(f%backward(i),N,NP)
    f%backward(i)%N=n
    f%backward(i)%NP=np
    f%backward(i)%no=no
    call read_tree_element(f%backward(i),inf)
  enddo
close(inf)
 endif
   
    if(f%kind==kind3.or.f%kind==kind5) then   
        IF(.not.ASSOCIATED(f%B_SOL)) ALLOCATE(f%B_SOL);
       F%B_SOL=ele0%B_SOL
    endif

   
endif
endif
end subroutine el_el0



  subroutine print_ElementLIST(el,dir,mf)
    implicit none
    type(element), pointer :: el
    integer i
     logical(lp),optional ::  dir
     integer,optional :: mf
    character*255 line
!     logical(lp),optional ::  ch
!    logical(lp) change
!    change=.false.
!    if(present(ch)) change=ch
    select case(el%kind)
    CASE(KIND0,KIND1,kind2,kind6,kind7,kind8,kind9,KIND11:KIND15,kind17)
  case(kind3)
     call thin3_thin30(el,dir,mf)
    case(kind4)
        call cav4_cav40(EL,dir,mf)
    case(kind5)
        call sol5_sol50(EL,dir,mf)
    case(kind10)
        call tp10_tp100(EL,dir,mf)
 
    case(kindabell)
        call ab_ab0(EL,dir,mf)

    case(kind16,kind20)
        call k16_k160(EL,dir,mf)

    case(kind18)
!       WRITE(MF,*) " RCOLLIMATOR HAS AN INTRINSIC APERTURE "
  !     CALL  ap_aplist(el,dir,mf) 
    case(kind19)
!       WRITE(MF,*) " ECOLLIMATOR HAS AN INTRINSIC APERTURE "
!       CALL print_aperture(EL%ECOL19%A,mf)
    case(kindhel)
        call hel_hel0(EL,dir,mf)
    case(kind21)
        call tcav4_tcav40(EL,dir,mf)
!       WRITE(MF,*) el%cav21%PSI,el%cav21%DPHAS,el%cav21%DVDS
    case(KINDWIGGLER)
        call wig_wig0(el,dir,mf)
 !      call print_wig(el%wi,mf)
    case(KINDpa)
       if(tpsa_pancake) then
        if(dir) call print_pancake_field(el%pa,el%vorname)
       else
        if(dir) call print_pancake_field_integer(el%pa,el%vorname)
       endif
    case default
       write(6,*) " not supported in print_specific_element",el%kind
 !      stop 101
    end select
 
    CALL  ap_aplist(el,dir,mf) 

 

  END SUBROUTINE print_ElementLIST

   
subroutine  cav4_cav40(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

 cav0%f=0.0_dp
 cav0%PH=0.0_dp   
 cav0%N_BESSEL=F%c4%N_BESSEL
 cav0%NF=F%c4%NF
 cav0%CAVITY_TOTALPATH=F%c4%CAVITY_TOTALPATH
 cav0%phase0=F%c4%phase0
 cav0%t=F%c4%t
 cav0%always_on=F%c4%always_on
 cav0%xprime=F%c4%xprime
 cav0%f(1:F%c4%NF)=F%c4%f
 cav0%PH(1:F%c4%NF)=F%c4%PH
 cav0%A=F%c4%A
 cav0%R=F%c4%R
    if(present(mf)) then
     write(mf,NML=CAVname)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=CAVname)
    endif   
N_CAV4_F=cav0%NF
CALL SETFAMILY(f)
 F%c4%N_BESSEL=cav0%N_BESSEL
 F%c4%NF =cav0%NF


 F%c4%CAVITY_TOTALPATH=cav0%CAVITY_TOTALPATH
 F%c4%phase0=cav0%phase0
 F%c4%t=cav0%t
 F%c4%always_on=cav0%always_on
 F%c4%xprime=cav0%xprime
 F%c4%f=cav0%f(1:F%c4%NF)

 
 F%c4%PH=cav0%PH(1:F%c4%NF)
 F%c4%A=cav0%A
 F%c4%R=cav0%R  
endif
endif
end subroutine cav4_cav40

subroutine  hel_hel0(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG


hel0%fake_shift=F%he22%fake_shift
hel0%N_BESSEL=F%he22%N_BESSEL

    if(present(mf)) then
     write(mf,NML=helname)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=helname)
    endif   
 
CALL SETFAMILY(f)
 F%he22%N_BESSEL=hel0%N_BESSEL
 F%he22%fake_shift=hel0%fake_shift

endif
endif
end subroutine hel_hel0

subroutine  wig_wig0(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf
integer n,ne,nmax,met,nst
logical change
n=0
ne=0
if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG


 wig0%internal=F%wi%internal
 wig0%offset=F%wi%w%offset
 wig0%ex=f%wi%w%ex
 wig0%ey=f%wi%w%ey
wig0%n=0
 if(associated(f%wi%w%a)) wig0%n=size(f%wi%w%a)
 n=wig0%n
wig0%a=0.0_dp
wig0%f=0.0_dp
wig0%form=0.0_dp
wig0%k=0.0_dp
if(n>0) then
 wig0%a(1:n)=f%wi%w%a(1:n)
 wig0%f(1:n)=f%wi%w%f(1:n)
 wig0%form(1:n)=f%wi%w%form(1:n)
 wig0%k(1:3,1:n)=f%wi%w%k(1:3,1:n)
else
 wig0%a(1:n)=0
 wig0%f(1:n)=0
 wig0%form(1:n)=0
 wig0%k(1:3,1:n)=0
endif

wig0%ne=0
 if(associated(f%wi%w%ae)) wig0%ne=size(f%wi%w%ae)
 ne=wig0%ne
if(ne>0) then
 wig0%ae(1:ne)=f%wi%w%ae(1:ne)
 wig0%fe(1:ne)=f%wi%w%fe(1:ne)
 wig0%forme(1:ne)=f%wi%w%forme(1:ne)
 wig0%ke(1:3,1:ne)=f%wi%w%ke(1:3,1:ne)
else
 wig0%ae(1:ne)=0
 wig0%fe(1:ne)=0
 wig0%forme(1:ne)=0
 wig0%ke(1:3,1:ne)=0
endif

    if(present(mf)) then
     write(mf,NML=wigname)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=wigname)
    endif   
    if(.not.associated(f%wi%internal)) allocate(f%wi%internal(6))
  F%wi%internal=wig0%internal 
  N=wig0%N
  ne=wig0%Ne
  nmax=max(ne,n)
    call pointers_w(f%wi%w,N,ne)

 F%wi%w%offset=wig0%offset
 F%wi%w%ex=wig0%ex
 F%wi%w%ey=wig0%ey

if(n>0) then
 F%wi%w%a(1:N)=wig0%a(1:N)
 F%wi%w%f(1:N)=wig0%f(1:N)
 F%wi%w%form(1:N)=wig0%form(1:N)
 F%wi%w%k(1:3,1:N)=wig0%k(1:3,1:N)
endif

if(ne>0) then
 F%wi%w%ae(1:Ne)=wig0%ae(1:Ne)
 F%wi%w%fe(1:Ne)=wig0%fe(1:Ne)
 F%wi%w%forme(1:Ne)=wig0%forme(1:Ne)
 F%wi%w%ke(1:3,1:Ne)=wig0%ke(1:3,1:Ne)
endif
 
f%p%method=metwig
f%p%nst=nstwig
limit_int0_new=limit_int0_new*nmax
 call against_the_method(f%p%method,f%p%nst,met,nst,f%kind,change)
limit_int0_new=limit_int0_new/nmax
 
endif
endif


end subroutine wig_wig0

subroutine  tcav4_tcav40(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG


 tcav0%PSI_DPHAS_DVDS_phase0(1)=F%cav21%psi
 tcav0%PSI_DPHAS_DVDS_phase0(2)=F%cav21%dphas
 tcav0%PSI_DPHAS_DVDS_phase0(3)=F%cav21%dvds
 tcav0%PSI_DPHAS_DVDS_phase0(4)=F%cav21%phase0
 tcav0%always_on=F%cav21%always_on
 tcav0%implicit=F%cav21%implicit

    if(present(mf)) then
     write(mf,NML=tCAVname)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=tCAVname)
    endif   

 F%cav21%psi=tcav0%PSI_DPHAS_DVDS_phase0(1)
 F%cav21%dphas=tcav0%PSI_DPHAS_DVDS_phase0(2)
 F%cav21%dvds=tcav0%PSI_DPHAS_DVDS_phase0(3)
 F%cav21%phase0=tcav0%PSI_DPHAS_DVDS_phase0(4)
 F%cav21%always_on=tcav0%always_on
 F%cav21%implicit=tcav0%implicit

endif

endif
end subroutine tcav4_tcav40

subroutine  sol5_sol50(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

    


 sol50%dx_dy_pitchx_pitchy(1)=F%s5%dx
 sol50%dx_dy_pitchx_pitchy(2)=F%s5%dy
 sol50%dx_dy_pitchx_pitchy(3)=F%s5%pitch_x
 sol50%dx_dy_pitchx_pitchy(4)=F%s5%pitch_y

    if(present(mf)) then
     write(mf,NML=sol50name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=sol50name)
    endif   
       IF(.NOT.ASSOCIATED(f%B_SOL)) then
        ALLOCATE(f%B_SOL)
          f%B_SOL=0.0_dp
       endif
       CALL SETFAMILY(f) 

 F%s5%dx= sol50%dx_dy_pitchx_pitchy(1)
 F%s5%dy= sol50%dx_dy_pitchx_pitchy(2)
 F%s5%pitch_x= sol50%dx_dy_pitchx_pitchy(3)
 F%s5%pitch_y= sol50%dx_dy_pitchx_pitchy(4)
endif
endif
end subroutine sol5_sol50


subroutine  thin3_thin30(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

    

 thin30%thin_h_foc=F%k3%thin_h_foc
 thin30%thin_v_foc=F%k3%thin_v_foc
 thin30%thin_h_angle=F%k3%thin_h_angle
 thin30%thin_v_angle=F%k3%thin_v_angle
 thin30%hf=F%k3%hf
 thin30%vf=F%k3%vf
 thin30%patch=F%k3%patch
 thin30%ls=F%k3%ls 
  thin30%dx_dy_pitchx_pitchy(1)=F%k3%dx
  thin30%dx_dy_pitchx_pitchy(2)=F%k3%dy
  thin30%dx_dy_pitchx_pitchy(3)=F%k3%pitch_x
  thin30%dx_dy_pitchx_pitchy(4)=F%k3%pitch_y
    if(present(mf)) then
     write(mf,NML=thin30name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=thin30name)
    endif   
       IF(.NOT.ASSOCIATED(f%B_SOL)) then
        ALLOCATE(f%B_SOL)
          f%B_SOL=0.0_dp
       endif
       CALL SETFAMILY(f) 
 f%k3%thin_h_foc=thin30%thin_h_foc
 f%k3%thin_v_foc=thin30%thin_v_foc
 f%k3%thin_h_angle=thin30%thin_h_angle
 f%k3%thin_v_angle=thin30%thin_v_angle
 f%k3%hf=thin30%hf
 f%k3%vf=thin30%vf
 f%k3%patch=thin30%patch 
 f%k3%ls=thin30%ls
 F%k3%dx= thin30%dx_dy_pitchx_pitchy(1)
 F%k3%dy= thin30%dx_dy_pitchx_pitchy(2)
 F%k3%pitch_x= thin30%dx_dy_pitchx_pitchy(3)
 F%k3%pitch_y= thin30%dx_dy_pitchx_pitchy(4)
endif
endif
end subroutine thin3_thin30

subroutine  tp10_tp100(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf
!     logical(lp),optional ::  ch
!    logical(lp) change
!    change=.false.
!    if(present(ch)) change=ch 

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

 tp100%ae=0.0_dp
 tp100%be=0.0_dp
 
 tp100%DRIFTKICK=F%tp10%DRIFTKICK
 if(f%electric) then
  tp100%ae(1:size(F%tp10%ae))=F%tp10%ae
  tp100%be(1:size(F%tp10%be))=F%tp10%be
 endif

     if(present(mf)) then
     write(mf,NML=tp100name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=tp100name)
    endif   
   CALL SETFAMILY(f)
 if(f%electric) then
  F%tp10%ae=tp100%ae(1:sector_nmul_max) 
  F%tp10%be=tp100%be(1:sector_nmul_max)
  call GETAEBE(f%TP10)
 endif

 F%tp10%DRIFTKICK=tp100%DRIFTKICK
  if(switch_to_drift_kick) then
   F%tp10%DRIFTKICK=.true.
!    if(lielib_print(17)==1)write(6,*) "TP10 changed to drift-kick ",f%name
   endif
endif
endif
end subroutine tp10_tp100


subroutine  ab_ab0(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf
 

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG

 ab0%N_m=0 
 ab0%dz_t_te=0 
 ab0%b=0 
 ab0%e=0 
 ab0%scale_angc_xc_vc_dc_hc=0 
 
 ab0%n_m(1)= F%ab%n
 ab0%n_m(2)= F%ab%m
 ab0%b(1:ab0%n_m(2) ,1:ab0%n_m(2))= F%ab%b
 ab0%dz_t_te(2*ab0%n_m(2)+3:3*ab0%n_m(2)+3)=F%ab%te
ab0%dz_t_te(ab0%n_m(2)+2:2*ab0%n_m(2)+2)= F%ab%t
ab0%dz_t_te(1:ab0%n_m(2)+1)= F%ab%dz
ab0%scale_angc_xc_vc_dc_hc(1)=F%ab%SCALE 
 ab0%scale_angc_xc_vc_dc_hc(2)=F%ab%angc
 ab0%scale_angc_xc_vc_dc_hc(3)=F%ab%xc
 ab0%scale_angc_xc_vc_dc_hc(4)=F%ab%vc
 ab0%scale_angc_xc_vc_dc_hc(5)=F%ab%dc
 ab0%scale_angc_xc_vc_dc_hc(6)=F%ab%hc
 

! if(f%electric) then
!  tp100%ae(1:size(F%tp10%ae))=F%tp10%ae
!  tp100%be(1:size(F%tp10%be))=F%tp10%be
! endif

     if(present(mf)) then
     write(mf,NML=ab0name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=ab0name)
    endif  
   n_abell=ab0%n_m(1) 
   m_abell=ab0%n_m(2)
   CALL SETFAMILY(f)
 !if(f%electric) then
 ! F%tp10%ae=tp100%ae(1:sector_nmul_max) 
 ! F%tp10%be=tp100%be(1:sector_nmul_max)
 ! call GETAEBE(f%TP10)
 !endif
 F%ab%n=ab0%n_m(1) 
 F%ab%m=ab0%n_m(2) 
 F%ab%b=ab0%b(1:ab0%n_m(2) ,1:ab0%n_m(2))
 F%ab%te=ab0%dz_t_te(2*ab0%n_m(2)+3:3*ab0%n_m(2)+3)
 F%ab%t=ab0%dz_t_te(ab0%n_m(2)+2:2*ab0%n_m(2)+2)
 F%ab%dz=ab0%dz_t_te(1:ab0%n_m(2)+1)
F%ab%SCALE= ab0%scale_angc_xc_vc_dc_hc(1)
F%ab%angc= ab0%scale_angc_xc_vc_dc_hc(2)
F%ab%xc= ab0%scale_angc_xc_vc_dc_hc(3)
F%ab%vc= ab0%scale_angc_xc_vc_dc_hc(4)
F%ab%dc= ab0%scale_angc_xc_vc_dc_hc(5)
F%ab%hc= ab0%scale_angc_xc_vc_dc_hc(6)
 
endif
endif
end subroutine ab_ab0

subroutine  k16_k160(f,dir,mf)
implicit none
type(element), target :: f
logical(lp),optional ::  dir
integer,optional :: mf
!     logical(lp),optional ::  ch
!    logical(lp) change
!    change=.false.
!    if(present(ch)) change=ch 

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
 
 k160%DRIFTKICK=F%k16%DRIFTKICK
 k160%LIKEMAD=.true.  !F%k16%LIKEMAD
     if(present(mf)) then
     write(mf,NML=k160name)
    endif   
 
 else
    if(present(mf)) then
     read(mf,NML=k160name)
    endif   
!!! bug by etienne since this line was commented off
 F%k16%DRIFTKICK=k160%DRIFTKICK

 if(switch_to_drift_kick) then
   F%k16%DRIFTKICK=.true.
    if(lielib_print(17)==1) write(6,*) "K16 changed to drift-kick ",f%name
  endif
! F%k16%LIKEMAD=k160%LIKEMAD
endif
endif
end subroutine k16_k160

subroutine  ap_aplist(f,dir,mf)
implicit none
type(element), target :: f
type(MADX_APERTURE), pointer :: a
logical(lp),optional ::  dir
integer,optional :: mf
logical(LP) :: here
CHARACTER(120) LINE

here=associated(f%p%APERTURE)

if(present(dir)) then
if(dir) then   !BETA0,GAMMA0I,GAMBET,MASS ,AG
 
if(here) then
 a=>f%p%aperture

    aplist%kind=a%kind
    aplist%r=a%r
    aplist%x=a%x
    aplist%y=a%y
    aplist%dx=a%dx
    aplist%dy=a%dy
     if(present(mf)) then
     Write(mf,*) " APERTURE "  
     write(mf,NML=aperturename)
    endif   
else
 Write(mf,*) " NO APERTURE "
endif
 
 else   ! dir=false
 !  if(present(mf)) then     
 !     READ(MF,'(a120)') LINE 
 !     CALL CONTEXT(LINE)
 !  ENDIF
    IF(aplist%on) THEN
       IF(.NOT.HERE) THEN
          CALL alloc(f%p%aperture)
           a=>f%p%aperture
       ENDIF
    
    
    if(present(mf)) then     
     read(mf,NML=aperturename)
    endif   
      a%kind=aplist%kind
      a%r=aplist%r
      a%x=aplist%x
      a%y=aplist%y
      a%dx=aplist%dx
      a%dy=aplist%dy      
    ENDIF
endif
endif
end subroutine ap_aplist


subroutine  r_ap_aplist(mf)
implicit none
integer  mf
CHARACTER(120) LINE
   
      READ(MF,'(a120)') LINE 
      CALL CONTEXT(LINE)
    IF(LINE(1:2)/='NO') THEN
        read(mf,NML=aperturename)
       aplist%on=.true.
   !    write(6,nml=aperturename)
    else
     aplist%on=.false.
    endif
    
end subroutine r_ap_aplist

!!!!!!!!!!!!!!!!!!!!!    pointed lattices !!!!!!!!!!!!!!!!!!!!!

subroutine  print_universe_pointed(ud,un,filename,com)
! call print_universe_pointed(m_u,m_t,'junk3.txt')
implicit none
type(mad_universe),target :: ud,un
type(layout),pointer :: r
character(6), optional ::com
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
character(*) filename
integer i,j,i0,j0,i1,j1,jb,MF
character (6) comt
logical(lp) before,just
logical, allocatable :: dna(:)
character(nlp) name

allocate(dna(ud%n))

comt='REWIND'
if(present(com)) comt=com

call kanalnummer(mf)
open(unit=mf,file=filename,position=comt) !,recl=4000)

call TIE_MAD_UNIVERSE(ud)


r=>un%start

write(mf,*) un%n, ud%n," trackable and DNA Layouts"


do i=1,un%n
dna=.false.
write(mf,'(a120)') r%name
p0=>r%start
p=>p0
call locate_in_universe(p,i0,j0)
jb=j0
j1=j0
!   write(mf,*) i,r%n," New "
   write(mf,*)  r%n," Number elements in pointed layout "

dna(i0)=.true.

  call Print_initial_chart(p,mf)
 name=p%mag%name
 if(.not.old_name_vorname) then
  call context(name,dollar=my_true)
 endif
  write(mf,'(1x,i4,1x,i8,1x,i2,1x,a24,1x,i8)') i0,p%dir*j0,p%patch%patch,name, 1
  call fib_fib0(p,my_true,mf)
  before=my_false
  just=my_false
 if(p%patch%patch/=0) call patch_patch0(p%patch,my_true,mf)

 do j=2,r%n
  p=>p%next
    jb=j1
   call locate_in_universe(p,i1,j1)
name=p%mag%name
 if(.not.old_name_vorname) then
  call context(name,dollar=my_true)
 endif
 write(mf,'(1x,i4,1x,i8,1x,i2,1x,a24,1x,i8)') i1,p%dir*j1,p%patch%patch,name, j
    dna(i1)=my_true      
       just=my_false
       if(before) then 
        call fib_fib0(p,my_true,mf)
        just=my_true
        before=my_false
       endif
    if(p%patch%patch/=0) then
       if(.not.just) call fib_fib0(p,my_true,mf)  
       call patch_patch0(p%patch,my_true,mf)   
       before=my_true   
    endif
 enddo
 write(mf,*) " DNA SUMMARY "
do i0=1,ud%n
 write(mf,*) i0,dna(i0)
enddo
write(mf,*) " !!!!!!! End of Pointed Layout !!!!!!!"

 r=>r%next

enddo

close(mf)

deallocate(dna)

end subroutine  print_universe_pointed


subroutine  read_universe_pointed(ud,un,filename)
implicit none
type(mad_universe),target :: ud,un
type(layout),pointer :: r,rd
type(fibre),pointer :: p,p0,ps
type(element),pointer :: m,m0
character(*) filename
integer i,j,i0,MF,n,n_u,k(3),mypause,ipause,n_d,size_dna
integer pos
character(120) line,name1
logical(lp) doneit,first
character(nlp) name
logical, allocatable :: dna(:)
 
call kanalnummer(mf,filename,old=.true.)
!open(unit=mf,file=filename)


call TIE_MAD_UNIVERSE(ud)


read(mf,*)n_u,n_d

allocate(dna(n_d))

do i=1,n_u
   read(mf,'(a120)') name1 
!   read(mf,*) i0,n 
   read(mf,*)  n 

  read(mf,'(a120)') line
 call read_initial_chart(mf)
  read(mf,'(a120)') line

call append_empty_layout(un) 
!call set_up(un%end)  !
 
       R => un%end
r%name=name1

 do j=1,n 
  !     read(mf,'(1x,i4,1x,i8,1x,i2,1x,a24)' ) k    ,name
       read(mf,*) k    ,name
!        write(6,*)  k  
!        write(6,*)   name
!pause 234
if(j==1.or.k(3)>0) then
 read(mf,NML=FIBRENAME)
else
 if(r%end%patch%patch>0) read(mf,NML=FIBRENAME)    ! previous had a patch
endif

       call MOVE_TO_LAYOUT_I( ud,rd,k(1) )
       pos=iabs(k(2))   
       call move_to_p_safe( rd,p,POS)
       call append_point(r, p)  
       if( p%mag%name/=name) then
         write(6,*) j," serious error in read_universe_pointed "
         write(6,*)  k ,name
         write(6,*) i,p%mag%name ,pos
           read(mf,'(a120)') line
           write(6,'(a120)') line
           read(mf,'(a120)') line
           write(6,'(a120)') line
           read(mf,'(a120)') line
           write(6,'(a120)') line
          ipause=mypause(666)
         stop 666   
       endif

        if(k(3)/=0) then
        read(mf,NML=patchname)
        call patch_patch0(r%end%patch,my_false)
        endif
        r%end%patch%patch=k(3)

       call fib_fib0(r%end,my_false) 
       if(k(2)>0) then   ! because fib0 could have the wrong dir
        r%end%dir=1
       else
        r%end%dir=-1
       endif
 enddo

    write(6,*) r%index,ud%n,un%n
    r%closed=my_true
    doneit=my_true
    call ring_l(r,doneit)

  p=>r%start
  p%chart%f%a=a_
  p%chart%f%ent=ent_
! else
  p%chart%f%b=b_
  p%chart%f%exi=exi_
call survey(r,pi=1,fi=r%n,a=r%start%chart%f%a,ent=r%start%chart%f%ent)
!call survey(INJECTION_and_RING,pi=1,fi=INJECTION_and_RING%n,ent=p_start_ring%t1%ent )
read(mf,'(a120)') line

size_dna=0
do j=1,n_d
 read(mf,*) i0,dna(j)
 if(dna(j)) size_dna=size_dna+1
enddo

 read(mf,'(a10)') line(1:10)

! Setting up the dna
 allocate(r%dna(size_dna))

  i0=0
 do j=1,n_d
  if(dna(j))  then
      i0=i0+1
      call MOVE_TO_LAYOUT_I( ud,rd,j )    
     r%dna(i0)%L=>rd
  endif
 enddo 
 if(i0/=size_dna) then
  write(6,*) i0,size_dna, "error in read_universe_pointed " 
  stop 567
 endif

enddo



close(mf)

deallocate(dna)

end subroutine  read_universe_pointed


 subroutine Print_initial_chart(f,mf)
 implicit none
 type(fibre), target :: f
 integer mf

 write(mf,*) " $$$$$$$$$$$$$$$$$ INITIAL CHART $$$$$$$$$$$$$$$$$"
!  IF(F%DIR==1) THEN
   write(mf,*) f%chart%f%A
   write(mf,*) f%chart%f%ENT(1,1:3)
   write(mf,*) f%chart%f%ENT(2,1:3)
   write(mf,*) f%chart%f%ENT(3,1:3)
!  ELSE
   write(mf,*) f%chart%f%B
   write(mf,*) f%chart%f%EXI(1,1:3)
   write(mf,*) f%chart%f%EXI(2,1:3)
   write(mf,*) f%chart%f%EXI(3,1:3)
!  ENDIF
 write(mf,*) " $$$$$$$$$$$$$$$$$ END OF INITIAL CHART $$$$$$$$$$$$$$$$$"

 end subroutine Print_initial_chart

 subroutine read_initial_chart(mf)
 implicit none
 integer mf
 
   read(mf,*) A_
   read(mf,*) ENT_(1,1:3)
   read(mf,*) ENT_(2,1:3)
   read(mf,*) ENT_(3,1:3)

   read(mf,*) B_
   read(mf,*) EXI_(1,1:3)
   read(mf,*) EXI_(2,1:3)
   read(mf,*) EXI_(3,1:3)

 end subroutine read_initial_chart


subroutine create_dna(ud,ut)
implicit none
type(mad_universe), target :: ut,ud
type(layout), pointer :: r,rd
type(fibre), pointer :: ps
integer i,j,j1,j2,k,kn
logical(lp), allocatable :: here(:)



call TIE_MAD_UNIVERSE(ut)

allocate(here(ud%n))

r=>ut%start
do i=1,ut%n

here=my_false

ps => r%start

do j=1,r%n
 call locate_in_universe(ps,j1,j2)
 here(j1)=my_true
ps=>ps%next
enddo

kn=0
do k=1,size(here)
 if(here(k)) kn=kn+1
enddo

if(associated(r%DNA)) then 
deallocate(r%DNA)
write(6,*) "deallocated DNA"
endif
allocate(r%DNA(kn))

rd=>ud%start
kn=0
do k=1,size(here)
 if(here(k)) then
    kn=kn+1
    r%dna(kn)%L=>rd
    r%dna(kn)%counter=k
  endif
 rd=>rd%next
enddo


 r=>r%next

enddo

deallocate(here)

end subroutine create_dna

subroutine zero_ele0
implicit none

    ele0%name_vorname=' ' 
	ele0%L=0;ele0%B_SOL=0;
	ele0%an=0;ele0%bn=0;
    ele0%VOLT_FREQ_PHAS=0
    ele0%THIN=.false. 
    ele0%fint_hgap_h1_h2_va_vs=0
	ele0%recut_even_electric_MIS=.false.
    ele0%slow_ac=0
    ele0%usebf_do1bf=.false.
    ele0%skipptcbf=0
    ele0%filef=' '
    ele0%fileb=' '


!     logical(lp) usebf_do1bf(4)!
!	 integer skipptcbf(2)

end subroutine zero_ele0

subroutine zero_fib0
implicit none

    fib0%GAMMA0I_GAMBET_MASS_AG=0 
    fib0%CHARGE=0
    fib0%DIR=0
    fib0%patch=0

end subroutine zero_fib0

subroutine zero_CHART0
implicit none
    CHART0%D_IN=0; CHART0%D_OUT=0;CHART0%ANG_IN=0; CHART0%ANG_OUT=0;

end subroutine zero_CHART0

subroutine zero_MAGL0
implicit none
 
    MAGL0%LC_LD_B0_P0=0  ! LC LD B0 P0C
    MAGL0%TILTD_EDGE=0   ! TILTD EDGE
    MAGL0%KIN_KEX_BENDFRINGE_EXACT=.false. ! KILL_ENT_FRINGE, KILL_EXI_FRINGE, bend_fringe,EXACT
	MAGL0%METHOD_NST_NMUL_permfringe_highest=0  ! METHOD,NST,NMUL,permfringr, highest_fringe
    MAGL0%kill_spin=.false. 

end subroutine zero_MAGL0

subroutine zero_patch0
implicit none
 
 
     patch0%A_X1=0;patch0%A_X2=0;patch0%B_X1=0;patch0%B_X2=0
     patch0%A_D=0
     patch0%B_D=0
     patch0%A_ANG=0
     patch0%B_ANG=0 
     patch0%A_T=0
     patch0%B_T=0
     patch0%A_L=0
     patch0%B_L=0
     patch0%ENERGY=0
     patch0%TIME=0
     patch0%GEOMETRY=0
     patch0%track=.true.

end subroutine zero_patch0

subroutine specify_element_type(f,mf)
implicit none
type(fibre), target ::f
integer mf

if(associated(f%mag%D0)) then

Write(mf,*) f%mag%kind, "TYPE(DRIFT1), POINTER :: D0"
return
endif
if(associated(f%mag%K2)) then

Write(mf,*) f%mag%kind,"TYPE(DKD2), POINTER :: K2 "
return
endif
if(associated(f%mag%K3)) then

Write(mf,*) f%mag%kind,"TYPE(KICKT3), POINTER :: K3"
return
endif

if(associated(f%mag%C4)) then

Write(mf,*) f%mag%kind,"TYPE(CAV4), POINTER :: C4 "
return
endif

if(associated(f%mag%S5)) then

Write(mf,*) f%mag%kind,"TYPE(SOL5), POINTER :: S5 "
return
endif

if(associated(f%mag%T6)) then

Write(mf,*) f%mag%kind,"TYPE(KTK), POINTER :: T6 "
return

endif
if(associated(f%mag%T6)) then

Write(mf,*) f%mag%kind,"TYPE(KTKP), POINTER :: T7 "
return
endif

if(associated(f%mag%T7)) then

Write(mf,*) f%mag%kind,"TYPE(TKTF), POINTER :: T7 "
return
endif

if(associated(f%mag%S8)) then

Write(mf,*) f%mag%kind,"TYPE(NSMI), POINTER :: S8"
return
endif

if(associated(f%mag%S9)) then

Write(mf,*) f%mag%kind,"TYPE(SSMI), POINTER :: S9"
return
endif

if(associated(f%mag%TP10)) then

Write(mf,*) f%mag%kind,"TYPE(TEAPOT), POINTER :: TP10 "
return
endif

if(associated(f%mag%MON14)) then

Write(mf,*) f%mag%kind,"TYPE(MON), POINTER :: MON14"
return
endif

if(associated(f%mag%SEP15)) then

Write(mf,*) f%mag%kind,"TYPE(ESEPTUM), POINTER :: SEP15"
return
endif

if(associated(f%mag%K16)) then

Write(mf,*) f%mag%kind,"TYPE(STREX), POINTER :: K16 "
return
endif

if(associated(f%mag%ENGE17)) then

Write(mf,*) f%mag%kind,"TYPE(ENGE), POINTER :: ENGE17"
return
endif

if(associated(f%mag%RCOL18)) then

Write(mf,*) f%mag%kind,"TYPE(RCOL), POINTER :: RCOL18"
return
endif

if(associated(f%mag%ECOL19)) then

Write(mf,*) f%mag%kind,"TYPE(ECOL), POINTER :: ECOL19"
return
endif

if(associated(f%mag%CAV21)) then

Write(mf,*) f%mag%kind,"TYPE(CAV_TRAV), POINTER :: CAV21"
return
endif

if(associated(f%mag%WI)) then

Write(mf,*) f%mag%kind,"TYPE(CAV_TRAV), POINTER :: WI"
return
endif
 
if(associated(f%mag%PA)) then

Write(mf,*) f%mag%kind,"TYPE(PANCAKE), POINTER :: PA"
return
endif

if(associated(f%mag%HE22)) then

Write(mf,*) f%mag%kind,"TYPE(HELICAL_DIPOLE), POINTER :: HE22"
return
endif

if(associated(f%mag%SDR)) then

Write(mf,*) f%mag%kind,"TYPE(HELICAL_DIPOLE), POINTER :: SDR"
return
endif
 
end subroutine specify_element_type 

end module madx_keywords

