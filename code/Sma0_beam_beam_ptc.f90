module beam_beam_ptc
  !  use orbit_ptc
  use S_TRACKING
  ! use madx_ptc_module
  implicit none
  public
  private BBKICKP, BBKICKR, BBKICKnP, BBKICKnR,PATCH_BBR,PATCH_BBP
  private ccperrfP, ccperrfr ,ccperrf,BBKICKxR,BBKICKxP,get_field_BBr,get_field_BBp
  private PATCH_BBxR,PATCH_BBxP,BBPATCH_SPINr,BBPATCH_SPINp
  private rot_spin_xr,rot_spin_xp,rot_spin_zr,rot_spin_zp  !,rot_spin_z,rot_spin_x,
  private rot_spin_yr,rot_spin_yp,get_omega_spin_bb_r,get_omega_spin_bb_p,RAD_SPIN_bb_PROBEr,RAD_SPIN_bb_PROBEp
    !,rot_spin_y 
  !  private TRACK_NODE_LAYOUT_FLAG_R,TRACK_NODE_LAYOUT_FLAG_P
  private imax
  integer ::imax=1000

  INTERFACE ccperrf
     MODULE PROCEDURE ccperrfP
     MODULE PROCEDURE ccperrfr
  END INTERFACE

  INTERFACE BBKICKx
     MODULE PROCEDURE BBKICKxR
     MODULE PROCEDURE BBKICKxP
  END INTERFACE

  INTERFACE BBKICKn
     MODULE PROCEDURE BBKICKnR
     MODULE PROCEDURE BBKICKnP
  END INTERFACE


  INTERFACE BBKICK
     MODULE PROCEDURE BBKICKR
     MODULE PROCEDURE BBKICKP
  END INTERFACE

  INTERFACE PATCH_BB
     MODULE PROCEDURE PATCH_BBR
     MODULE PROCEDURE PATCH_BBP
  END INTERFACE

  INTERFACE BBPATCH_SPIN
     MODULE PROCEDURE BBPATCH_SPINr
     MODULE PROCEDURE BBPATCH_SPINp
  END INTERFACE

  INTERFACE RAD_SPIN_bb_PROBE
     MODULE PROCEDURE RAD_SPIN_bb_PROBEr
     MODULE PROCEDURE RAD_SPIN_bb_PROBEp
  END INTERFACE


 

  INTERFACE PATCH_BBx
     MODULE PROCEDURE PATCH_BBxR
     MODULE PROCEDURE PATCH_BBxP
  END INTERFACE


  INTERFACE get_field_BB
     MODULE PROCEDURE get_field_BBr
     MODULE PROCEDURE get_field_BBp
  END INTERFACE


  INTERFACE rot_spin_x
     MODULE PROCEDURE rot_spin_xr
     MODULE PROCEDURE rot_spin_xp
  END INTERFACE

  INTERFACE rot_spin_y
     MODULE PROCEDURE rot_spin_yr
     MODULE PROCEDURE rot_spin_yp
  END INTERFACE

  INTERFACE rot_spin_z
     MODULE PROCEDURE rot_spin_zr
     MODULE PROCEDURE rot_spin_zp
  END INTERFACE

  INTERFACE get_omega_spin_bb
     MODULE PROCEDURE get_omega_spin_bb_r
     MODULE PROCEDURE get_omega_spin_bb_p
  END INTERFACE



  logical(lp), target :: do_beam_beam= my_false

contains


  subroutine rot_spin_yr(P,ang)
    implicit none
    TYPE(PROBE),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si,st
    INTEGER I
    type(quaternion) dq
    if(p%use_q) then
     dq%x(0)=COS(ang/2)
     dq%x(2)=sin(ang/2)
     dq%x(1)=0
     dq%x(3)=0
     p%q=dq*p%q
    else
    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(1)+SI *P%S(I)%X(3)
       P%S(I)%X(3)=CO *P%S(I)%X(3)-SI *P%S(I)%X(1)
       P%S(I)%X(1)=ST
    ENDDO
    endif
  END subroutine rot_spin_yr

  subroutine rot_spin_Xr(P,ang)
    implicit none
    TYPE(PROBE),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si,st
    INTEGER I
    type(quaternion) dq

    if(p%use_q) then

     dq%x(0)=COS(ang/2)
     dq%x(1)=-sin(ang/2)
     dq%x(2)=0
     dq%x(3)=0
     p%q=dq*p%q

    else
    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(2)+SI *P%S(I)%X(3)
       P%S(I)%X(3)=CO *P%S(I)%X(3)-SI *P%S(I)%X(2)
       P%S(I)%X(2)=ST
    ENDDO
   endif
  END subroutine rot_spin_Xr

  subroutine rot_spin_zr(P,ang)
    implicit none
    TYPE(PROBE),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si,st
    INTEGER I
    type(quaternion) dq

    if(p%use_q) then
  
     dq%x(0)=COS(ang/2)
     dq%x(3)=-sin(ang/2)
     dq%x(1)=0
     dq%x(2)=0
     p%q=dq*p%q
    else
    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(1)+SI *P%S(I)%X(2)
       P%S(I)%X(2)=CO *P%S(I)%X(2)-SI *P%S(I)%X(1)
       P%S(I)%X(1)=ST
    ENDDO
    endif

  END subroutine rot_spin_zr


  subroutine rot_spin_yp(P,ang)
    implicit none
    type(PROBE_8),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si
    type(real_8) st
    !type(real_8) co,si,st
    INTEGER I
    type(quaternion_8) dq
    if(p%use_q) then
     call alloc(dq)
     dq%x(0)=COS(ang/2)
     dq%x(2)=sin(ang/2)
     dq%x(1)=0.0_dp
     dq%x(3)=0.0_dp
     p%q=dq*p%q
     call kill(dq)
    else
    call alloc(st)

    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0r,ISPIN1r
       ST=  CO *P%S(I)%X(1)+SI *P%S(I)%X(3)
       P%S(I)%X(3)=CO *P%S(I)%X(3)-SI *P%S(I)%X(1)
       P%S(I)%X(1)=ST
    ENDDO

    call kill(st)
    endif
  END subroutine rot_spin_yp

  subroutine rot_spin_xp(P,ang)
    implicit none
    type(PROBE_8),INTENT(INOUT) :: P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si
    type(real_8) st
    INTEGER I
    type(quaternion_8) dq

    if(p%use_q) then
     call alloc(dq)
     dq%x(0)=COS(ang/2)
     dq%x(1)=-sin(ang/2)
     dq%x(2)=0.0_dp
     dq%x(3)=0.0_dp
     p%q=dq*p%q
     call kill(dq)
    else
    call alloc(st)

    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(2)+SI *P%S(I)%X(3)
       P%S(I)%X(3)=CO *P%S(I)%X(3)-SI *P%S(I)%X(2)
       P%S(I)%X(2)=ST
    ENDDO

    call kill(st)
    endif
  END subroutine rot_spin_xp

  subroutine rot_spin_zp(P,ang)
    implicit none
    type(PROBE_8),INTENT(INOUT) ::  P
    REAL(DP), INTENT(IN) :: ang
    REAL(DP) co,si
    type(real_8) st
    INTEGER I

    type(quaternion_8) dq

    if(p%use_q) then
     call alloc(dq)
     dq%x(0)=COS(ang/2)
     dq%x(3)=-sin(ang/2)
     dq%x(1)=0.0_dp
     dq%x(2)=0.0_dp
     p%q=dq*p%q
     call kill(dq)
    else
    call alloc(st)

    CO =COS(ang)
    SI =sin(ang)

    DO I=ISPIN0R,ISPIN1R
       ST=  CO *P%S(I)%X(1)+SI *P%S(I)%X(2)
       P%S(I)%X(2)=CO *P%S(I)%X(2)-SI *P%S(I)%X(1)
       P%S(I)%X(1)=ST
    ENDDO

    call kill(st)
    endif
  END subroutine rot_spin_zp


  SUBROUTINE PATCH_BBR(Beam,p,k,BETA0,exact,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(BEAM_BEAM_NODE),TARGET,INTENT(INOUT):: Beam
    type(probe), INTENT(INOUT):: p
    logical(lp),INTENT(IN):: exact,ENTERING
    REAL(DP),INTENT(IN):: BETA0
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    call PATCH_BBx(Beam,p%x,k,BETA0,exact,ENTERING)
    if(k%spin) call BBPATCH_SPIN(Beam,p,ENTERING)

  end SUBROUTINE PATCH_BBR

  SUBROUTINE PATCH_BBP(Beam,P,k,BETA0,exact,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(BEAM_BEAM_NODE),TARGET,INTENT(INOUT):: Beam
    type(probe_8), INTENT(INOUT):: p
    logical(lp),INTENT(IN):: exact,ENTERING
    REAL(DP),INTENT(IN):: BETA0
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K
 
    call PATCH_BBx(Beam,p%x,k,BETA0,exact,ENTERING)
    if(k%spin) call BBPATCH_SPIN(Beam,p,ENTERING)


  END SUBROUTINE PATCH_BBP



  SUBROUTINE BBPATCH_SPINR(Beam,P,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(BEAM_BEAM_NODE),TARGET,INTENT(INOUT):: Beam
    TYPE(PROBE), INTENT(INOUT):: P
    !    real(dp), INTENT(INOUT):: s(3)
    logical(lp),INTENT(IN):: ENTERING
    real(dp) da,A(3)

     IF(ENTERING) THEN
       da=Beam%A(1)+((Beam%A_X1-1.0_dp)/2)*pi
       call rot_spin_x(P,da)
       call rot_spin_y(P,Beam%A(2)) ! 2016_5_9
       call rot_spin_z(P,Beam%A(3))
       da=((Beam%A_X2-1.0_dp)/2.0_dp)*pi
       call rot_spin_x(P,da)
    ELSE
       A=-Beam%A
       da=A(1)+((Beam%A_X2-1.0_dp)/2.0_dp)*pi
       call rot_spin_x(P,da)
       ! error etienne
       !      call rot_spin_y(P,C%PATCH%A_ANG(2))
       call rot_spin_y(P,A(2)) ! 2016_5_9
       call rot_spin_z(P,A(3))
       da=((Beam%A_X1-1.0_dp)/2.0_dp)*pi
       call rot_spin_x(P,da)
    ENDIF
 
 
  END SUBROUTINE BBPATCH_SPINR

  SUBROUTINE BBPATCH_SPINp(Beam,P,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(BEAM_BEAM_NODE),TARGET,INTENT(INOUT):: Beam
    TYPE(PROBE_8), INTENT(INOUT):: P
    !    real(dp), INTENT(INOUT):: s(3)
    logical(lp),INTENT(IN):: ENTERING
    real(dp) da,A(3)

     IF(ENTERING) THEN
       da=Beam%A(1)+((Beam%A_X1-1.0_dp)/2)*pi
       call rot_spin_x(P,da)
       call rot_spin_y(P,Beam%A(2)) ! 2016_5_9
       call rot_spin_z(P,Beam%A(3))
       da=((Beam%A_X2-1.0_dp)/2.0_dp)*pi
       call rot_spin_x(P,da)
    ELSE
       A=-Beam%A
       da=A(1)+((Beam%A_X2-1.0_dp)/2.0_dp)*pi
       call rot_spin_x(P,da)
       ! error etienne
       !      call rot_spin_y(P,C%PATCH%A_ANG(2))
       call rot_spin_y(P,A(2)) ! 2016_5_9
       call rot_spin_z(P,A(3))
       da=((Beam%A_X1-1.0_dp)/2.0_dp)*pi
       call rot_spin_x(P,da)
    ENDIF
 
  END SUBROUTINE BBPATCH_SPINp


  SUBROUTINE PATCH_BBxR(Beam,X,k,BETA0,exact,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(BEAM_BEAM_NODE),TARGET,INTENT(INOUT):: Beam
    real(dp), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: exact,ENTERING
    REAL(DP),INTENT(IN):: BETA0
    REAL(DP) A(3),D(3)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(ENTERING) THEN
       X(3)=Beam%A_X1*X(3);X(4)=Beam%A_X1*X(4);
       CALL ROT_YZ(Beam%A(1),X,BETA0,exact,k%TIME)
       CALL ROT_XZ(Beam%A(2),X,BETA0,exact,k%TIME)
       CALL ROT_XY(Beam%A(3),X)  !,exact)
       CALL TRANS(Beam%D,X,BETA0,exact,k%TIME)
       X(3)=Beam%A_X2*X(3);X(4)=Beam%A_X2*X(4);
    ELSE
       A=-Beam%A
       D=-Beam%D
       X(3)=Beam%A_X2*X(3);X(4)=Beam%A_X2*X(4);
       CALL TRANS(D,X,BETA0,exact,k%TIME)
       CALL ROT_XY(A(3),X)  !,exact)
       CALL ROT_XZ(A(2),X,BETA0,exact,k%TIME)
       CALL ROT_YZ(A(1),X,BETA0,exact,k%TIME)
       X(3)=Beam%A_X1*X(3);X(4)=Beam%A_X1*X(4);
    ENDIF


  END SUBROUTINE PATCH_BBxR

  SUBROUTINE PATCH_BBxP(Beam,X,k,BETA0,exact,ENTERING)
    implicit none
    ! MISALIGNS REAL FIBRES IN PTC ORDER FOR FORWARD AND BACKWARD FIBRES
    TYPE(BEAM_BEAM_NODE),TARGET,INTENT(INOUT):: Beam
    TYPE(REAL_8), INTENT(INOUT):: X(6)
    logical(lp),INTENT(IN):: exact,ENTERING
    REAL(DP),INTENT(IN):: BETA0
    REAL(DP) A(3),D(3)
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(ENTERING) THEN
       X(3)=Beam%A_X1*X(3);X(4)=Beam%A_X1*X(4);
       CALL ROT_YZ(Beam%A(1),X,BETA0,exact,k%TIME)
       CALL ROT_XZ(Beam%A(2),X,BETA0,exact,k%TIME)
       CALL ROT_XY(Beam%A(3),X)  !,exact)
       CALL TRANS(Beam%D,X,BETA0,exact,k%TIME)
       X(3)=Beam%A_X2*X(3);X(4)=Beam%A_X2*X(4);
    ELSE
       A=-Beam%A
       D=-Beam%D
       X(3)=Beam%A_X2*X(3);X(4)=Beam%A_X2*X(4);
       CALL TRANS(D,X,BETA0,exact,k%TIME)
       CALL ROT_XY(A(3),X)  !,exact)
       CALL ROT_XZ(A(2),X,BETA0,exact,k%TIME)
       CALL ROT_YZ(A(1),X,BETA0,exact,k%TIME)
       X(3)=Beam%A_X1*X(3);X(4)=Beam%A_X1*X(4);
    ENDIF


  END SUBROUTINE PATCH_BBxP

 
  subroutine BBKICKR(c,P,k,beta0,exact,time)
  implicit none
  TYPE(integration_node), target::c
  TYPE(BEAM_BEAM_NODE), pointer::Beam
  type(probe), INTENT(INOUT) :: P
  logical, intent(in) ::  exact,time
  real(dp), intent(in) :: beta0
  real(dp) d(3),lh
  real(dp) kicks(3) 
  type(internal_state) k
  integer i
  Beam=> c%bb
  
  if(Beam%n>1) then
  d=0
   lh=(Beam%s(Beam%n)-Beam%s(1))/2
   d(3)=-lh
       CALL TRANS(d,p%X,BETA0,exact,TIME)
       call BBKICKn(Beam,P,1,kicks)
       call RAD_SPIN_bb_PROBE(c,p,k,kicks)
     do i=2,Beam%n
       d(3)=(Beam%s(i)-Beam%s(i-1))    !/2

       CALL TRANS(d,p%X,BETA0,exact,TIME)
       call BBKICKn(Beam,P,i,kicks)
       call RAD_SPIN_bb_PROBE(c,p,k,kicks)

     enddo
   d(3)=-lh
       CALL TRANS(d,p%X,BETA0,exact,TIME)
  else
       call BBKICKn(Beam,P,1,kicks)
       call RAD_SPIN_bb_PROBE(c,p,k,kicks)

  endif
  end subroutine BBKICKR

  subroutine BBKICKP(c,P,k,beta0,exact,time)
  TYPE(integration_node), target::c
  TYPE(BEAM_BEAM_NODE), pointer::Beam
  type(probe_8), INTENT(INOUT) :: P
  logical, intent(in) ::  exact,time
  real(dp), intent(in) :: beta0
  real(dp) d(3),lh,dh
  type(real_8)  kicks(3)
  type(internal_state) k
  integer i
  Beam=> c%bb
  call alloc(kicks)
  if(Beam%n>1) then
  d=0
   lh=(Beam%s(Beam%n)-Beam%s(1))/2
   d(3)=-lh
       CALL TRANS(d,p%X,BETA0,exact,TIME)
       call BBKICKn(Beam,P,1,kicks)
       call RAD_SPIN_bb_PROBE(c,p,k,kicks)

     do i=2,Beam%n
       d(3)=(Beam%s(i)-Beam%s(i-1))  !/2

       CALL TRANS(d,p%X,BETA0,exact,TIME)
       call BBKICKn(Beam,P,i,kicks)
       call RAD_SPIN_bb_PROBE(c,p,k,kicks)

     enddo
   d(3)=-lh
       CALL TRANS(d,p%X,BETA0,exact,TIME)
  else
       call BBKICKn(Beam,P,1,kicks)
       call RAD_SPIN_bb_PROBE(c,p,k,kicks)

  endif
  call kill(kicks)
  end subroutine BBKICKP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_field_BBr(EL,KICK,B,E,X,k)
    implicit none
    TYPE(ELEMENT), POINTER::EL

    REAL(DP),INTENT(INOUT) :: B(3),E(3),KICK(3)
    REAL(DP),INTENT(INOUT) :: X(6)

    REAL(DP) beta0,beta,ff !,p0c
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

 
    beta0=el%PARENT_FIBRE%beta0
  !  p0c=el%parent_fibre%mag%p%p0c
 

    beta=root(1.0_dp+2.0_dp*x(5)/beta0+x(5)**2)/(1.0_dp/BETA0 + x(5))  ! replaced
    ff=1.0_dp+beta**2

 !!!  
       E(1)=kick(1)/ff/EL%P%CHARGE 
       E(2)=kick(2)/ff/EL%P%CHARGE  
       E(3)=0.0_dp


       B(1)= E(2)*beta
       B(2)=-E(1)*beta
       B(3)=0.0_dp


  end subroutine get_field_BBr

  subroutine get_field_BBp(EL,KICK,B,E,X,k)
    implicit none
    TYPE(ELEMENTp), POINTER::EL

    type(real_8),INTENT(INOUT) :: B(3),E(3),KICK(3)
    type(real_8),INTENT(INOUT) :: X(6)

    REAL(DP) beta0 !,p0c
    type(real_8) beta,ff
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

     call alloc(beta,ff)
 
    beta0=el%PARENT_FIBRE%beta0
 !   p0c=el%parent_fibre%mag%p%p0c
 

      beta=sqrt(1.0_dp+2.0_dp*x(5)/beta0+x(5)**2)/(1.0_dp/BETA0 + x(5))  ! replaced
    ff=1.0_dp+beta**2

 !!!  
       E(1)=kick(1)/ff/EL%P%CHARGE 
       E(2)=kick(2)/ff/EL%P%CHARGE  
       E(3)=0.0_dp


       B(1)= E(2)*beta
       B(2)=-E(1)*beta
       B(3)=0.0_dp

     call kill(beta,ff)

  end subroutine get_field_BBp


  subroutine get_omega_spin_bb_r(c,kick,OM,B2,dlds,XP,X,k,ED,B)
    implicit none
    TYPE(integration_node), target::c
    TYPE(ELEMENT), POINTER::EL
    TYPE(MAGNET_CHART), POINTER::P
    REAL(DP),INTENT(INOUT) :: X(6),OM(3),B2,XP(2),DLDS,B(3),ED(3),kick(3)
    REAL(DP)  BPA(3),BPE(3),D1,D2,GAMMA,EB(3),EFD(3),beta,e(3)
    REAL(DP) BETA0,GAMMA0I,XPA(2),del  !,phi,z
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(.NOT.CHECK_STABLE) return
    el=>c%parent_fibre%mag
    P=>EL%P
    P%DIR    => C%PARENT_FIBRE%DIR
    P%beta0  => C%PARENT_FIBRE%beta0
    P%GAMMA0I=> C%PARENT_FIBRE%GAMMA0I
    P%GAMBET => C%PARENT_FIBRE%GAMBET
    P%MASS => C%PARENT_FIBRE%MASS
    P%ag => C%PARENT_FIBRE%ag
    P%CHARGE=>C%PARENT_FIBRE%CHARGE
    ! DLDS IS  REALLY D(CT)/DS * (1/(ONE/BETA0+X(5)))
    OM(2)=0.0_dp
    EB=0.0_dp
    BPA=0.0_dp
    BPE=0.0_dp
    B=0.0_dp
    E=0.0_dp
    ED=0.0_dp
    EFD=0.0_dp
  !  phi=0.0_dp

    xp(1)=x(2)
    xp(2)=x(4)   !  to prevent a crash in monitors, etc... CERN june 2010
    dlds=0.0_dp
    del=x(5)
    CALL get_field_BB(EL,KICK,B,E,X,k)

 
    CALL B_PARA_PERP(k,EL,X,B,BPA,BPE,XP,XPA,ed,E,EB,EFD)

 

       IF(k%TIME) THEN
          DLDS=1.0_dp/root(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/root((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF


    IF(.not.k%TIME) THEN
      del=(2*del+del**2)/(root(1.0_dp/p%beta0**2+2.0_dp*del+del**2)+1.0_dp/p%beta0)
    endif

    !  MUST ALWAYS COMPUTER GAMMA EVEN IF TIME=FALSE.
    GAMMA=P%BETA0/P%GAMMA0I*( 1.0_dp/P%BETA0 + del )

    OM(1)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(1) + (1.0_dp+p%AG)*BPA(1) )
    OM(2)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(2) + (1.0_dp+p%AG)*BPA(2) )+OM(2)
    OM(3)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(3) + (1.0_dp+p%AG)*BPA(3) )


    beta=root(1.0_dp+2.0_dp*del/p%beta0+del**2)/(1.0_dp/P%BETA0 + del)  ! replaced


    DO I=1,3
       OM(I)=OM(I)+a_spin_scale*DLDS*beta*gamma*(p%AG+1.0_dp/(1.0_dp+GAMMA))*EB(I)
    ENDDO

   beta=root(1.0_dp+2.0_dp*x(5)/p%beta0+x(5)**2)*P%BETA0/P%GAMMA0I  ! replace  this

    om(1)=-DLDS*0.5_dp*e_muon*beta*(ed(2)*BPE(3)-ed(3)*BPE(2)) +  om(1)
    om(2)=-DLDS*0.5_dp*e_muon*beta*(ed(3)*BPE(1)-ed(1)*BPE(3)) +  om(2)
    om(3)=-DLDS*0.5_dp*e_muon*beta*(ed(1)*BPE(2)-ed(2)*BPE(1)) +  om(3)

    DO I=1,3
       OM(I)=OM(I)-DLDS*0.5_dp*e_muon*(GAMMA*E(I)+(1-GAMMA)*EFD(I))
    ENDDO

    !IF(.not.k%TIME) THEN
    !   del=(2.0_dp*del/p%beta0+del**2)/(sqrt(1.0_dp+2.0_dp*del/p%beta0+del**2)+1.0_dp)
    !endif

    if((k%radiation.or.k%envelope)) then
       !      if(P%RADIATION) then
       B2=BPE(1)**2+BPE(2)**2+BPE(3)**2
       !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    ENDIF

  end subroutine get_omega_spin_bb_r


  subroutine get_omega_spin_bb_p(c,kick,OM,B2,dlds,XP,X,k,ED,B)
    implicit none
    TYPE(integration_node), target::c
    TYPE(ELEMENTp), POINTER::EL
    TYPE(MAGNET_CHART), POINTER::P
    TYPE(REAL_8), INTENT(INOUT) :: X(6),OM(3),B2,XP(2),B(3),ED(3),kick(3)
    TYPE(REAL_8)  BPA(3),BPE(3),DLDS,D1,D2,GAMMA,EB(3),efd(3),XPA(2),e(3),beta,del !,phi,z
    REAL(DP) BETA0,GAMMA0I
    INTEGER I
    TYPE(INTERNAL_STATE) k !,OPTIONAL :: K

    IF(.NOT.CHECK_STABLE) return
    CALL ALLOC(E,3)
    CALL ALLOC(efd,3)
    CALL ALLOC(beta,del)
    CALL ALLOC(EB,3)
    CALL ALLOC(BPA,3)
    CALL ALLOC(BPE,3)
    CALL ALLOC(XPA,2)
    CALL ALLOC(D1,D2,GAMMA)

    el=>c%parent_fibre%magp

    P=>EL%P
    P%DIR    => C%PARENT_FIBRE%DIR
    P%beta0  => C%PARENT_FIBRE%beta0
    P%GAMMA0I=> C%PARENT_FIBRE%GAMMA0I
    P%GAMBET => C%PARENT_FIBRE%GAMBET
    P%MASS => C%PARENT_FIBRE%MASS
    P%ag => C%PARENT_FIBRE%ag
    P%CHARGE=>C%PARENT_FIBRE%CHARGE
    ! DLDS IS  REALLY D(CT)/DS * (1/(ONE/BETA0+X(5)))
    OM(2)=0.0_dp
    EB=0.0_dp
    BPA=0.0_dp
    BPE=0.0_dp
    B=0.0_dp
    E=0.0_dp
    ED=0.0_dp
    EFD=0.0_dp
  !  phi=0.0_dp

    xp(1)=x(2)
    xp(2)=x(4)   !  to prevent a crash in monitors, etc... CERN june 2010
    dlds=0.0_dp
    del=x(5)
    CALL get_field_BB(EL,KICK,B,E,X,k)

 
    CALL B_PARA_PERP(k,EL,X,B,BPA,BPE,XP,XPA,ed,E,EB,EFD)

 

       IF(k%TIME) THEN
          DLDS=1.0_dp/sqrt(1.0_dp+2.0_dp*del/P%BETA0+del**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ELSE
          DLDS=1.0_dp/sqrt((1.0_dp+del)**2-XPA(2)**2-XPA(1)**2)*(1.0_dp+P%b0*X(1))
       ENDIF


    IF(.not.k%TIME) THEN
      del=(2*del+del**2)/(sqrt(1.0_dp/p%beta0**2+2.0_dp*del+del**2)+1.0_dp/p%beta0)
    endif

    !  MUST ALWAYS COMPUTER GAMMA EVEN IF TIME=FALSE.
    GAMMA=P%BETA0/P%GAMMA0I*( 1.0_dp/P%BETA0 + del )

    OM(1)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(1) + (1.0_dp+p%AG)*BPA(1) )
    OM(2)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(2) + (1.0_dp+p%AG)*BPA(2) )+OM(2)
    OM(3)=-DLDS*a_spin_scale*( (1.0_dp+p%AG*GAMMA)*BPE(3) + (1.0_dp+p%AG)*BPA(3) )


    beta=sqrt(1.0_dp+2.0_dp*del/p%beta0+del**2)/(1.0_dp/P%BETA0 + del)  ! replaced


    DO I=1,3
       OM(I)=OM(I)+a_spin_scale*DLDS*beta*gamma*(p%AG+1.0_dp/(1.0_dp+GAMMA))*EB(I)
    ENDDO

   beta=sqrt(1.0_dp+2.0_dp*x(5)/p%beta0+x(5)**2)*P%BETA0/P%GAMMA0I  ! replace  this

    om(1)=-DLDS*0.5_dp*e_muon*beta*(ed(2)*BPE(3)-ed(3)*BPE(2)) +  om(1)
    om(2)=-DLDS*0.5_dp*e_muon*beta*(ed(3)*BPE(1)-ed(1)*BPE(3)) +  om(2)
    om(3)=-DLDS*0.5_dp*e_muon*beta*(ed(1)*BPE(2)-ed(2)*BPE(1)) +  om(3)

    DO I=1,3
       OM(I)=OM(I)-DLDS*0.5_dp*e_muon*(GAMMA*E(I)+(1-GAMMA)*EFD(I))
    ENDDO

    !IF(.not.k%TIME) THEN
    !   del=(2.0_dp*del/p%beta0+del**2)/(sqrt(1.0_dp+2.0_dp*del/p%beta0+del**2)+1.0_dp)
    !endif

    if((k%radiation.or.k%envelope)) then
       !      if(P%RADIATION) then
       B2=BPE(1)**2+BPE(2)**2+BPE(3)**2
       !        B2=-CRADF(EL%P)*(one+X(5))**3*B2*DLDS
    ENDIF
    CALL KILL(E,3)
    CALL KILL(efd,3)
    CALL KILL(beta,del)
    CALL KILL(EB,3)
    CALL KILL(BPA,3)
    CALL KILL(BPE,3)
    CALL KILL(XPA,2)
    CALL KILL(D1,D2,GAMMA)

  end subroutine get_omega_spin_bb_p


SUBROUTINE RAD_SPIN_bb_PROBER(c,p,k,kick)
    type(probe), INTENT(INOUT) :: p
    TYPE(integration_node),target :: c
    REAL(DP)  B(3),XP(2),XPA(2),ed(3),kick(3)
    REAL(DP) om(3),b2,dlds,FAC
     TYPE(INTERNAL_STATE) k 
 
     integer i
     FAC=0.5_dp
 
 
    call get_omega_spin_bb(c,kick,OM,B2,dlds,XP,p%X,k,ED,B)
 
    !if((k%radiation.or.k%envelope)) then
    !   call radiate_2_probe(c,DS,FAC,P,b2,dlds,k,pos)
    !endif

   if(k%spin) then   
   do i=1,3
     om(i)=om(i)/2.0_dp
    enddo
   
    
   call push_quaternion(p,om)
  endif
 !   if((k%radiation.or.k%envelope)) then
 !call radiate_2_probe(c,DS,FAC,P,b2,dlds,k,pos)
!
!    endif

 end SUBROUTINE RAD_SPIN_bb_PROBER


SUBROUTINE RAD_SPIN_bb_PROBEp(c,p,k,kick)
    type(probe_8), INTENT(INOUT) :: p
    TYPE(integration_node),target :: c
    TYPE(REAL_8)  B(3),XP(2),XPA(2),ed(3),om(3),b2,dlds,kick(3)
    REAL(DP)  FAC
     TYPE(INTERNAL_STATE) k 
     integer i
      CALL alloc(B);CALL alloc(XP);CALL alloc(XPA);CALL alloc(ed);
     CALL alloc(OM);CALL alloc(B2);CALL alloc(DLDS); 
 
 
    call get_omega_spin_bb(c,kick,OM,B2,dlds,XP,p%X,k,ED,B)
 
    !if((k%radiation.or.k%envelope)) then
    !   call radiate_2_probe(c,DS,FAC,P,b2,dlds,k,pos)
    !endif

   if(k%spin) then   
   do i=1,3
     om(i)=om(i)/2.0_dp
    enddo
   
    
   call push_quaternion(p,om)
  endif
 !   if((k%radiation.or.k%envelope)) then
 !call radiate_2_probe(c,DS,FAC,P,b2,dlds,k,pos)
!
!    endif
      CALL kill(B);CALL kill(XP);CALL kill(XPA);CALL kill(ed);
     CALL kill(OM);CALL kill(B2);CALL kill(DLDS); 

 end SUBROUTINE RAD_SPIN_bb_PROBEp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine BBKICKnR(BB,P,n,kick)

    implicit none
    !----------------------------------------------------------------------*
    ! purpose:                                                             *
    !   track a set of particle through a beam-beam interaction region.    *
    !   see mad physicist's manual for the formulas used.                  *
    !input:                                                                *
    ! input/output:                                                        *
    !   b%x(:,6)(double)  track coordinates: (x, px, y, py,  pt,t).        *
    !   b%n    (integer) number of tracks.                                 *
    !----------------------------------------------------------------------*
    real(dp) sx2,sy2,xs,ys,rho2,fk,tk,phix,phiy,rk,xb,yb,crx,cry,xr,yr,r,r2,&
         cbx,cby,ten3m,explim,sx,sy,xm,ym
    TYPE(BEAM_BEAM_NODE), INTENT(INOUT) ::BB
    type(probe), INTENT(INOUT) :: P
    real(dp) X(6)
    real(dp), optional ::  kick(3)
    parameter(ten3m=1.0e-3_dp,explim=150.0_dp)
    integer n

    SX=BB%SX(n)
    SY=BB%SY(n)
    XM=BB%XM(n)
    YM=BB%YM(n)
    FK=BB%FK(n)
    !    write(6,*) "bb%FK = " ,bb%FK
    kick=0.0_dp
    if (fk == 0.0_dp)  return
    x=p%x
    kick(1)=-x(2)   !  negative for substraction later on exit
    kick(2)=-x(4)

    sx2 = sx*sx
    sy2 = sy*sy
    !---- limit formulae for sigma(x) = sigma(y).
    if (abs(sx2 - sy2) .le. ten3m * (sx2 + sy2)) then !ten3m = 1.0d-3
       xs = x(1) - xm
       ys = x(3) - ym
       rho2 = xs * xs + ys * ys
       tk = rho2 / (sx2 + sy2)
       if (tk .gt. explim) then
          phix = xs * fk / rho2
          phiy = ys * fk / rho2
       else if (rho2 .ne. 0.0_dp) then
          phix = xs * fk / rho2 * (1.0_dp - exp(-tk) )
          phiy = ys * fk / rho2 * (1.0_dp - exp(-tk) )
       else
          phix = 0.0_dp
          phiy = 0.0_dp
       endif
       phix = phix - bb%bbk(n,1) ! subtract horizontal bb kick
       phiy = phiy - bb%bbk(n,2) ! subtract vertical co
       x(2) = x(2) + phix
       x(4) = x(4) + phiy

       !---- case sigma(x) > sigma(y).
    else if (sx2 > sy2) then
       r2 = 2.0_dp * (sx2 - sy2)
       r  = sqrt(r2)
       rk = fk * sqrt(pi) / r
       xs = x(1) - xm
       ys = x(3) - ym
       xr = abs(xs) / r
       yr = abs(ys) / r
       call ccperrf(xr, yr, crx, cry)
       tk = (xs * xs / sx2 + ys * ys / sy2) / 2.0_dp
       if (tk .gt. explim) then
          phix = rk * cry
          phiy = rk * crx
       else
          xb = (sy / sx) * xr
          yb = (sx / sy) * yr
          call ccperrf(xb, yb, cbx, cby)
          phix = rk * (cry - exp(-tk) * cby)
          phiy = rk * (crx - exp(-tk) * cbx)
       endif
       x(2) = x(2) + phix * sign(1.0_dp,xs)
       x(4) = x(4) + phiy * sign(1.0_dp,ys)
       x(2) = x(2) - BB%bbk(n,1)
       x(4) = x(4) - BB%bbk(n,2)

       !---- case sigma(x) < sigma(y).
    else
       r2 = 2.0_dp * (sy2 - sx2)
       r  = sqrt(r2)
       rk = fk * sqrt(pi) / r
       !       do itrack = 1, b%n
       !         IF(B%U(itrack)) CYCLE
       !          xs = track(1,itrack) - xm
       !          ys = track(3,itrack) - ym
       xs = x(1) - xm
       ys = x(3) - ym
       xr = abs(xs) / r
       yr = abs(ys) / r

       call ccperrf(yr, xr, cry, crx)


       tk = (xs * xs / sx2 + ys * ys / sy2) / 2.0_dp
       if (tk .gt. explim) then
          phix = rk * cry
          phiy = rk * crx
       else
          xb  = (sy / sx) * xr
          yb  = (sx / sy) * yr
          call ccperrf(yb, xb, cby, cbx)

          phix = rk * (cry - exp(-tk) * cby)
          phiy = rk * (crx - exp(-tk) * cbx)
       endif

       x(2) = x(2) + phix * sign(1.0_dp,xs)
       x(4) = x(4) + phiy * sign(1.0_dp,ys)

       x(2) = x(2) - BB%bbk(n,1)
       x(4) = x(4) - BB%bbk(n,2)

    endif
    !    IF(DZ/=ZERO) CALL DRIFT_BEAM_BACK_TO_POSITION(th,DZ,B)

    kick(1)=x(2)+kick(1)   !  negative for substraction later on exit
    kick(2)=x(4)+kick(2) 

    p%x=x

  end subroutine BBKICKnR

  subroutine ccperrfr(xx, yy, wx, wy)
    implicit none
    !----------------------------------------------------------------------*
    ! purpose:                                                             *
    !   modification of wwerf, double precision complex error function,    *
    !   written at cern by K. Koelbig.                                     *
    ! input:                                                               *
    !   xx, yy    (double)    real + imag argument                         *
    ! output:                                                              *
    !   wx, wy    (double)    real + imag function result                  *
    !----------------------------------------------------------------------*
    integer n,nc,nu
    real(dp), INTENT(INOUT):: xx,yy,wx,wy
    real(dp) x,y,q,h,xl,xh,yh,tx,ty,tn,sx,sy,saux,  &
         rx(33),ry(33),cc,xlim,ylim,fac1,fac2,fac3
    parameter(cc=1.12837916709551_dp,        &
         xlim=5.33_dp,ylim=4.29_dp,fac1=3.2_dp,fac2=23.0_dp,fac3=21.0_dp)

    x = abs(xx)
    y = abs(yy)

    if (y .lt. ylim  .and.  x .lt. xlim) then
       q  = (1.0_dp - y / ylim) * sqrt(1.0_dp - (x/xlim)**2)
       h  = 1.0_dp / (fac1 * q)
       nc = 7 + int(fac2*q)
       xl = h**(1 - nc)
       xh = y + 0.5_dp/h
       yh = x
       nu = 10 + int(fac3*q)
       rx(nu+1) = 0.0_dp
       ry(nu+1) = 0.0_dp

       do n = nu, 1, -1
          tx = xh + n * rx(n+1)
          ty = yh - n * ry(n+1)
          tn = tx*tx + ty*ty
          rx(n) = 0.5_dp * tx / tn
          ry(n) = 0.5_dp * ty / tn
       enddo

       sx = 0.0_dp
       sy = 0.0_dp

       do n = nc, 1, -1
          saux = sx + xl
          sx = rx(n) * saux - ry(n) * sy
          sy = rx(n) * sy + ry(n) * saux
          xl = h * xl
       enddo

       wx = cc * sx
       wy = cc * sy
    else
       xh = y
       yh = x
       rx(1) = 0.0_dp
       ry(1) = 0.0_dp

       do n = 9, 1, -1
          tx = xh + n * rx(1)
          ty = yh - n * ry(1)
          tn = tx*tx + ty*ty
          rx(1) = 0.5_dp * tx / tn
          ry(1) = 0.5_dp * ty / tn
       enddo

       wx = cc * rx(1)
       wy = cc * ry(1)
    endif

    !      if(y .eq. zero) wx = exp(-x**2)
    if(yy .lt. 0.0_dp) then
       wx =   2.0_dp * exp(y*y-x*x) * cos(2.0_dp*x*y) - wx
       wy = - 2.0_dp * exp(y*y-x*x) * sin(2.0_dp*x*y) - wy
       if(xx .gt. 0.0_dp) wy = -wy
    else
       if(xx .lt. 0.0_dp) wy = -wy
    endif

  end SUBROUTINE ccperrfr



  subroutine BBKICKnP(BB,P,n,kick)

    implicit none
    !----------------------------------------------------------------------*
    ! purpose:                                                             *
    !   track a set of particle through a beam-beam interaction region.    *
    !   see mad physicist's manual for the formulas used.                  *
    !input:                                                                *
    ! input/output:                                                        *
    !   b%x(:,6)(double)  track coordinates: (x, px, y, py,  pt,t).      *
    !   b%n    (integer) number of tracks.                              *
    !----------------------------------------------------------------------*
    integer it
    TYPE(REAL_8) xs,ys,tk,phix,phiy,xb,yb,crx,cry
    TYPE(REAL_8) xr,yr,cbx,cby,rho2
    TYPE(REAL_8), optional ::  kick(3)
    REAL(DP) sx2,sy2,sx,sy,xm,ym,fk,ten3m,explim,xn1,xn2,xs1,xs2,arglim,rk
    REAL(DP) r,r2 ,nr
    TYPE(BEAM_BEAM_NODE), INTENT(INOUT) ::BB
    TYPE(REAL_8) X(6)
    TYPE(probe_8), INTENT(INOUT) ::p
    parameter(ten3m=1.0e-3_dp,arglim=1.0e-2_dp,explim=150.0_dp)
    integer n

       kick=0.0_dp
    if (BB%fk(n) == 0.0_dp)  return



    CALL ALLOC(xr,yr,cbx,cby,rho2)
    CALL ALLOC(xs,ys,tk,phix,phiy)
    CALL ALLOC(xb,yb,crx,cry)
    call alloc(x)
    x=p%x
    kick(1)=-x(2)   !  negative for substraction later on exit
    kick(2)=-x(4)

    SX=BB%SX(n)
    SY=BB%SY(n)
    XM=BB%XM(n)
    YM=BB%YM(n)
    FK=BB%FK(n)


    sx2 = sx*sx
    sy2 = sy*sy
    !---- limit formulae for sigma(x) = sigma(y).
    xn1=abs(sx2 - sy2)
    xn2= ten3m * (sx2 + sy2)
    xs1=sx2
    xs2= sy2
    if (xn1 .le.xn2 ) then !ten3m = 1.0d-3
       xs = x(1) - xm
       ys = x(3) - ym
       rho2 = xs * xs + ys * ys
       tk = rho2 / (sx2 + sy2)
       if (tk .gt. explim) then
          phix = xs * fk / rho2
          phiy = ys * fk / rho2
       else if (tk > arglim) then
          phix = xs * fk / rho2 * (1.0_dp - exp(-tk) )
          phiy = ys * fk / rho2 * (1.0_dp - exp(-tk) )
       else

          xr=1.0_dp
          yr=1.0_dp

          nr=mybig
          do it=1,imax
             xr=-xr*tk/(it+1)
             yr=yr+xr
             if(it>10)nr=full_abs(xr)
             if(nr<=puny) exit
          enddo
          if(it>imax-2) then
             write(6,*) it,nr
             write(6,*) " Stopped in Beam-Beam "
          endif
          phix = xs * fk / (2.0_dp * sx2) * YR ! fudge
          phiY = Ys * fk / (2.0_dp * sx2) * YR ! fudge
       endif

       phix = phix - bb%bbk(n,1)
       phiy = phiy - bb%bbk(n,2)

       x(2) = x(2) + phix
       x(4) = x(4) + phiy


       !---- case sigma(x) > sigma(y).
    else


       xs = x(1) - xm
       ys = x(3) - ym
       tk = (xs * xs / sx2 + ys * ys / sy2) / 2.0_dp

       if(xs1 > xs2) then
          r2 = 2.0_dp * (sx2 - sy2)
          r  = sqrt(r2)
          rk = fk * sqrt(pi) / r
          xr = xs / r   !
          yr = ys / r   !
          
!          if( (xr.sub.'0') < 0) then
          if( (xr.sub.0) < 0) then    ! should same but faster perhaps
            xr = -xr
          endif
!          if( (yr.sub.'0') < 0) then
          if( (yr.sub.0) < 0) then   ! should same but faster perhaps
            yr = -yr
          endif
          
          call ccperrf(xr, yr, crx, cry)
          if (tk .gt. explim) then
             phix = rk * cry
             phiy = rk * crx
          else
             xb = (sy / sx) * xr
             yb = (sx / sy) * yr
             call ccperrf(xb, yb, cbx, cby)
             phix = rk * (cry - exp(-tk) * cby)
             phiy = rk * (crx - exp(-tk) * cbx)
          endif
          !          if (.NOT.bborbit)  then
          if(xs<0) then
           x(2) = x(2) - phix
          else
           x(2) = x(2) + phix
          endif
          if(ys<0) then
          x(4) = x(4) - phiy
          else
          x(4) = x(4) + phiy
          endif
          x(2) = x(2) - BB%bbk(n,1)
          x(4) = x(4) - BB%bbk(n,2)
          !          endif
          !       enddo

          !---- case sigma(x) < sigma(y).
       else
          r2 = 2.0_dp * (sy2 - sx2)
          r  = sqrt(r2)
          rk = fk * sqrt(pi) / r

          !       xs = x(1) - xm
          !       ys = x(3) - ym
          xr = xs / r !abs
          yr = ys / r !abs

          if( (xr.sub.'0') < 0) then
            xr = -xr
          endif
          if( (yr.sub.'0') < 0) then
            yr = -yr
          endif


          call ccperrf(yr, xr, cry, crx)

          !       tk = (xs * xs / sx2 + ys * ys / sy2) / two
          if (tk .gt. explim) then
             phix = rk * cry
             phiy = rk * crx
          else
             xb  = (sy / sx) * xr
             yb  = (sx / sy) * yr
             call ccperrf(yb, xb, cby, cbx)
             phix = rk * (cry - exp(-tk) * cby)
             phiy = rk * (crx - exp(-tk) * cbx)
          endif
          !          track(2,itrack) = track(2,itrack) + phix * sign(one,xs)
          !          track(4,itrack) = track(4,itrack) + phiy * sign(one,ys)
          if(xs<0) then
           x(2) = x(2) - phix
          else
           x(2) = x(2) + phix
          endif
          if(ys<0) then
          x(4) = x(4) - phiy
          else
          x(4) = x(4) + phiy
          endif

          !          if (.NOT.bborbit)  then
          !--- subtract closed orbit kick
          !            track(2,itrack) = track(2,itrack) - bb_kick(1,ipos)
          !            track(4,itrack) = track(4,itrack) - bb_kick(2,ipos)
          x(2) = x(2) - BB%bbk(n,1)
          x(4) = x(4) - BB%bbk(n,2)
          !          endif
          !       enddo
       endif
    endif

    kick(1)=x(2)+kick(1)   !  negative for substraction later on exit
    kick(2)=x(4)+kick(2) 

    p%x=x

    CALL KILL(xr,yr,cbx,cby,rho2)
    CALL KILL(xs,ys,tk,phix,phiy)
    CALL KILL(xb,yb,crx,cry)
    call KILL(x)



  end subroutine BBKICKnP

  subroutine ccperrfP(xx, yy, wx, wy)
    implicit none
    !----------------------------------------------------------------------*
    ! purpose:                                                             *
    !   modification of wwerf, double precision complex error function,    *
    !   written at cern by K. Koelbig.                                     *
    ! input:                                                               *
    !   xx, yy    (double)    real + imag argument                         *
    ! output:                                                              *
    !   wx, wy    (double)    real + imag function result                  *
    !----------------------------------------------------------------------*
    TYPE(REAL_8),INTENT(INOUT):: xx,yy,wx,wy
    TYPE(complex_8) z,zt,w
    complex(dp) z0,w0,w1,wt0
    real(dp) xx0, yy0, wx0, wy0
    integer i
    call alloc( z,zt,w)

    z=xx+i_*yy
    z0=z
    z=z-z0

    xx0=real(z0)
    yy0=aimag(z0)
    call ccperrf(xx0, yy0, wx0, wy0)

    w0=wx0+i_*wy0

    w1=-2.0_dp*z0*w0+2.0_dp*i_/sqrt(pi)

    w=w0+w1*z

    zt=z

    do i=2,c_%no

       zt=z*zt
       wt0=  -2.0_dp*(w0+z0*w1)/i

       w=w+wt0*zt

       w0=w1
       w1=wt0

    enddo

    wx=real(w)
    wy=aimag(w)


    call kill( z,zt,w)

  end SUBROUTINE ccperrfP

!!!!!!!!!!!!!!!! for  old X(6) routines no spin !!!!!!!!!!!!!!!!!!!

  subroutine BBKICKxR(BB,X,beta0,exact,time)
  implicit none
  TYPE(BEAM_BEAM_NODE), INTENT(INOUT) ::BB
  real(dp), INTENT(INOUT) :: x(6)
  type(probe) P
  logical, intent(in) ::  exact,time
  real(dp), intent(in) :: beta0
  real(dp) d(3),lh 
  integer i

  p=x

  if(bb%n>1) then
  d=0
   lh=(bb%s(bb%n)-bb%s(1))/2
   d(3)=-lh
       CALL TRANS(d,p%X,BETA0,exact,TIME)
       call BBKICKn(BB,P,1)
     do i=2,bb%n
       d(3)=(bb%s(i)-bb%s(i-1))    !/2

       CALL TRANS(d,p%X,BETA0,exact,TIME)
       call BBKICKn(BB,P,i)
     enddo
   d(3)=-lh
       CALL TRANS(d,p%X,BETA0,exact,TIME)
  else
       call BBKICKn(BB,P,1)
  endif
  x=p%x
  end subroutine BBKICKxR

  subroutine BBKICKxP(BB,x,beta0,exact,time)
  implicit none
  TYPE(BEAM_BEAM_NODE), INTENT(INOUT) ::BB
  type(real_8), INTENT(INOUT) :: x(6)
  type(probe_8) p
  logical, intent(in) ::  exact,time
  real(dp), intent(in) :: beta0
  real(dp) d(3),lh,dh
  integer i
  call alloc(p)
  p%x=x
  if(bb%n>1) then
  d=0
   lh=(bb%s(bb%n)-bb%s(1))/2
   d(3)=-lh
       CALL TRANS(d,p%X,BETA0,exact,TIME)
       call BBKICKn(BB,P,1)
     do i=2,bb%n
       d(3)=(bb%s(i)-bb%s(i-1))  !/2

       CALL TRANS(d,p%X,BETA0,exact,TIME)
       call BBKICKn(BB,P,i)
     enddo
   d(3)=-lh
       CALL TRANS(d,p%X,BETA0,exact,TIME)
  else
       call BBKICKn(BB,P,1)
  endif
  x=p%x
  call kill(p)

  end subroutine BBKICKxP


end module  beam_beam_ptc
