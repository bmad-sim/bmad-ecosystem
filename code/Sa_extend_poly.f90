!The Polymorphic Tracking Code
!Copyright (C) Etienne Forest and CERN

module S_extend_poly
  !USE tree_element_MODULE
   use c_TPSA
   use precision_constants
  IMPLICIT NONE
  public
  logical(lp), target :: ALWAYS_knobs=.false.
  type(real_8) e_muon_scale,a_spin_scale
   private  daddsco,scdaddo,daddsc,scdadd
   private equal_real8_cmap,equal_cmap_real8,EQUAL_c_map_RAY8,EQUAL_RAY8_c_map
  integer,private,parameter::ndd=6
  ! LD: 22.03.2019 (see Sc_euclidean.f90, Sh_def_kinf.f90 and Sr_spin.f90)
  character(len=150) :: ELEM_NAME = "UNKNOWN"
  integer            :: MAPDUMP = 0 ! 0: no dump, 1: dump no=0, 2: dump no=1

  INTERFACE OPERATOR (+)

     MODULE PROCEDURE daddsco   !# c_damap + real(6)
     MODULE PROCEDURE scdaddo   !# real(6) + c_damap
     MODULE PROCEDURE daddsc    !# c_damap + probe_8
     MODULE PROCEDURE scdadd    !# probe_8 + c_damap

  END INTERFACE

  INTERFACE assignment (=)

     MODULE PROCEDURE equal_real8_cmap  !# put cmap in real_8(6)
     MODULE PROCEDURE equal_cmap_real8  ! put real_8(6) in cmap
     MODULE PROCEDURE EQUAL_c_map_RAY8   !#  c_damap=probe_8
     MODULE PROCEDURE EQUAL_RAY8_c_map   !#  probe_8=c_damap

    end INTERFACE
CONTAINS

  
  FUNCTION scdadd( S2,S1  )
    implicit none
    TYPE (probe_8) scdadd
    TYPE (c_damap), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,i,j,nd2h,dc
    type(taylor) d
    integer rf,nd2t,nd2

    call get_rf(rf,nd2t,nd2)

    call alloc(d)

    !   call ass(scdadd)
    scdadd%u=my_false
    scdadd%E_ij=0.0_dp
    scdadd%nac=s2%nac
    scdadd%use_q=s2%use_q  
    if(doing_ac_modulation_in_ptc) then
       dc=2*rf
    else
       dc=0
    endif

    if(c_%ndpt==0) then    ! 1                    
      nd2h=nd2t
       do i=1,nd2h                 !   from 1-4 or 1-6 (if ndpt=0)
          localmaster=master
          call ass(scdadd%x(i))
          !       scdadd%x(i)=s1%m%v(i)+s2%x(i)

          d=s1%V(i)+s2%x(i)-(s1%V(i).sub.'0')
          scdadd%x(i)=d
          master=localmaster
       enddo
  
       do i=nd2h+1,6
          localmaster=master
          call ass(scdadd%x(i))
!  !       if((c_%npara==5+dc).AND.I==5) then   ! npr
          if(((c_%npara==5+dc).AND.I==5+ndpt_bmad).or.((c_%npara==3+dc).AND.I==5+ndpt_bmad)) then   ! npr
             scdadd%x(i)=s2%x(i)+(1.0_dp.mono.c_%npara)
          else
             scdadd%x(i)=s2%x(i)
          endif
          master=localmaster
       enddo

    else        ! 1
      nd2h=nd2t+2
       do i=1,nd2h                 !   from 1-4 or 1-6 (if ndpt=0)
          localmaster=master
          call ass(scdadd%x(i))
          !       scdadd%x(i)=s1%m%v(i)+s2%x(i)

          d=s1%V(i)+s2%x(i)-(s1%V(i).sub.'0')
          scdadd%x(i)=d
          master=localmaster
       enddo
    endif       ! 1


 do i=1,scdadd%nac
    localmaster=master
    call ass(scdadd%AC(i)%x(1))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
     j=2*scdadd%nac-(2*i-1)
     d=s1%V(C_%ND2-j) +addclock*s2%AC(i)%x(1)
    scdadd%ac(i)%x(1)=d
    master=localmaster
    localmaster=master
     j=2*scdadd%nac-(2*i)
    call ass(scdadd%AC(i)%x(2))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    d=s1%V(C_%ND2-j) +addclock*s2%AC(i)%x(2)
    scdadd%ac(i)%x(2)=d
    master=localmaster
    localmaster=master
    call ass(scdadd%AC(i)%om)
!    call ass(scdadd%AC%t)
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    scdadd%AC(i)%om=s2%AC(i)%om
    scdadd%AC(i)%t=s2%AC(i)%t
!    scdadd%AC(i)%f=s2%AC(i)%f
!    scdadd%AC(i)%phase=s2%AC(i)%phase
    master=localmaster
enddo
    !    endif


 !   if(use_quaternion)   THEN
       DO J=0,3
          localmaster=master
          call ass(scdadd%q%x(J))
          d=S1%q%x(j)
          scdadd%q%x(J)=d
          master=localmaster
       ENDDO
!else
    DO I=1,3
       !          call ass(scdadd%s%x(i))
       DO J=1,3
          localmaster=master
          call ass(scdadd%s(J)%x(i))
          d=S1%S%s(I,J)
          scdadd%s(J)%x(i)=d
          master=localmaster
       ENDDO

    ENDDO

!endif
!    scdadd%damps=s1%damps
!    scdadd%d_spin=s1%d_spin
!    scdadd%b_kin=s1%b_kin
   
    scdadd%e_ij=s1%e_ij
    scdadd%x0(1:6)=s2%x
 
    call kill(d)

  END FUNCTION scdadd 


  FUNCTION daddsc( S1,S2  )
    implicit none
    TYPE (probe_8) daddsc
    TYPE (c_damap), INTENT (IN) :: S1
    type(probe) , INTENT (IN) :: S2
    integer localmaster,i,j,nd2h,dc
    type(taylor) d
    integer rf,nd2t,nd2

    call get_rf(rf,nd2t,nd2)

    call alloc(d)

    !   call ass(daddsc)
    daddsc%u=my_false
    daddsc%E_ij=0.0_dp
    daddsc%nac=s2%nac

    daddsc%nac=s2%nac
    daddsc%use_q=s2%use_q
    if(doing_ac_modulation_in_ptc) then
       dc=2*rf
    else
       dc=0
    endif

    if(c_%ndpt==0) then    ! 1                    
      nd2h=nd2t
       do i=1,nd2h                 !   from 1-4 or 1-6 (if ndpt=0)
          localmaster=master
          call ass(daddsc%x(i))
          !       scdadd%x(i)=s1%m%v(i)+s2%x(i)

          d=s1%V(i)+s2%x(i)-(s1%V(i).sub.'0')
          daddsc%x(i)=d
          master=localmaster
       enddo
  
       do i=nd2h+1,6
          localmaster=master
          call ass(daddsc%x(i))
!  !       if((c_%npara==5+dc).AND.I==5) then   ! npr
          if(((c_%npara==5+dc).AND.I==5+ndpt_bmad).or.((c_%npara==3+dc).AND.I==5+ndpt_bmad)) then   ! npr
             daddsc%x(i)=s2%x(i)+(1.0_dp.mono.c_%npara)
          else
             daddsc%x(i)=s2%x(i)
          endif
          master=localmaster
       enddo

    else        ! 1
      nd2h=nd2t+2
       do i=1,nd2h                 !   from 1-4 or 1-6 (if ndpt=0)
          localmaster=master
          call ass(daddsc%x(i))
          !       scdadd%x(i)=s1%m%v(i)+s2%x(i)

          d=s1%V(i)+s2%x(i)-(s1%V(i).sub.'0')
          daddsc%x(i)=d
          master=localmaster
       enddo
    endif       ! 1


do i=1,daddsc%nac
    localmaster=master
    call ass(daddsc%AC(i)%x(1))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
     j=2*daddsc%nac-(2*i-1)
     d=s1%V(C_%ND2-j) +addclock*s2%AC(i)%x(1)
    daddsc%ac(i)%x(1)=d
    master=localmaster
    localmaster=master
     j=2*daddsc%nac-(2*i)
    call ass(daddsc%AC(i)%x(2))
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    d=s1%V(C_%ND2-j) +addclock*s2%AC(i)%x(2)
    daddsc%ac(i)%x(2)=d
    master=localmaster
    localmaster=master
    call ass(daddsc%AC(i)%om)
!    call ass(daddsc%AC%t)
    !       scdadd%x(i)=s1%m%v(i)+s2%x(i)
    daddsc%AC(i)%om=s2%AC(i)%om
    daddsc%AC(i)%t=s2%AC(i)%t
!    daddsc%AC(i)%f=s2%AC(i)%f
!    daddsc%AC(i)%phase=s2%AC(i)%phase
    master=localmaster
    !    endif
enddo
    DO I=1,3
       !          call ass(scdadd%s%x(i))
       DO J=1,3
          localmaster=master
          call ass(daddsc%s(J)%x(i))
          d=S1%S%s(I,J)
          daddsc%s(J)%x(i)=d
          master=localmaster
       ENDDO

    ENDDO

       !          call ass(scdadd%s%x(i))
       DO J=0,3
          localmaster=master
          call ass(daddsc%q%x(J))
          d=S1%q%x(j)
          daddsc%q%x(J)=d
          master=localmaster
       ENDDO

!    daddsc%damps=s1%damps
!    daddsc%d_spin=s1%d_spin
!    daddsc%b_kin=s1%b_kin
    daddsc%e_ij=s1%e_ij
    daddsc%x0(1:6)=s2%x
    call kill(d)

  END FUNCTION daddsc 

 

  FUNCTION daddsco( S1, S2 )
    implicit none
    TYPE (real_8) daddsco(ndd)
    TYPE (c_damap), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2(ndd)
    type(taylor) t
    integer localmaster,nd2h,i
     integer rf,nd2t,nd2

    call get_rf(rf,nd2t,nd2)

     nd2h=nd2-2*rf 
 
 !   do i=1,nd2
 !      localmaster=master
 !      call ass(daddsco(i))
 !      daddsco(i)=s1%v(i)+s2(i)-(s1%v(i).sub.'0')
 !      master=localmaster
 !   enddo
 !   do i=nd2+1,ndd
 !      localmaster=master
 !      call ass(daddsco(i))
 !      if(nd2<=4.and.(c_%npara==3.or.c_%npara==5.or.c_%npara==8).and.i==5+ndpt_bmad) then
 !         If(ndpt_bmad==0) then
 !          if(nd2==4) daddsco(i)=s2(i)+(1.0_dp.mono.'00001')
 !          if(nd2==2) daddsco(i)=s2(i)+(1.0_dp.mono.'001')
 !         endif
 !      else
 !         daddsco(i)=s2(i)
 !      endif
 !      master=localmaster
 !   enddo


    call alloc(t)

    do i=1,nd2h
       localmaster=master
       call ass(daddsco(i))
       t= s1%v(i)+s2(i)-(s1%v(i).sub.'0')
       daddsco(i)=t
       master=localmaster
    enddo
    do i=nd2h+1,ndd
       localmaster=master
       call ass(daddsco(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).and.i==5+ndpt_bmad) then
          daddsco(i)=s2(i)+(1.0_dp.mono.'00001')
       elseif(nd2==2.and.(c_%npara==3.or.c_%npara==6).and.i==5+ndpt_bmad) then
          daddsco(i)=s2(i)+(1.0_dp.mono.'001')
       else
          daddsco(i)=s2(i)
       endif
       master=localmaster
    enddo


    call kill(t)

  END FUNCTION daddsco

  FUNCTION scdaddo( S2,S1  )
    implicit none
    TYPE (real_8) scdaddo(ndd)
    TYPE (c_damap), INTENT (IN) :: S1
    real(dp) , INTENT (IN) :: S2(ndd)
    integer localmaster  ,nd2h,i
    type(taylor) t
    integer rf,nd2t,nd2

    call get_rf(rf,nd2t,nd2)

      nd2h=nd2-2*rf 

    call alloc(t)
    do i=1,nd2h
       localmaster=master
       call ass(scdaddo(i))
       t=s1%v(i)+s2(i)-(s1%v(i).sub.'0')
       scdaddo(i)=t
       master=localmaster
    enddo
    do i=nd2h+1,ndd
       localmaster=master
       call ass(scdaddo(i))
       if(nd2==4.and.(c_%npara==5.or.c_%npara==8).and.i==5+ndpt_bmad) then
          scdaddo(i)=s2(i)+(1.0_dp.mono.'00001')
       elseif(nd2==2.and.(c_%npara==3.or.c_%npara==6).and.i==5+ndpt_bmad) then
          scdaddo(i)=s2(i)+(1.0_dp.mono.'001')
       else
          scdaddo(i)=s2(i)
       endif
       master=localmaster
    enddo

    call kill(t)

  END FUNCTION scdaddo


  SUBROUTINE  equal_real8_cmap(S2,S1)
!*
    implicit none
    type (real_8),INTENT(inOUT)::S2(:)
    type (c_damap),INTENT(IN)::S1
    type(taylor) ct
    integer i
    integer rf,nd2t,nd2

    call get_rf(rf,nd2t,nd2)

    call c_check_snake

    call alloc(ct)

    do i=1,nd2
     ct=s1%v(i)

      s2(i)=ct 
    enddo
 
    call kill(ct)

 end SUBROUTINE  equal_real8_cmap

  SUBROUTINE  equal_cmap_real8(S1,S2)
!*
    implicit none
    type (real_8),INTENT(in)::S2(:)
    type (c_damap),INTENT(inOUT)::S1
    type(taylor) ct
    integer i
    integer rf,nd2t,nd2

    call get_rf(rf,nd2t,nd2)

    call check_snake

    call alloc(ct)

    do i=1,nd2
     ct=s2(i)

      s1%v(i)=ct 

    enddo
 
    call kill(ct)

 end SUBROUTINE  equal_cmap_real8


  subroutine EQUAL_c_map_RAY8(DS,R)
!*
    implicit none
    TYPE(probe_8), INTENT(IN) :: R
    TYPE(c_damap), INTENT(INOUT) :: DS
    real(dp) m(6,6)
    type(taylor) t

    INTEGER I,J,nd2t1
    logical(lp) rad_in
    integer rf,nd2t,nd2

    call get_rf(rf,nd2t,nd2)

    call alloc(t)

     nd2t1=c_%nd2-2*rf 
 !   nd2t1=C_%ND2
 !   if(doing_ac_modulation_in_ptc) then
 !      nd2t1=C_%ND2-2
 !   endif



    DO I=1,nd2t1
       t=R%X(I)
       DS%V(I)=t
    ENDDO

 !   DO I=nd2t1+1,C_%ND2
 !      t=R%ac%x(i-nd2t1)
 !      DS%V(I)=t
 !   ENDDO

    j=1
    DO I=nd2t1+1,C_%ND2,2
       t=R%ac(j)%x(1)
       DS%V(I)=t
       t=R%ac(j)%x(2)
       DS%V(I+1)=t
       j=j+1
    ENDDO

! quaternion
    if(use_quaternion)   THEN
    DO I=0,3

          t=r%q%x(i)
          DS%q%x(i)=t

    ENDDO
else
    DO I=1,3
       DO J=1,3
          t=R%S(J)%X(I)
          DS%S%s(I,J)=t
       ENDDO
    ENDDO
endif

    call check_rad(r%e_ij,rad_in)
    ds%e_ij=0.0_dp
    if(rad_in) then
       m=ds
       ds%e_ij=matmul(matmul(m,r%e_ij),transpose(m))
    endif

!ds%damps=r%damps

!ds%d_spin=r%d_spin
!ds%b_kin=r%b_kin

DS%x0=0
ds%x0(1:6)=r%x0
!if(use_quaternion) then
! call makeso3(DS%q,ds%sm)
!endif
    call kill(t)

  END subroutine EQUAL_c_map_RAY8

  subroutine EQUAL_RAY8_c_map(R,DS)
!*
    implicit none
    TYPE(probe_8), INTENT(INOUT) :: R
    TYPE(c_damap), INTENT(IN) :: DS
    real(dp) m(6,6)
    logical(lp) rad_in
    INTEGER I,J,nd2t1
    type(taylor) t
    type(c_damap) mm
    integer rf,nd2t,nd2

    call get_rf(rf,nd2t,nd2)
    call alloc(t)

    nd2t1=c_%nd2-2*rf

    DO I=1,nd2t1
       t=DS%V(I)
       R%X(I)= t  
    ENDDO
 
    
 !     do i=nd2t1+1,c_%nd2
 !      t=DS%V(I)
 !      r%ac%x(i-nd2t1) =  t
 !     enddo

    j=1
    DO I=nd2t1+1,C_%ND2,2
       t=DS%V(I)
       R%ac(j)%x(1)=t
       t=DS%V(I+1) 
       R%ac(j)%x(2)=t
       j=j+1
    ENDDO

! quaternion
if(use_quaternion)   THEN
    DO I=0,3
          t=DS%q%x(i)
          r%q%x(i)=t
    ENDDO
else
    DO J=1,3
       DO I=1,3
          t=DS%S%s(I,J) 
          R%S(J)%X(I)=t     
       ENDDO
    ENDDO
endif
    call c_check_rad(ds%e_ij,rad_in)
        r%e_ij=0.0_dp
    if(rad_in) then
     call alloc(mm)
        mm=ds
        m=mm**(-1)
         call kill(mm)
       r%e_ij=matmul(matmul(m,ds%e_ij),transpose(m))
    endif
 
   call kill(t)

  END subroutine EQUAL_RAY8_c_map

 
  ! LD: 03.04.2019
  SUBROUTINE PRTP1(S, X)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN):: S
    TYPE(REAL_8), INTENT(IN):: X

    ! cancel all PRTP
    if (MAPDUMP .eq. 0) return

    if (X%KIND /= 1) then
      ! @@ + elem + func + 7 columns
      WRITE(*, '(a,a15,a,a15,7E25.16)') '@@ ', ELEM_NAME, ' ', S, X.sub.'000000'&
                                , X.sub.'100000', X.sub.'010000', X.sub.'001000'&
                                , X.sub.'000100',-X.sub.'000001', X.sub.'000010'
    else
      ! @@ + elem + func + 1 columns
      WRITE(*, '(a,a15,a,a15,1E25.16)') '@@ ', ELEM_NAME, ' ', S, X%R
    endif
  END SUBROUTINE PRTP1

  ! LD: 22.03.2019
  SUBROUTINE PRTP(S, X)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN):: S
    TYPE(REAL_8), OPTIONAL, INTENT(IN):: X(6)

    ! cancel all PRTP
    if (MAPDUMP .eq. 0) return

    ! special case: display only string without X
    if (.not. PRESENT(X)) then
      WRITE(*, '(a,a)') '@@ ', S
      return
    endif

    ! @@ + elem + func + 6 columns
    if (MAPDUMP .eq. 1) then
      WRITE(*, '(a,a15,a,a15,6E25.16)') '@@ ', ELEM_NAME, ' ', S &
        , X(1).sub.'000000', X(2).sub.'000000', X(3).sub.'000000', X(4).sub.'000000',-X(6).sub.'000000', X(5).sub.'000000'
      return
    endif

    ! @@ + elem + func + 42 columns
    WRITE(*, '(a,a15,a,a15,42E25.16)') '@@ ', ELEM_NAME, ' ', S &
      , X(1).sub.'000000', X(2).sub.'000000', X(3).sub.'000000', X(4).sub.'000000',-X(6).sub.'000000', X(5).sub.'000000'&
      , X(1).sub.'100000', X(1).sub.'010000', X(1).sub.'001000', X(1).sub.'000100',-X(1).sub.'000001', X(1).sub.'000010'&
      , X(2).sub.'100000', X(2).sub.'010000', X(2).sub.'001000', X(2).sub.'000100',-X(2).sub.'000001', X(2).sub.'000010'&
      , X(3).sub.'100000', X(3).sub.'010000', X(3).sub.'001000', X(3).sub.'000100',-X(3).sub.'000001', X(3).sub.'000010'&
      , X(4).sub.'100000', X(4).sub.'010000', X(4).sub.'001000', X(4).sub.'000100',-X(4).sub.'000001', X(4).sub.'000010'&
      ,-X(6).sub.'100000',-X(6).sub.'010000',-X(6).sub.'001000',-X(6).sub.'000100', X(6).sub.'000001',-X(6).sub.'000010'&
      , X(5).sub.'100000', X(5).sub.'010000', X(5).sub.'001000', X(5).sub.'000100',-X(5).sub.'000001', X(5).sub.'000010'
  END SUBROUTINE PRTP

  SUBROUTINE ANALYSE_APERTURE_FLAG(I,R)
    IMPLICIT NONE
    INTEGER I,B,K
    INTEGER :: R(:)

    K=I
    B=1
    r=-1
    DO WHILE (K>0.AND.B<=SIZE(R))
       R(B)=MOD(K,2)
       IF(MOD(K,2)==1) THEN
          K=(K-1)/2
       ELSE
          K=K/2
       ENDIF
       B=B+1
    ENDDO

  END   SUBROUTINE ANALYSE_APERTURE_FLAG




  REAL(DP) FUNCTION  SINEHX_X(X) ! REPLACES SINH(X)/X
    IMPLICIT NONE
    REAL(DP),INTENT(IN)::X
    IF(.NOT.c_%CHECK_STABLE) then
       sinehx_x=1.0_dp
       return
    endif

    IF((ABS(X)>hyperbolic_aperture).AND.ROOT_CHECK) THEN
       SINEHX_X=0.0_dp
       CHECK_STABLE=.FALSE.
       messagelost="Sa_extend_poly.f90 SINEXHX_X : argument out of range" !CERN
    ELSEIF(ABS(X)<=hyperbolic_aperture) THEN
       sinehx_x = sinhx_x(x)
    ELSE      !  IF X IS NOT A NUMBER
       sinehx_x=1.0_dp
       CHECK_STABLE=.FALSE.
       messagelost="Sa_extend_poly.f90 SINEXHX_X : should never happen" !CERN
    ENDIF

  END FUNCTION SINEHX_X


  ! Some polymorphism


  ! End of Some polymorphism

end module S_extend_poly



