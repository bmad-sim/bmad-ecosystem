module elliptic_integral_mod

use, intrinsic :: iso_fortran_env

implicit none


integer, parameter, private :: sp = REAL32
integer, parameter, private :: dp = REAL64

contains


!----------------------------------------------------------------------
! subroutine ellipinc(phi, m, ellipkinc, ellipeinc)
!
! Calculates the first and second incomplete elliptic integrals,
! using methods from T. Fukushima, (2011, 2018)
!
! Uses classical transformations to handle negative m.
! This package needs a function for the third kind to use the new 2018 transformations.
!
! Input
! -----
! phi : real
!     amplitude of the elliptic integral
! m   : real
!        parameter of the elliptic integral. Requires m < 1
!
! Output
! ------
! ellipkinc : real
!            incomplete elliptic integral of the first kind K(phi | m)
!            also called F(phi | m)
!
! ellipeinc :real 
!            incomplete elliptic integral of the second kind F(phi | m),
!
! References
! ----------
! T. Fukushima, (2011), J. Comp. Appl. Math., 235, 4140-4148
!        "Precise and Fast Computation of General Incomplete Elliptic Integral
!         of Second Kind by Half and Double Argument Transformations"
! https://doi.org/10.1016/j.cam.2011.03.004
!
! T. Fukushima, New variable transformation of elliptic integrals (2018)
!
! https://arxiv.org/abs/math/9409227
! Scipy modifications for negative m use different methods:
! https://github.com/scipy/scipy/blob/3f017716532001349ccddb05498408d711bacdde/scipy/special/cephes/ellik.c#L160
!
!----------------------------------------------------------------------
elemental subroutine ellipinc(phi, m, ellipkinc, ellipeinc)

real(dp), intent(in)  :: phi, m
real(dp), intent(out) :: ellipkinc, ellipeinc
real(dp) :: m1, phi1, d, elb,eld


! Negative parameters need a special transformation, 
! 'classically called the imaginary modulus transformation' - Fukushima 2018
! Use his notation, except with original m = -m. 
! Note there is a typo in Eq. 8. All phi on the r.h.s should be phi1
! 
! Alternatively, Appendix B. from Fukushima could be used. 
if (m<0) then
  m1 = m/(m-1) 
  d = sqrt(1-m) ! denominator in Eq. 7
  phi1 = asin(d*sin(phi)/sqrt(1-m*sin(phi)**2)) ! Eq. 10

  ! Actual calc takes the complementary parameter mc = 1-m
  call gelbd(phi1,(1-m1),elb,eld)
  
  ellipkinc = (elb + eld)/d  
  ellipeinc = d*( (elb + (1-m1)*eld) - m1*sin(phi1)*cos(phi1)/sqrt(1-m1*sin(phi1)**2))
    
else

  ! Calculate the associate incomplete integrals B(phi|m), D(phi|m)
  call gelbd(phi,(1-m),elb,eld)
  
  ! The usual incomplete elliptic integrals are:
  ! K(phi|m) = B(phi|m) + D(phi|m)
  ! 
  ellipkinc = elb + eld
  
  ! E(phi|m) =  = B(phi|m) + (1-m)*D(phi|m)
  ellipeinc = elb + (1-m)*eld

endif

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! T. Fukushima code below
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine test_xgelbd()
!
! test driver for gelbd
!
implicit none

real(dp) ::  dmc,mc,m,dphi,phi,b,d
real(sp) ::  rmc,rphi,rb,rd,rdb,rdd
integer jend,iend,j,i
!
jend=10
iend=2
!
dmc=1.d0/dble(jend)
dphi=1.d0/dble(iend)
write(*,'(1x,2a10,2a25,2a10)') "m","phi","B","D","rB-B","rD-D"
!
do j=1,jend
    mc=dble(j)*dmc
  m=1.d0-mc
    rmc=real(mc, sp)
    do i=0,iend*10
      phi=dphi*dble(i)
        rphi=real(phi, sp)
      call gelbd(phi,mc,b,d)
      call rgelbd(rphi,rmc,rb,rd)
        rdb=real((rb-b)/(1.d-16+abs(b)), sp)
        rdd=real((rd-d)/(1.d-16+abs(d)), sp)
        write(*,'(1x,0p2f10.5,1p2e25.15,1p2e10.2)') m,phi,b,d,rdb,rdd
  enddo
    write(*,'(1x)')
enddo
end subroutine test_xgelbd




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elemental subroutine gelbd(phi,mc,elb,eld)
!
! Double precision associate incomplete elliptic integrals of the second kind
!
!   For general amplitude: -infty < phi < infty
!
!     Reference: T. Fukushima, (2011), J. Comp. Appl. Math., 235, 4140-4148
!        "Precise and Fast Computation of General Incomplete Elliptic Integral
!         of Second Kind by Half and Double Argument Transformations"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Used subprograms: elbd, celbd
!
!     Inputs: phi = amplitude, mc = complementary parameter, 0 <= m < 1
!
!     Outputs: elb = B(phi|m), eld = D(phi|m)
!
real(dp), intent(in) ::  phi,mc
real(dp), intent(out) ::  elb,eld
real(dp) ::  m,phix,celb,celd,phic
real(dp) ::  PIHALF,PI,TWOPI,PIH3
integer n
parameter (PI=3.141592653589793238462d0)
parameter (PIHALF=PI*0.5d0,TWOPI=PI*2.d0,PIH3=PI*1.5d0)
!
if(mc.le.0.d0) then
!  write(*,*) "(gelbd) too small parameter: mc=",mc
    return
endif
m=1.d0-mc
phix=abs(phi)
call celbd(mc,celb,celd)
n=int(phix/TWOPI)
phix=phix-dble(n)*TWOPI
if(phix.lt.PIHALF) then
    phic=PIHALF-phix; call elbd(phix,phic,mc,elb,eld)
elseif(phix.lt.PI) then
    phix=PI-phix
    phic=PIHALF-phix; call elbd(phix,phic,mc,elb,eld)
    elb=2.d0*celb-elb
    eld=2.d0*celd-eld
elseif(phix.lt.PIH3) then
    phix=phix-PI
    phic=PIHALF-phix; call elbd(phix,phic,mc,elb,eld)
    elb=elb+2.d0*celb
    eld=eld+2.d0*celd
else
    phix=TWOPI-phix
    phic=PIHALF-phix; call elbd(phix,phic,mc,elb,eld)
    elb=4.d0*celb-elb
    eld=4.d0*celd-eld
endif
if(n.gt.0) then
    elb=elb+dble(4*n)*celb
    eld=eld+dble(4*n)*celd
endif
if(phi.lt.0.d0) then
    elb=-elb
    eld=-eld
endif
return
end subroutine gelbd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rgelbd(phi,mc,elb,eld)
!
! Single precision associate incomplete elliptic integrals of the second kind
!
!   For general amplitude: -infty < phi < infty
!
!     Reference: T. Fukushima, (2011), J. Comp. Appl. Math., 235, 4140-4148
!        "Precise and Fast Computation of General Incomplete Elliptic Integral
!         of Second Kind by Half and Double Argument Transformations"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Used subprograms: relbd, rcelbd
!
!     Inputs: phi = amplitude, mc = complementary parameter, 0 <= m < 1
!
!     Outputs: elb = B(phi|m), eld = D(phi|m)
!
real(sp) ::  phi,mc,elb,eld
real(sp) ::  m,phix,celb,celd,phic
real(sp) ::  PIHALF,PI,TWOPI,PIH3
integer n
parameter (PI=3.14159265)
parameter (PIHALF=PI*0.5,TWOPI=PI*2.0,PIH3=PI*1.5)
!
if(mc.le.0.0) then
  write(*,*) "(rgelbd) too small parameter: mc=",mc
    return
endif
m=1.0-mc
phix=abs(phi)
call rcelbd(mc,celb,celd)
n=int(phix/TWOPI)
phix=phix-real(n)*TWOPI
if(phix.lt.PIHALF) then
    phic=PIHALF-phix; call relbd(phix,phic,mc,elb,eld)
elseif(phix.lt.PI) then
    phix=PI-phix
    phic=PIHALF-phix; call relbd(phix,phic,mc,elb,eld)
    elb=2.0*celb-elb
    eld=2.0*celd-eld
elseif(phix.lt.PIH3) then
    phix=phix-PI
    phic=PIHALF-phix; call relbd(phix,phic,mc,elb,eld)
    elb=elb+2.0*celb
    eld=eld+2.0*celd
else
    phix=TWOPI-phix
    phic=PIHALF-phix; call relbd(phix,phic,mc,elb,eld)
    elb=4.0*celb-elb
    eld=4.0*celd-eld
endif
if(n.gt.0) then
    elb=elb+real(4*n)*celb
    eld=eld+real(4*n)*celd
endif
if(phi.lt.0.0) then
    elb=-elb
    eld=-eld
endif
return
end subroutine rgelbd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elemental subroutine elbd(phi,phic,mc,b,d)
!
! Double precision associate incomplete elliptic integrals of the second kind
!
!   For standard amplitude: 0 <= phi < Pi/2
!
!     Reference: T. Fukushima, (2011), J. Comp. Appl. Math., 235, 4140-4148
!        "Precise and Fast Computation of General Incomplete Elliptic Integral
!         of Second Kind by Half and Double Argument Transformations"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Used subprograms: celbd,elsbd,elcbd
!
!     Inputs: phi = amplitude, mc = complementary parameter, 0 <= m < 1
!
!     Outputs: elb = B(phi|m), eld = D(phi|m)
!
real(dp), intent(in) :: phi,phic,mc
real(dp), intent(out) :: b,d
real(dp) ::  m,c,x,d2,z,bc,dc,sz,v,t2

if(phi.lt.1.25d0) then
    call elsbd(sin(phi),mc,b,d)
else
    m=1.d0-mc
    c=sin(phic)
    x=c*c
    d2=mc+m*x
    if(x.lt.0.9d0*d2) then
        z=c/sqrt(d2)
        call elsbd(z,mc,b,d)
        call celbd(mc,bc,dc)
    sz=z*sqrt(1.d0-x)
    b=bc-(b-sz)
    d=dc-(d+sz)
  else
    v=mc*(1.d0-x)
    if(v.lt.x*d2) then
            call elcbd(c,mc,b,d)
        else
            t2=(1.d0-x)/d2
            call elcbd(sqrt(mc*t2),mc,b,d)
            call celbd(mc,bc,dc)
      sz=c*sqrt(t2)
      b=bc-(b-sz)
      d=dc-(d+sz)
    endif
  endif
endif
return
end subroutine elbd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elemental subroutine elsbd(s0,mc,b,d)

real(dp), intent(in) ::  s0,mc
real(dp), intent(out) ::  b,d
real(dp) ::  m,del,s,y,sy
real(dp) ::  yy(11),ss(11)
integer i, j

! write(*,*) "(elsbd) s0,m=",s0,m

m=1.d0-mc
del=0.04094d0-0.00652d0*m ! F9  Optimum
s=s0
y=s*s
if(y.lt.del) then
  call serbd(y,m,b,d)
  b=s*b
  d=s*y*d
  return
endif
ss(1)=s
do j=1,10
    y=y/((1.d0+sqrt(1.d0-y))*(1.d0+sqrt(1.d0-m*y)))
  yy(j+1)=y
  ss(j+1)=sqrt(y)
  if(y.lt.del) then
    goto 1
  endif
enddo
!write(*,*) "(elsbd) too many iterations: s0,m=",s0,m
1 continue
! write(*,*) 'j=',j
call serbd(y,m,b,d)
b=ss(j+1)*b
d=ss(j+1)*y*d
do i=1,j
  sy=ss(j-i+1)*yy(j-i+2)
  b=b+(b-sy)
  d=d+(d+sy)
enddo
return
end subroutine elsbd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elemental subroutine elcbd(c0,mc,b,dx)

real(dp), intent(in)  ::  c0,mc
real(dp), intent(out) ::  b,dx
real(dp) ::  c,x,y,s,m,d,sy
real(dp) ::  yy(11),ss(11)
integer j,i

! write(*,*) "(elcbd) c0,mc=",c0,mc

c=c0
x=c*c
y=1.d0-x
s=sqrt(y)
if(x.gt.0.1d0) then
! write(*,*) "(elcbd) elsbd"
  call elsbd(s,mc,b,dx)
    return
endif
m=1.d0-mc
ss(1)=s
! write(*,*) "(elcbd) j,y,s=",1,y,ss(1)
do j=1,10
    d=sqrt(mc+m*x)
  x=(c+d)/(1.d0+d)
  y=1.d0-x
  yy(j+1)=y
  ss(j+1)=sqrt(y)
! write(*,*) "(elcbd) j,y,s=",j+1,yy(j+1),ss(j+1)
    if(x.gt.0.1d0) then
        goto 1
  endif
  c=sqrt(x)
enddo
!write(*,*) "(elcbd) too many iterations: c0,mc=",c0,mc
1 continue
s=ss(j+1)
call elsbd(s,mc,b,dx)
do i=1,j
  sy=ss(j-i+1)*yy(j-i+2)
  b=b+(b-sy)
  dx=dx+(dx+sy)
enddo
return
end subroutine elcbd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elemental subroutine celbd(mc,elb,eld)
!
! Double precision associate complete elliptic integral of the second kind
!
!     Reference: T. Fukushima, (2011), Math. Comp., 80, 1725-1743
!        "Precise and Fast Computation of General Complete Elliptic Integral
!         of Second Kind"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Inputs: mc = complementary parameter, 0 < mc <= 1
!
!     Outputs: elb = B(m), eld = D(m)
!
real(dp), intent(in)  ::  mc
real(dp), intent(out) ::  elb,eld
real(dp) ::  m,nome,dkkc,dddc,mx,kkc,logq2,elk,dele,elk1,delb

real(dp) ::  PIQ,PIHALF
parameter (PIQ=0.78539816339744830961566084581988d0)
parameter (PIHALF=1.5707963267948966192313216916398d0)

real(dp) ::  mcold,elbold,eldold
! save mcold,elbold,eldold

real(dp) ::  Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16
parameter (Q1=1.d0/16.d0,Q2=1.d0/32.d0,Q3=21.d0/1024.d0)
parameter (Q4=31.d0/2048.d0,Q5=6257.d0/524288.d0)
parameter (Q6=10293.d0/1048576.d0,Q7=279025.d0/33554432.d0)
parameter (Q8=483127.d0/67108864.d0)
parameter (Q9=435506703.d0/68719476736.d0)
parameter (Q10=776957575.d0/137438953472.d0)
parameter (Q11=22417045555.d0/4398046511104.d0)
parameter (Q12=40784671953.d0/8796093022208.d0)
parameter (Q13=9569130097211.d0/2251799813685248.d0)
parameter (Q14=17652604545791.d0/4503599627370496.d0)
parameter (Q15=523910972020563.d0/144115188075855872.d0)
parameter (Q16=976501268709949.d0/288230376151711744.d0)

real(dp) ::  K1,K2,K3,K4,K5,K6,K7
parameter (K1=1.d0/4.d0)
parameter (K2=9.d0/64.d0)
parameter (K3=25.d0/256.d0)
parameter (K4=1225.d0/16384.d0)
parameter (K5=3969.d0/65536.d0)
parameter (K6=53361.d0/1048576.d0)
parameter (K7=184041.d0/4194304.d0)

real(dp) ::  B1,B2,B3,B4,B5,B6,B7,B8
parameter (B1=1.d0/2.d0)
parameter (B2=1.d0/16.d0)
parameter (B3=3.d0/128.d0)
parameter (B4=25.d0/2048.d0)
parameter (B5=245.d0/32768.d0)
parameter (B6=1323.d0/262144.d0)
parameter (B7=7623.d0/2097152.d0)
parameter (B8=184041.d0/67108864.d0)

real(dp) ::  D1,D2,D3,D4,D5,D6,D7,D8
parameter (D1=1.d0/2.d0)
parameter (D2=3.d0/16.d0)
parameter (D3=15.d0/128.d0)
parameter (D4=175.d0/2048.d0)
parameter (D5=2205.d0/32768.d0)
parameter (D6=14553.d0/262144.d0)
parameter (D7=99099.d0/2097152.d0)
parameter (D8=2760615.d0/67108864.d0)

! logical first/.TRUE./

!if(first) then
!    first=.FALSE.
  mcold=1.d0
  elbold=PIQ
  eldold=PIQ
!endif
m=1.d0-mc
if(abs(mc-mcold).lt.1.11d-16*mc) then
  elb=elbold
  eld=eldold
elseif(m.lt.1.11d-16) then
    elb=PIQ
  eld=PIQ
elseif(mc.lt.1.11d-16) then
    elb=1.d0
  eld=0.3862943611198906188344642429164d0-0.5d0*log(mc)
elseif(mc.lt.0.1d0) then
    nome=mc*(Q1+mc*(Q2+mc*(Q3+mc*(Q4+mc*(Q5+mc*(Q6 &
        +mc*(Q7+mc*(Q8+mc*(Q9+mc*(Q10+mc*(Q11+mc*(Q12 &
        +mc*(Q13+mc*(Q14+mc*(Q15+mc*Q16))))))))))))))) 
    if(mc.lt.0.01d0) then
        dkkc=mc*(K1+mc*(K2+mc*(K3+mc*(K4+mc*(K5+mc*(K6+mc*K7))))))
        dddc=mc*(D1+mc*(D2+mc*(D3+mc*(D4+mc*(D5+mc*(D6+mc*D7))))))
    else
        mx=mc-0.05d0

! (K'-1)/(pi/2)

        dkkc=    0.01286425658832983978282698630501405107893d0 &
            +mx*(0.26483429894479586582278131697637750604652d0 &
            +mx*(0.15647573786069663900214275050014481397750d0 &
            +mx*(0.11426146079748350067910196981167739749361d0 &
            +mx*(0.09202724415743445309239690377424239940545d0 &
            +mx*(0.07843218831801764082998285878311322932444d0 &
            +mx*(0.06935260142642158347117402021639363379689d0 &
            +mx*(0.06293203529021269706312943517695310879457d0 &
            +mx*(0.05821227592779397036582491084172892108196d0 &
            +mx*(0.05464909112091564816652510649708377642504d0 &
            +mx*(0.05191068843704411873477650167894906357568d0 &
            +mx*(0.04978344771840508342564702588639140680363d0 &
            +mx*(0.04812375496807025605361215168677991360500d0 &
            ))))))))))))

! (K'-E')/(pi/2)

        dddc=    0.02548395442966088473597712420249483947953d0 &
            +mx*(0.51967384324140471318255255900132590084179d0 &
            +mx*(0.20644951110163173131719312525729037023377d0 &
            +mx*(0.13610952125712137420240739057403788152260d0 &
            +mx*(0.10458014040566978574883406877392984277718d0 &
            +mx*(0.08674612915759188982465635633597382093113d0 &
            +mx*(0.07536380269617058326770965489534014190391d0 &
            +mx*(0.06754544594618781950496091910264174396541d0 &
            +mx*(0.06190939688096410201497509102047998554900d0 &
            +mx*(0.05771071515451786553160533778648705873199d0 &
            +mx*(0.05451217098672207169493767625617704078257d0 &
            +mx*(0.05204028407582600387265992107877094920787d0 &
            +mx*(0.05011532514520838441892567405879742720039d0 &
            ))))))))))))
    endif
  kkc=1.d0+dkkc
  logq2=-0.5d0*log(nome)
  elk=kkc*logq2
  dele=-dkkc/kkc+logq2*dddc
  elk1=elk-1.d0
  delb=(dele-mc*elk1)/m
  elb=1.d0+delb
  eld=elk1-delb
elseif(m.le.0.01d0) then
  elb=PIHALF*(B1+m*(B2+m*(B3+m*(B4+m*(B5+m*(B6+m*(B7+m*B8)))))))
  eld=PIHALF*(D1+m*(D2+m*(D3+m*(D4+m*(D5+m*(D6+m*(D7+m*D8)))))))
elseif(m.le.0.1d0) then
  mx=0.95d0-mc
  elb=     0.790401413584395132310045630540381158921005d0 &
      +mx*(0.102006266220019154892513446364386528537788d0 &
      +mx*(0.039878395558551460860377468871167215878458d0 &
      +mx*(0.021737136375982167333478696987134316809322d0 &
      +mx*(0.013960979767622057852185340153691548520857d0 &
      +mx*(0.009892518822669142478846083436285145400444d0 &
      +mx*(0.007484612400663335676130416571517444936951d0 &
      +mx*(0.005934625664295473695080715589652011420808d0 &
      +mx*(0.004874249053581664096949448689997843978535d0 &
      +mx*(0.004114606930310886136960940893002069423559d0 &
      +mx*(0.003550452989196176932747744728766021440856d0 &
      +mx*(0.003119229959988474753291950759202798352266d0 &
      )))))))))))
  eld=     0.800602040206397047799296975176819811774784d0 &
      +mx*(0.313994477771767756849615832867393028789057d0 &
      +mx*(0.205913118705551954501930953451976374435088d0 &
      +mx*(0.157744346538923994475225014971416837073598d0 &
      +mx*(0.130595077319933091909091103101366509387938d0 &
      +mx*(0.113308474489758568672985167742047066367053d0 &
      +mx*(0.101454199173630195376251916342483192174927d0 &
      +mx*(0.0929187842072974367037702927967784464949434d0 &
      +mx*(0.0865653801481680871714054745336652101162894d0 &
      +mx*(0.0817279846651030135350056216958053404884715d0 &
      +mx*(0.0779906657291070378163237851392095284454654d0 &
      +mx*(0.075080426851268007156477347905308063808848d0 &
      )))))))))))
elseif(m.le.0.2d0) then
  mx=0.85d0-mc
  elb=     0.80102406445284489393880821604009991524037d0 &
      +mx*(0.11069534452963401497502459778015097487115d0 &
      +mx*(0.047348746716993717753569559936346358937777d0 &
      +mx*(0.028484367255041422845322166419447281776162d0 &
      +mx*(0.020277811444003597057721308432225505126013d0 &
      +mx*(0.015965005853099119442287313909177068173564d0 &
      +mx*(0.013441320273553634762716845175446390822633d0 &
      +mx*(0.011871565736951439501853534319081030547931d0 &
      +mx*(0.010868363672485520630005005782151743785248d0 &
      +mx*(0.010231587232710564565903812652581252337697d0 &
      +mx*(0.009849585546666211201566987057592610884309d0 &
      +mx*(0.009656606347153765129943681090056980586986d0 &
      )))))))))))
  eld=     0.834232667811735098431315595374145207701720d0 &
      +mx*(0.360495281619098275577215529302260739976126d0 &
      +mx*(0.262379664114505869328637749459234348287432d0 &
      +mx*(0.223723944518094276386520735054801578584350d0 &
      +mx*(0.206447811775681052682922746753795148394463d0 &
      +mx*(0.199809440876486856438050774316751253389944d0 &
      +mx*(0.199667451603795274869211409350873244844882d0 &
      +mx*(0.204157558868236842039815028663379643303565d0 &
      +mx*(0.212387467960572375038025392458549025660994d0 &
      +mx*(0.223948914061499360356873401571821627069173d0 &
      +mx*(0.238708097425597860161720875806632864507536d0 &
      +mx*(0.256707203545463755643710021815937785120030d0 &
      )))))))))))
elseif(m.le.0.3d0) then
  mx=0.75d0-mc
  elb=     0.81259777291992049322557009456643357559904d0 &
      +mx*(0.12110961794551011284012693733241967660542d0 &
      +mx*(0.057293376831239877456538980381277010644332d0 &
      +mx*(0.038509451602167328057004166642521093142114d0 &
      +mx*(0.030783430301775232744816612250699163538318d0 &
      +mx*(0.027290564934732526869467118496664914274956d0 &
      +mx*(0.025916369289445198731886546557337255438083d0 &
      +mx*(0.025847203343361799141092472018796130324244d0 &
      +mx*(0.026740923539348854616932735567182946385269d0 &
      +mx*(0.028464314554825704963640157657034405579849d0 &
      +mx*(0.030995446237278954096189768338119395563447d0 &
      +mx*(0.034384369179940975864103666824736551261799d0 &
      +mx*(0.038738002072493935952384233588242422046537d0 &
      ))))))))))))
  eld=     0.873152581892675549645633563232643413901757d0 &
      +mx*(0.420622230667770215976919792378536040460605d0 &
      +mx*(0.344231061559450379368201151870166692934830d0 &
      +mx*(0.331133021818721761888662390999045979071436d0 &
      +mx*(0.345277285052808411877017306497954757532251d0 &
      +mx*(0.377945322150393391759797943135325823338761d0 &
      +mx*(0.427378012464553880508348757311170776829930d0 &
      +mx*(0.494671744307822405584118022550673740404732d0 &
      +mx*(0.582685115665646200824237214098764913658889d0 &
      +mx*(0.695799207728083164790111837174250683834359d0 &
      +mx*(0.840018401472533403272555302636558338772258d0 &
      +mx*(1.023268503573606060588689738498395211300483d0 &
      +mx*(1.255859085136282496149035687741403985044122d0 &
      ))))))))))))
elseif(m.le.0.4d0) then
  mx=0.65d0-mc
  elb=     0.8253235579835158949845697805395190063745d0 &
      +mx*(0.1338621160836877898575391383950840569989d0 &
      +mx*(0.0710112935979886745743770664203746758134d0 &
      +mx*(0.0541784774173873762208472152701393154906d0 &
      +mx*(0.0494517449481029932714386586401273353617d0 &
      +mx*(0.0502221962241074764652127892365024315554d0 &
      +mx*(0.0547429131718303528104722303305931350375d0 &
      +mx*(0.0627462579270016992000788492778894700075d0 &
      +mx*(0.0746698810434768864678760362745179321956d0 &
      +mx*(0.0914808451777334717996463421986810092918d0 &
      +mx*(0.1147050921109978235104185800057554574708d0 &
      +mx*(0.1465711325814398757043492181099197917984d0 &
      +mx*(0.1902571373338462844225085057953823854177d0 &
      ))))))))))))
  eld=     0.9190270392420973478848471774160778462738d0 &
      +mx*(0.5010021592882475139767453081737767171354d0 &
      +mx*(0.4688312705664568629356644841691659415972d0 &
      +mx*(0.5177142277764000147059587510833317474467d0 &
      +mx*(0.6208433913173031070711926900889045286988d0 &
      +mx*(0.7823643937868697229213240489900179142670d0 &
      +mx*(1.0191145350761029126165253557593691585239d0 &
      +mx*(1.3593452027484960522212885423056424704073d0 &
      +mx*(1.8457173023588279422916645725184952058635d0 &
      +mx*(2.5410717031539207287662105618152273788399d0 &
      +mx*(3.5374046552080413366422791595672470037341d0 &
      +mx*(4.9692960029774259303491034652093672488707d0 &
      +mx*(7.0338228700300311264031522795337352226926d0 &
      +mx*(10.020043225034471401553194050933390974016d0 &
      )))))))))))))
elseif(m.le.0.5d0) then
  mx=0.55d0-mc
  elb=     0.8394795702706129706783934654948360410325d0 &
      +mx*(0.1499164403063963359478614453083470750543d0 &
      +mx*(0.0908319358194288345999005586556105610069d0 &
      +mx*(0.0803470334833417864262134081954987019902d0 &
      +mx*(0.0856384405004704542717663971835424473169d0 &
      +mx*(0.1019547259329903716766105911448528069506d0 &
      +mx*(0.1305748115336160150072309911623351523284d0 &
      +mx*(0.1761050763588499277679704537732929242811d0 &
      +mx*(0.2468351644029554468698889593583314853486d0 &
      +mx*(0.3564244768677188553323196975301769697977d0 &
      +mx*(0.5270025622301027434418321205779314762241d0 &
      +mx*(0.7943896342593047502260866957039427731776d0 &
      +mx*(1.2167625324297180206378753787253096783993d0 &
      ))))))))))))
  eld=     0.9744043665463696730314687662723484085813d0 &
      +mx*(0.6132468053941609101234053415051402349752d0 &
      +mx*(0.6710966695021669963502789954058993004082d0 &
      +mx*(0.8707276201850861403618528872292437242726d0 &
      +mx*(1.2295422312026907609906452348037196571302d0 &
      +mx*(1.8266059675444205694817638548699906990301d0 &
      +mx*(2.8069345309977627400322167438821024032409d0 &
      +mx*(4.4187893290840281339827573139793805587268d0 &
      +mx*(7.0832360574787653249799018590860687062869d0 &
      +mx*(11.515088120557582942290563338274745712174d0 &
      +mx*(18.931511185999274639516732819605594455165d0 &
      +mx*(31.411996938204963878089048091424028309798d0 &
      +mx*(52.520729454575828537934780076768577185134d0 &
      +mx*(88.384854735065298062125622417251073520996d0 &
      +mx*(149.56637449398047835236703116483092644714d0 &
      +mx*(254.31790843104117434615624121937495622372d0 &
      )))))))))))))))
elseif(m.le.0.6d0) then
  mx=0.45d0-mc
  elb=     0.8554696151564199914087224774321783838373d0 &
      +mx*(0.1708960726897395844132234165994754905373d0 &
      +mx*(0.1213352290269482260207667564010437464156d0 &
      +mx*(0.1282018835749474096272901529341076494573d0 &
      +mx*(0.1646872814515275597348427294090563472179d0 &
      +mx*(0.2374189087493817423375114793658754489958d0 &
      +mx*(0.3692081047164954516884561039890508294508d0 &
      +mx*(0.6056587338479277173311618664015401963868d0 &
      +mx*(1.0337055615578127436826717513776452476106d0 &
      +mx*(1.8189884893632678849599091011718520567105d0 &
      +mx*(3.2793776512738509375806561547016925831128d0 &
      +mx*(6.0298883807175363312261449542978750456611d0 &
      +mx*(11.269796855577941715109155203721740735793d0 &
      +mx*(21.354577850382834496786315532111529462693d0 &
      )))))))))))))
  eld=     1.04345529511513353426326823569160142342838d0 &
      +mx*(0.77962572192850485048535711388072271372632d0 &
      +mx*(1.02974236093206758187389128668777397528702d0 &
      +mx*(1.62203722341135313022433907993860147395972d0 &
      +mx*(2.78798953118534762046989770119382209443756d0 &
      +mx*(5.04838148737206914685643655935236541332892d0 &
      +mx*(9.46327761194348429539987572314952029503864d0 &
      +mx*(18.1814899494276679043749394081463811247757d0 &
      +mx*(35.5809805911791687037085198750213045708148d0 &
      +mx*(70.6339354619144501276254906239838074917358d0 &
      +mx*(141.828580083433059297030133195739832297859d0 &
      +mx*(287.448751250132166257642182637978103762677d0 &
      +mx*(587.115384649923076181773192202238389711345d0 &
      +mx*(1207.06543522548061603657141890778290399603d0 &
      +mx*(2495.58872724866422273012188618178997342537d0 &
      +mx*(5184.69242939480644062471334944523925163600d0 &
      +mx*(10817.2133369041327524988910635205356016939d0 &
      ))))))))))))))))
elseif(m.le.0.7d0) then
  mx=0.35d0-mc
  elb=     0.8739200618486431359820482173294324246058d0 &
      +mx*(0.1998140574823769459497418213885348159654d0 &
      +mx*(0.1727696158780152128147094051876565603862d0 &
      +mx*(0.2281069132842021671319791750725846795701d0 &
      +mx*(0.3704681411180712197627619157146806221767d0 &
      +mx*(0.6792712528848205545443855883980014994450d0 &
      +mx*(1.3480084966817573020596179874311042267679d0 &
      +mx*(2.8276709768538207038046918250872679902352d0 &
      +mx*(6.1794682501239140840906583219887062092430d0 &
      +mx*(13.935686010342811497608625663457407447757d0 &
      +mx*(32.218929281059722026322932181420383764028d0 &
      +mx*(76.006962959226101026399085304912635262362d0 &
      +mx*(182.32144908775406957609058046006949657416d0 &
      +mx*(443.51507644112648158679360783118806161062d0 &
      +mx*(1091.8547229028388292980623647414961662223d0 &
      +mx*(2715.7658664038195881056269799613407111521d0 &
      )))))))))))))))
  eld=     1.13367833657573316571671258513452768536080d0 &
      +mx*(1.04864317372997039116746991765351150490010d0 &
      +mx*(1.75346504119846451588826580872136305225406d0 &
      +mx*(3.52318272680338551269021618722443199230946d0 &
      +mx*(7.74947641381397458240336356601913534598302d0 &
      +mx*(17.9864500558507330560532617743406294626849d0 &
      +mx*(43.2559163462326133313977294448984936591235d0 &
      +mx*(106.681534454096017031613223924991564429656d0 &
      +mx*(268.098486573117433951562111736259672695883d0 &
      +mx*(683.624114850289804796762005964155730439745d0 &
      +mx*(1763.49708521918740723028849567007874329637d0 &
      +mx*(4592.37475383116380899419201719007475759114d0 &
      +mx*(12053.4410190488892782190764838488156555734d0 &
      +mx*(31846.6630207420816960681624497373078887317d0 &
      +mx*(84621.2213590568080177035346867495326879117d0 &
      +mx*(225956.423182907889987641304430180593010940d0 &
      +mx*(605941.517281758859958050194535269219533685d0 &
      +mx*(1.63108259953926832083633544697688841456604d6 &
      )))))))))))))))))
elseif(m.le.0.8d0) then
  mx=0.25d0-mc
  elb=     0.895902820924731621258525533131864225704d0 &
      +mx*(0.243140003766786661947749288357729051637d0 &
      +mx*(0.273081875594105531575351304277604081620d0 &
      +mx*(0.486280007533573323895498576715458103274d0 &
      +mx*(1.082747437228230914750752674136983406683d0 &
      +mx*(2.743445290986452500459431536349945437824d0 &
      +mx*(7.555817828670234627048618342026400847824d0 &
      +mx*(22.05194082493752427472777448620986154515d0 &
      +mx*(67.15640644740229407624192175802742979626d0 &
      +mx*(211.2722537881770961691291434845898538537d0 &
      +mx*(681.9037843053270682273212958093073895805d0 &
      +mx*(2246.956231592536516768812462150619631201d0 &
      +mx*(7531.483865999711792004783423815426725079d0 &
      +mx*(25608.51260130241579018675054866136922157d0 &
      +mx*(88140.74740089604971425934283371277143256d0 &
      +mx*(306564.4242098446591430938434419151070722d0 &
      +mx*(1.076036077811072193752770590363885180738d6 &
      +mx*(3.807218502573632648224286313875985190526d6 &
      +mx*(1.356638224422139551020110323739879481042d7 &
      ))))))))))))))))))
  eld=     1.26061282657491161418014946566845780315983d0 &
      +mx*(1.54866563808267658056930177790599939977154d0 &
      +mx*(3.55366941187160761540650011660758187283401d0 &
      +mx*(9.90044467610439875577300608183010716301714d0 &
      +mx*(30.3205666174524719862025105898574414438275d0 &
      +mx*(98.1802586588830891484913119780870074464833d0 &
      +mx*(329.771010434557055036273670551546757245808d0 &
      +mx*(1136.65598974289039303581967838947708073239d0 &
      +mx*(3993.83433574622979757935610692842933356144d0 &
      +mx*(14242.7295865552708506232731633468180669284d0 &
      +mx*(51394.7572916887209594591528374806790960057d0 &
      +mx*(187246.702914623152141768788258141932569037d0 &
      +mx*(687653.092375389902708761221294282367947659d0 &
      +mx*(2.54238553565398227033448846432182516906624d6 &
      +mx*(9.45378121934749027243313241962076028066811d6 &
      +mx*(3.53283630179709170835024033154326126569613d7 &
      +mx*(1.32593262383393014923560730485845833322771d8 &
      +mx*(4.99544968184054821463279808395426941549833d8 &
      +mx*(1.88840934729443872364972817525484292678543d9 &
      +mx*(7.16026753447893719179055010636502508063102d9 &
      +mx*(2.72233079469633962247554894093591262281929d10 &
        ))))))))))))))))))))
elseif(m.le.0.85d0) then
  mx=0.175d0-mc
  elb=     0.915922052601931494319853880201442948834592d0 &
      +mx*(0.294714252429483394379515488141632749820347d0 &
      +mx*(0.435776709264636140422971598963772380161131d0 &
      +mx*(1.067328246493644238508159085364429570207744d0 &
      +mx*(3.327844118563268085074646976514979307993733d0 &
      +mx*(11.90406004445092906188837729711173326621810d0 &
      +mx*(46.47838820224626393512400481776284680677096d0 &
      +mx*(192.7556002578809476962739389101964074608802d0 &
      +mx*(835.3356299261900063712302517586717381557137d0 &
      +mx*(3743.124548343029102644419963712353854902019d0 &
      +mx*(17219.07731004063094108708549153310467326395d0 &
      +mx*(80904.60401669850158353080543152212152282878d0 &
      +mx*(386808.3292751742460123683674607895217760313d0 &
      +mx*(1.876487670110449342170327796786290400635732d6 &
      +mx*(9.216559908641567755240142886998737950775910d6 &
      ))))))))))))))
  eld=     1.402200569110579095046054435635136986038164d0 &
      +mx*(2.322205897861749446534352741005347103992773d0 &
      +mx*(7.462158366466719682730245467372788273333992d0 &
      +mx*(29.43506890797307903104978364254987042421285d0 &
      +mx*(128.1590924337895775262509354898066132182429d0 &
      +mx*(591.0807036911982326384997979640812493154316d0 &
      +mx*(2830.546229607726377048576057043685514661188d0 &
      +mx*(13917.76431889392229954434840686741305556862d0 &
      +mx*(69786.10525163921228258055074102587429394212d0 &
      +mx*(355234.1420341879634781808899208309503519936d0 &
      +mx*(1.830019186413931053503912913904321703777885d6 &
      +mx*(9.519610812032515607466102200648641326190483d6 &
      +mx*(4.992086875574849453986274042758566713803723d7 &
      +mx*(2.635677009826023473846461512029006874800883d8 &
      +mx*(1.399645765120061118824228996253541612110338d9 &
      +mx*(7.469935792837635004663183580452618726280406d9 &
      +mx*(4.004155595835610574316003488168804738481448d10 &
      +mx*(2.154630668144966654449602981243932210695662d11 &
      )))))))))))))))))
else
    mx=0.125d0-mc
  elb=     0.931906061029524827613331428871579482766771d0 &
      +mx*(0.348448029538453860999386797137074571589376d0 &
      +mx*(0.666809178846938247558793864839434184202736d0 &
      +mx*(2.210769135708128662563678717558631573758222d0 &
      +mx*(9.491765048913406881414290930355300611703187d0 &
      +mx*(47.09304791027740853381457907791343619298913d0 &
      +mx*(255.9200460211233087050940506395442544885608d0 &
      +mx*(1480.029532675805407554800779436693505109703d0 &
      +mx*(8954.040904734313578374783155553041065984547d0 &
      +mx*(56052.48220982686949967604699243627330816542d0 &
      +mx*(360395.7241626000916973524840479780937869149d0 &
      +mx*(2.367539415273216077520928806581689330885103d6 &
      +mx*(1.582994957277684102454906900025484391190264d7 &
      +mx*(1.074158093278511100137056972128875270067228d8 &
      +mx*(7.380585460239595691878086073095523043390649d8 &
      +mx*(5.126022002555101496684687154904781856830296d9 &
      +mx*(3.593534065502416588712409180013118409428367d10 &
      +mx*(2.539881257612812212004146637239987308133582d11 &
      +mx*(1.808180007145359569674767150594344316702507d12 &
        ))))))))))))))))))
  eld=     1.541690112721819084362258323861459983048179d0 &
      +mx*(3.379176214579645449453938918349243359477706d0 &
      +mx*(14.94058385670236671625328259137998668324435d0 &
      +mx*(81.91773929235074880784578753539752529822986d0 &
      +mx*(497.4900546551479866036061853049402721939835d0 &
      +mx*(3205.184010234846235275447901572262470252768d0 &
      +mx*(21457.32237355321925571253220641357074594515d0 &
      +mx*(147557.0156564174712105689758692929775004292d0 &
      +mx*(1.035045290185256525452269053775273002725343d6 &
      +mx*(7.371922334832212125197513363695905834126154d6 &
      +mx*(5.314344395142401141792228169170505958906345d7 &
      +mx*(3.868823475795976312985118115567305767603128d8 &
      +mx*(2.839458401528033778425531336599562337200510d9 &
      +mx*(2.098266122943898941547136470383199468548861d10 &
      +mx*(1.559617754017662417944194874282275405676282d11 &
      +mx*(1.165096220419884791236699872205721392201682d12 &
      +mx*(8.742012983013913804987431275193291316808766d12 &
      +mx*(6.584725462672366918676967847406180155459650d13 &
      +mx*(4.976798737062434393396993620379481464465749d14 &
      +mx*(3.773018634056605404718444239040628892506293d15 &
      +mx*(2.868263194837819660109735981973458220407767d16 &
        ))))))))))))))))))))
endif

! write(*,"(1x,a20,0p5f20.15)") "(celbd) mc,b,d=",mc,elb,eld

mcold=mc
elbold=elb
eldold=eld

return
end subroutine celbd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elemental subroutine serbd(y,m,b,d)

real(dp), intent(in)  :: y,m
real(dp), intent(out) :: b,d
real(dp) ::  F1,F2,F3,F4,F10,F20,F21,F30,F31,F40,F41,F42
real(dp) ::  F5,F50,F51,F52,F6,F60,F61,F62,F63
real(dp) ::  F7,F70,F71,F72,F73,F8,F80,F81,F82,F83,F84
real(dp) ::  F9,F90,F91,F92,F93,F94
real(dp) ::  FA,FA0,FA1,FA2,FA3,FA4,FA5
real(dp) ::  FB,FB0,FB1,FB2,FB3,FB4,FB5
parameter (F10=1.d0/6.d0)
parameter (F20=3.d0/40.d0)
parameter (F21=2.d0/40.d0)
parameter (F30=5.d0/112.d0)
parameter (F31=3.d0/112.d0)
parameter (F40=35.d0/1152.d0)
parameter (F41=20.d0/1152.d0)
parameter (F42=18.d0/1152.d0)
parameter (F50=63.d0/2816.d0)
parameter (F51=35.d0/2816.d0)
parameter (F52=30.d0/2816.d0)
parameter (F60=231.d0/13312.d0)
parameter (F61=126.d0/13312.d0)
parameter (F62=105.d0/13312.d0)
parameter (F63=100.d0/13312.d0)
parameter (F70=429.d0/30720.d0)
parameter (F71=231.d0/30720.d0)
parameter (F72=189.d0/30720.d0)
parameter (F73=175.d0/30720.d0)
parameter (F80=6435.d0/557056.d0)
parameter (F81=3432.d0/557056.d0)
parameter (F82=2722.d0/557056.d0)
parameter (F83=2520.d0/557056.d0)
parameter (F84=2450.d0/557056.d0)
parameter (F90=12155.d0/1245184.d0)
parameter (F91=6435.d0/1245184.d0)
parameter (F92=5148.d0/1245184.d0)
parameter (F93=4620.d0/1245184.d0)
parameter (F94=4410.d0/1245184.d0)
parameter (FA0=46189.d0/5505024.d0)
parameter (FA1=24310.d0/5505024.d0)
parameter (FA2=19305.d0/5505024.d0)
parameter (FA3=17160.d0/5505024.d0)
parameter (FA4=16170.d0/5505024.d0)
parameter (FA5=15876.d0/5505024.d0)
parameter (FB0=88179.d0/12058624.d0)
parameter (FB1=46189.d0/12058624.d0)
parameter (FB2=36465.d0/12058624.d0)
parameter (FB3=32175.d0/12058624.d0)
parameter (FB4=30030.d0/12058624.d0)
parameter (FB5=29106.d0/12058624.d0)

real(dp) ::  A1,A2,A3,A4,A5,A6,A7,A8,A9,AA,AB
parameter (A1=3.d0/5.d0)
parameter (A2=5.d0/7.d0)
parameter (A3=7.d0/9.d0)
parameter (A4=9.d0/11.d0)
parameter (A5=11.d0/13.d0)
parameter (A6=13.d0/15.d0)
parameter (A7=15.d0/17.d0)
parameter (A8=17.d0/19.d0)
parameter (A9=19.d0/21.d0)
parameter (AA=21.d0/23.d0)
parameter (AB=23.d0/25.d0)

real(dp) ::  B1,B2,B3,B4,B5,B6,B7,B8,B9,BA,BB
real(dp) ::  D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,DA,DB
parameter (D0=1.d0/3.d0)

! write(*,*) "(serbd) y,m=",y,m

F1=F10+m*F10
F2=F20+m*(F21+m*F20)
F3=F30+m*(F31+m*(F31+m*F30))
F4=F40+m*(F41+m*(F42+m*(F41+m*F40)))
F5=F50+m*(F51+m*(F52+m*(F52+m*(F51+m*F50))))
F6=F60+m*(F61+m*(F62+m*(F63+m*(F62+m*(F61+m*F60)))))
F7=F70+m*(F71+m*(F72+m*(F73+m*(F73+m*(F72+m*(F71+m*F70))))))
F8=F80+m*(F81+m*(F82+m*(F83+m*(F84+m*(F83+m*(F82+m*(F81+m*F80)))))))
F9=F90+m*(F91+m*(F92+m*(F93+m*(F94+m*(F94+m*(F93+m*(F92+m*(F91+m*F90))))))))
FA=FA0+m*(FA1+m*(FA2+m*(FA3+m*(FA4+m*(FA5+m*(FA4+m*(FA3+m*(FA2+m*(FA1+m*FA0)))))))))
FB=FB0+m*(FB1+m*(FB2+m*(FB3+m*(FB4+m*(FB5+m*(FB5+m*(FB4+m*(FB3+m*(FB2+m*(FB1+m*FB0))))))))))

D1=F1*A1
D2=F2*A2
D3=F3*A3
D4=F4*A4
D5=F5*A5
D6=F6*A6
D7=F7*A7
D8=F8*A8
D9=F9*A9
DA=FA*AA
DB=FB*AB

d=D0+y*(D1+y*(D2+y*(D3+y*(D4+y*(D5+y*(D6+y*(D7+y*(D8+y*(D9+y*(DA+y*DB))))))))))

B1=F1-D0
B2=F2-D1
B3=F3-D2
B4=F4-D3
B5=F5-D4
B6=F6-D5
B7=F7-D6
B8=F8-D7
B9=F9-D8
BA=FA-D9
BB=FB-DA

b=1.d0+y*(B1+y*(B2+y*(B3+y*(B4+y*(B5+y*(B6+y*(B7+y*(B8+y*(B9+y*(BA+y*BB))))))))))

! write(*,*) "(serbd) b,d=",b,d

return
end subroutine serbd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine relbd(phi,phic,mc,b,d)

real(sp) ::  phi,phic,mc,b,d
real(sp) ::  m,c,x,d2,z,bc,dc,sz,v,t2

if(phi.lt.1.250) then
! write(*,*) "relsbd"
    call relsbd(sin(phi),mc,b,d)
! write(*,*) "relsbd:b,d=",b,d
else
    m=1.0-mc
! write(*,'(a10,1pe25.17)') "m=",m
    c=sin(phic)
    x=c*c
    d2=mc+m*x
    if(x.lt.0.9*d2) then
        z=c/sqrt(d2)
! write(*,*) "relsbd"
        call relsbd(z,mc,b,d)
! write(*,*) "relsbd:b,d=",b,d
        call rcelbd(mc,bc,dc)
    sz=z*sqrt(1.0-x)
    b=bc-(b-sz)
    d=dc-(d+sz)
  else
    v=mc*(1.0-x)
    if(v.lt.x*d2) then
! write(*,*) "relcbd"
            call relcbd(c,mc,b,d)
! write(*,*) "relcbd:b,d=",b,d
        else
! write(*,*) "relcbd"
            t2=(1.0-x)/d2
! write(*,*) "relcbd"
            call relcbd(sqrt(mc*t2),mc,b,d)
! write(*,*) "relcbd:b,d=",b,d
            call rcelbd(mc,bc,dc)
      sz=c*sqrt(t2)
      b=bc-(b-sz)
      d=dc-(d+sz)
    endif
  endif
endif
return
end subroutine relbd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine relsbd(s0,mc,b,d)

real(sp) ::  s0,mc,b,d
real(sp) ::  m,del,s,y,sy
real(sp) ::  yy(11),ss(11)
integer i, j

! write(*,*) "(relsbd) s0,m=",s0,m

m=1.0-mc
del=0.1888-0.0378*m ! F6  Optimum
s=s0
y=s*s
if(y.lt.del) then
  call rserbd(y,m,b,d)
  b=s*b
  d=s*y*d
  return
endif
ss(1)=s
do j=1,10
    y=y/((1.0+sqrt(1.0-y))*(1.0+sqrt(1.0-m*y)))
  yy(j+1)=y
  ss(j+1)=sqrt(y)
  if(y.lt.del) then
    goto 1
  endif
enddo
write(*,*) "(relsbd) too many iterations: s0,m=",s0,m
1 continue
! write(*,*) 'j=',j
call rserbd(y,m,b,d)
b=ss(j+1)*b
d=ss(j+1)*y*d
do i=1,j
  sy=ss(j-i+1)*yy(j-i+2)
  b=b+(b-sy)
  d=d+(d+sy)
enddo
return
end subroutine relsbd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine relcbd(c0,mc,b,dx)

real(sp) ::  c0,mc,b,dx
real(sp) ::  c,x,y,s,m,d,sy
real(sp) ::  yy(11),ss(11)
integer j,i

! write(*,*) "(relcbd) c0,mc=",c0,mc

c=c0
x=c*c
y=1.0-x
s=sqrt(y)
if(x.gt.0.1) then
! write(*,*) "(relcbd) elsbd"
  call relsbd(s,mc,b,dx)
    return
endif
m=1.0-mc
ss(1)=s
! write(*,*) "(relcbd) j,y,s=",1,y,ss(1)
do j=1,10
    d=sqrt(mc+m*x)
  x=(c+d)/(1.0+d)
  y=1.0-x
  yy(j+1)=y
  ss(j+1)=sqrt(y)
! write(*,*) "(relcbd) j,y,s=",j+1,yy(j+1),ss(j+1)
    if(x.gt.0.1) then
        goto 1
  endif
  c=sqrt(x)
enddo
write(*,*) "(relcbd) too many iterations: c0,mc=",c0,mc
1 continue
s=ss(j+1)
call relsbd(s,mc,b,dx)
do i=1,j
  sy=ss(j-i+1)*yy(j-i+2)
  b=b+(b-sy)
  dx=dx+(dx+sy)
enddo
return
end subroutine relcbd





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rcelbd(mc,elb,eld)
!
! Single precision associate complete elliptic integral of the second kind
!
!     Reference: T. Fukushima, (2011), Math. Comp., 80, 1725-1743
!        "Precise and Fast Computation of General Complete Elliptic Integral
!         of Second Kind"
!
!     Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!     Inputs: mc = complementary parameter, 0 < mc <= 1
!
!     Outputs: elb = B(m), eld = D(m)
!
real(sp) ::  mc,elb,eld
real(sp) ::  m,nome,dkkc,dddc,mx,kkc,logq2,elk,dele,elk1,delb

real(sp) ::  PIQ,PIHALF
parameter (PIQ=0.78539816)
parameter (PIHALF=1.57079633)

real(sp) ::  mcold,elbold,eldold
!save mcold,elbold,eldold

real(sp) ::  Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8
parameter (Q1=1.0/16.0,Q2=1.0/32.0,Q3=21.0/1024.0)
parameter (Q4=31.0/2048.0,Q5=6257.0/524288.0)
parameter (Q6=10293.0/1048576.0,Q7=279025.0/33554432.0)
parameter (Q8=483127.0/67108864.0)

real(sp) ::  K1,K2,K3,K4
parameter (K1=1.0/4.0)
parameter (K2=9.0/64.0)
parameter (K3=25.0/256.0)
parameter (K4=1225.0/16384.0)

real(sp) ::  B1,B2,B3,B4,B5
parameter (B1=1.0/2.0)
parameter (B2=1.0/16.0)
parameter (B3=3.0/128.0)
parameter (B4=25.0/2048.0)
parameter (B5=245.0/32768.0)

real(sp) ::  D1,D2,D3,D4,D5
parameter (D1=1.0/2.0)
parameter (D2=3.0/16.0)
parameter (D3=15.0/128.0)
parameter (D4=175.0/2048.0)
parameter (D5=2205.0/32768.0)

!logical first/.TRUE./

!if(first) then
!    first=.FALSE.
  mcold=1.0
  elbold=PIQ
  eldold=PIQ
!endif
m=1.0-mc
if(abs(mc-mcold).lt.1.19e-7*mc) then
  elb=elbold
  eld=eldold
elseif(m.lt.1.19e-7) then
    elb=PIQ
  eld=PIQ
elseif(mc.lt.1.19e-7) then
    elb=1.0
  eld=0.386294361-0.5*log(mc)
elseif(mc.lt.0.1) then
    nome=mc*(Q1+mc*(Q2+mc*(Q3+mc*(Q4+mc*(Q5+mc*(Q6+mc*(Q7+mc*Q8))))))) 
    if(mc.lt.0.01) then
        dkkc=mc*(K1+mc*(K2+mc*(K3+mc*K4)))
        dddc=mc*(D1+mc*(D2+mc*(D3+mc*D4)))
    else
        mx=mc-0.05

! (K'-1)/(pi/2)

      dkkc=0.0128642566+mx*(0.2648342989+mx*(0.1564757379+mx*( &
           0.1142614608+mx*(0.0920272442+mx*(0.0784321883+mx*( &
           0.0693526014))))))

! (K'-E')/(pi/2)

    dddc=0.0254839544+mx*(0.5196738432+mx*(0.2064495111+mx*( &
             0.1361095213+mx*(0.1045801404+mx*(0.0867461292+mx*( &
             0.0753638027))))))
    endif
  kkc=1.0+dkkc
  logq2=-0.5*log(nome)
  elk=kkc*logq2
  dele=-dkkc/kkc+logq2*dddc
  elk1=elk-1.0
  delb=(dele-mc*elk1)/m
  elb=1.0+delb
  eld=elk1-delb
elseif(m.le.0.01) then
  elb=PIHALF*(B1+m*(B2+m*(B3+m*(B4+m*B5))))
  eld=PIHALF*(D1+m*(D2+m*(D3+m*(D4+m*D5))))
elseif(m.le.0.1) then
  mx=0.95-mc
  elb=0.7904014136+mx*(0.1020062662+mx*(0.03987839556+mx*( &
        0.02173713638+mx*(0.01396097977+mx*(0.009892518823)))))
  eld=0.8006020402+mx*(0.3139944778+mx*(0.2059131187+mx*( &
        0.1577443465+mx*(0.1305950773+mx*(0.1133084745)))))
elseif(m.le.0.2) then
  mx=0.85-mc
  elb=0.8010240645+mx*(0.1106953445+mx*(0.0473487467+mx*( &
        0.0284843673+mx*(0.0202778114+mx*(0.0159650059)))))
  eld=0.8342326678+mx*(0.3604952816+mx*(0.2623796641+mx*( &
        0.2237239445+mx*(0.2064478118+mx*(0.1998094409)))))
elseif(m.le.0.3) then
  mx=0.75-mc
  elb=0.8125977729+mx*(0.1211096179+mx*(0.05729337683+mx*( &
        0.03850945160+mx*(0.03078343030+mx*(0.02729056493+mx*( &
        0.02591636929))))))
  eld=0.8731525819+mx*(0.4206222307+mx*(0.3442310616+mx*( &
        0.3311330218+mx*(0.3452772851+mx*(0.3779453222+mx*( &
        0.4273780125))))))
elseif(m.le.0.4) then
  mx=0.65-mc
  elb=0.8253235580+mx*(0.1338621161+mx*(0.07101129360+mx*( &
        0.05417847742+mx*(0.04945174495+mx*(0.05022219622+mx*( &
        0.05474291317))))))
  eld=0.9190270392+mx*(0.5010021593+mx*(0.4688312706+mx*( &
        0.5177142278+mx*(0.6208433913+mx*(0.7823643938+mx*( &
        1.0191145351))))))
elseif(m.le.0.5) then
  mx=0.55-mc
  elb=0.8394795703+mx*(0.1499164403+mx*(0.09083193582+mx*( &
        0.08034703348+mx*(0.08563844050+mx*(0.1019547259+mx*( &
        0.1305748115))))))
  eld=0.9744043665+mx*(0.6132468054+mx*(0.6710966695+mx*( &
        0.8707276202+mx*(1.2295422312+mx*(1.8266059675+mx*( &
        2.8069345310+mx*(4.4187893291)))))))
elseif(m.le.0.6) then
  mx=0.45-mc
  elb=0.8554696152+mx*(0.1708960727+mx*(0.1213352290+mx*( &
        0.1282018836+mx*(0.1646872815+mx*(0.2374189087+mx*( &
        0.3692081047))))))
  eld=1.0434552951+mx*(0.7796257219+mx*(1.0297423609+mx*( &
        1.6220372234+mx*(2.7879895312+mx*(5.0483814874+mx*( &
        9.4632776119+mx*(18.181489949+mx*(35.580980591))))))))
elseif(m.le.0.7) then
  mx=0.35-mc
  elb=0.8739200618+mx*(0.1998140575+mx*(0.1727696159+mx*( &
        0.2281069133+mx*(0.3704681411+mx*(0.6792712529+mx*( &
        1.3480084967+mx*(2.8276709769)))))))
  eld=1.1336783366+mx*(1.0486431737+mx*(1.7534650412+mx*( &
        3.5231827268+mx*(7.7494764138+mx*(17.986450056+mx*( &
        43.255916346+mx*(106.68153445+mx*(268.09848657))))))))
elseif(m.le.0.8) then
  mx=0.25-mc
  elb=0.8959028209+mx*(0.2431400038+mx*(0.2730818756+mx*( &
        0.4862800075+mx*(1.0827474372+mx*(2.7434452910+mx*( &
        7.5558178287+mx*(22.051940825+mx*(67.156406447+mx*( &
        211.27225379)))))))))
  eld=1.2606128266+mx*(1.5486656381+mx*(3.5536694119+mx*( &
        9.9004446761+mx*(30.320566617+mx*(98.180258659+mx*( &
        329.77101043+mx*(1136.6559897+mx*(3993.8343357+mx*( &
        14242.729587+mx*(51394.757292))))))))))
elseif(m.le.0.85) then
  mx=0.175-mc
  elb=0.9159220526+mx*(0.2947142524+mx*(0.4357767093+mx*( &
        1.0673282465+mx*(3.3278441186+mx*(11.904060044+mx*( &
        46.478388202+mx*(192.75560026)))))))
  eld=1.4022005691+mx*(2.3222058979+mx*(7.4621583665+mx*( &
        29.435068908+mx*(128.15909243+mx*(591.08070369+mx*( &
        2830.5462296+mx*(13917.764319+mx*(69786.105252))))))))
else
    mx=0.125-mc
  elb=0.9319060610+mx*(0.3484480295+mx*(0.6668091788+mx*( &
        2.2107691357+mx*(9.4917650489+mx*(47.093047910+mx*( &
        255.92004602+mx*(1480.0295327+mx*(8954.0409047+mx*( &
        56052.482210)))))))))
  eld=1.5416901127+mx*(3.3791762146+mx*(14.940583857+mx*( &
        81.917739292+mx*(497.49005466+mx*(3205.1840102+mx*( &
        21457.322374+mx*(147557.01566+mx*(1.0350452902e6+mx*( &
        7.3719223348e6+mx*(5.3143443951e7))))))))))
endif

! write(*,"(1x,a20,0p5f15.7)") "(rcelbd) mc,b,d=",mc,elb,eld

mcold=mc
elbold=elb
eldold=eld

end subroutine rcelbd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rserbd(y,m,b,d)

real(sp) ::  y,m,b,d
real(sp) ::  F1,F2,F3,F4,F10,F20,F21,F30,F31,F40,F41,F42
real(sp) ::  F5,F50,F51,F52,F6,F60,F61,F62,F63
parameter (F10=1.0/6.0)
parameter (F20=3.0/40.0)
parameter (F21=2.0/40.0)
parameter (F30=5.0/112.0)
parameter (F31=3.0/112.0)
parameter (F40=35.0/1152.0)
parameter (F41=20.0/1152.0)
parameter (F42=18.0/1152.0)
parameter (F50=63.0/2816.0)
parameter (F51=35.0/2816.0)
parameter (F52=30.0/2816.0)
parameter (F60=231.0/13312.0)
parameter (F61=126.0/13312.0)
parameter (F62=105.0/13312.0)
parameter (F63=100.0/13312.0)

real(dp) ::  A1,A2,A3,A4,A5,A6
parameter (A1=3.0/5.0)
parameter (A2=5.0/7.0)
parameter (A3=7.0/9.0)
parameter (A4=9.0/11.0)
parameter (A5=11.0/13.0)
parameter (A6=13.0/15.0)

real(dp) ::  B1,B2,B3,B4,B5,B6
real(dp) ::  D0,D1,D2,D3,D4,D5,D6
parameter (D0=1.0/3.0)

! write(*,*) "(rserbd) y,m=",y,m

F1=F10+m*F10
F2=F20+m*(F21+m*F20)
F3=F30+m*(F31+m*(F31+m*F30))
F4=F40+m*(F41+m*(F42+m*(F41+m*F40)))
F5=F50+m*(F51+m*(F52+m*(F52+m*(F51+m*F50))))
F6=F60+m*(F61+m*(F62+m*(F63+m*(F62+m*(F61+m*F60)))))

D1=F1*A1
D2=F2*A2
D3=F3*A3
D4=F4*A4
D5=F5*A5
D6=F6*A6

d=real(D0+y*(D1+y*(D2+y*(D3+y*(D4+y*(D5+y*D6))))), sp)

B1=F1-D0
B2=F2-D1
B3=F3-D2
B4=F4-D3
B5=F5-D4
B6=F6-D5

b=real(1.0+y*(B1+y*(B2+y*(B3+y*(B4+y*(B5+y*B6))))), sp)

! write(*,*) "(rserbd) b,d=",b,d

end subroutine rserbd



end module elliptic_integral_mod
