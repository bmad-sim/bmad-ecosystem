!+
! Implements Bengtsson's summation RDT formula's for first and second order resonant driving terms.
!-

module srdt_mod

use twiss_and_track_mod

implicit none

type summation_rdt_struct
  complex(rp) h11001
  complex(rp) h00111
  complex(rp) h20001
  complex(rp) h00201
  complex(rp) h10002
  complex(rp) h21000
  complex(rp) h30000
  complex(rp) h10110
  complex(rp) h10020
  complex(rp) h10200
  !2nd order in K2 moments
  complex(rp) h31000 
  complex(rp) h40000 
  complex(rp) h20110 
  complex(rp) h11200 
  complex(rp) h20020 
  complex(rp) h20200 
  complex(rp) h00310 
  complex(rp) h00400 
  complex(rp) h22000 
  complex(rp) h00220 
  complex(rp) h11110 
  real(rp) nux_Jx
  real(rp) nuy_Jy
  real(rp) nux_Jy
end type

character(6), parameter :: srdt_first(8) = [ 'h20001', 'h00201', 'h10002', 'h21000', 'h30000', 'h10110', 'h10020', 'h10200' ]
character(6), parameter :: srdt_second(14) = [ 'h31000', 'h40000', 'h20110', 'h11200', 'h20020', 'h20200', 'h00310', \
                                               'h00400', 'h22000', 'h00220', 'h11110', 'nux_Jx', 'nuy_Jy', 'nux_Jy' ]

type sliced_eles_struct
  integer ix
  real(rp) k1l,k2l,s,l
  real(rp) eta_a, beta_a, beta_b, phi_a, phi_b
  logical good_k2
  type(summation_rdt_struct) srdt
  complex(rp) ea, eb, e2a, e2b, e3a
end type

private make_slices

contains

!+ 
! Subroutine srdt_calc(lat, srdt, order, n_slices_gen_opt, n_slices_sext_opt)
!
! Calculate summation RDT terms up to order=1 or order=2 while slicing sextupoles 
! n_slices_sext_opt times and all other elements n_slices_gen_opt times.
!
! These formulas are documented in "The Sextupole Scheme for the Swiss Light Source (SLS): An Analytic Approach"
! by Johan Bengtsson.  SLS Note 9/97.
!
! The 2nd order formulas are documented in "Second-order driving terms due to sextupoles and
! chromatic effects of quadrupoles" by Chun-xi Wang.  AOP-TN-2009-020.
!
! Input:
!   lat                -- lat_struct: lattice with Twiss parameters calculated.
!   order              -- integer: 1 to calculate only first order terms.  2 to also calculate 2nd order terms.
!   n_slices_gen_opt   -- integer, optional: number of times to slice elements other than sextupoles.  Default is 10.
!   n_slices_sext_opt  -- integer, optional: nubmer of times to slice sextupoles.  Default is 20.
!
! Output:
!   srdt               -- summation_rdt_struct: contains complex RDT strengths.
!-

subroutine srdt_calc (lat, srdt_sums, order, n_slices_sext_opt, n_slices_gen_opt, per_ele_out)

implicit none

type(lat_struct) lat
type(summation_rdt_struct) srdt_sums
integer order
integer, optional :: n_slices_sext_opt
integer, optional :: n_slices_gen_opt

type(sliced_eles_struct), allocatable :: eles(:)
type(summation_rdt_struct), allocatable, optional :: per_ele_out(:)

real(rp) pinux, pinuy
real(rp) dmux, dmuy, prod, sqrtprod, sgn

integer pass, w, i, j, ns
integer n_slices_sext
integer n_slices_gen
integer nK2

n_slices_gen = integer_option(10, n_slices_gen_opt)
n_slices_sext = integer_option(20, n_slices_sext_opt)

pinux = lat%ele(lat%n_ele_track)%a%phi/2.0d0
pinuy = lat%ele(lat%n_ele_track)%b%phi/2.0d0

call make_slices(lat, eles, n_slices_gen, n_slices_sext)

w = size(eles)

! Calculate first order terms.
eles(:)%srdt%h11001 = (eles(:)%k1l-2.0*eles(:)%k2l*eles(:)%eta_a)*eles(:)%beta_a
eles(:)%srdt%h00111 = (eles(:)%k1l-2.0*eles(:)%k2l*eles(:)%eta_a)*eles(:)%beta_b
eles(:)%srdt%h20001 = (eles(:)%k1l-2.0*eles(:)%k2l*eles(:)%eta_a)*eles(:)%beta_a * eles(:)%e2a
eles(:)%srdt%h00201 = (eles(:)%k1l-2.0*eles(:)%k2l*eles(:)%eta_a)*eles(:)%beta_b * eles(:)%e2b
eles(:)%srdt%h10002 = (eles(:)%k1l-eles(:)%k2l*eles(:)%eta_a)*eles(:)%eta_a*sqrt(eles(:)%beta_a) * eles(:)%ea
eles(:)%srdt%h21000 = eles(:)%k2l*eles(:)%beta_a**(3./2.) * eles(:)%ea
eles(:)%srdt%h30000 = eles(:)%k2l*eles(:)%beta_a**(3./2.) * eles(:)%e3a
eles(:)%srdt%h10110 = eles(:)%k2l*eles(:)%beta_a**(1./2.)*eles(:)%beta_b * eles(:)%ea
eles(:)%srdt%h10020 = eles(:)%k2l*eles(:)%beta_a**(1./2.)*eles(:)%beta_b * eles(:)%ea/eles(:)%e2b
eles(:)%srdt%h10200 = eles(:)%k2l*eles(:)%beta_a**(1./2.)*eles(:)%beta_b * eles(:)%ea*eles(:)%e2b

srdt_sums%h11001 = sum(eles(1:size(eles))%srdt%h11001) 
srdt_sums%h00111 = sum(eles(1:size(eles))%srdt%h00111) 
srdt_sums%h20001 = sum(eles(1:size(eles))%srdt%h20001) 
srdt_sums%h00201 = sum(eles(1:size(eles))%srdt%h00201) 
srdt_sums%h10002 = sum(eles(1:size(eles))%srdt%h10002) 
srdt_sums%h21000 = sum(eles(1:size(eles))%srdt%h21000) 
srdt_sums%h30000 = sum(eles(1:size(eles))%srdt%h30000) 
srdt_sums%h10110 = sum(eles(1:size(eles))%srdt%h10110) 
srdt_sums%h10020 = sum(eles(1:size(eles))%srdt%h10020) 
srdt_sums%h10200 = sum(eles(1:size(eles))%srdt%h10200) 

srdt_sums%h11001 =  srdt_sums%h11001 / 4.0
srdt_sums%h00111 = -srdt_sums%h00111 / 4.0
srdt_sums%h21000 = -srdt_sums%h21000 / 8.0
srdt_sums%h30000 = -srdt_sums%h30000 / 24.0
srdt_sums%h10110 =  srdt_sums%h10110 / 4.0
srdt_sums%h10020 =  srdt_sums%h10020 / 8.0
srdt_sums%h10200 =  srdt_sums%h10200 / 8.0
srdt_sums%h20001 =  srdt_sums%h20001 / 8.0
srdt_sums%h00201 = -srdt_sums%h00201 / 8.0
srdt_sums%h10002 =  srdt_sums%h10002 / 2.0


!Calculate second order terms
srdt_sums%h22000 = 0
srdt_sums%h00220 = 0
srdt_sums%h11110 = 0
srdt_sums%nux_Jy = 0
srdt_sums%nuy_Jy = 0
srdt_sums%nux_Jx = 0
srdt_sums%h31000 = 0
srdt_sums%h40000 = 0
srdt_sums%h20110 = 0
srdt_sums%h11200 = 0
srdt_sums%h20020 = 0
srdt_sums%h20200 = 0
srdt_sums%h00310 = 0
srdt_sums%h00400 = 0
if(order .ge. 2) then
  do i=1,w
    if( eles(i)%good_k2 ) then
      do j=1,w
        if( eles(j)%good_k2 ) then
          dmux = eles(i)%phi_a-eles(j)%phi_a
          dmuy = eles(i)%phi_b-eles(j)%phi_b
          prod = eles(i)%k2l*eles(j)%k2l

          !srdt_sums%nux_Jx = srdt_sums%nux_Jx + prod * (eles(i)%beta_a*eles(j)%beta_a)**(3./2.)*(3.0*res(i,j,1,0)+res(i,j,3,0))
          !srdt_sums%nux_Jy = srdt_sums%nux_Jy + prod * (eles(i)%beta_a*eles(j)%beta_a)**(1./2.)*eles(i)%beta_b*( &
          !                                              2.0*eles(j)%beta_a*res(i,j,1,0)-eles(j)%beta_b*res(i,j,1,2)+eles(j)%beta_b*res(i,j,1,-2))
          !srdt_sums%nuy_Jy = srdt_sums%nuy_Jy + prod * (eles(i)%beta_a*eles(j)%beta_a)**(1./2.)*eles(i)%beta_b*eles(j)%beta_b*( &
          !                                              4.0*res(i,j,1,0)+res(i,j,1,2)+res(i,j,1,-2) )
          if(i /= j) then
            if( i < j ) then
              sgn = 1
            elseif( i > j ) then
              sgn = -1
            endif
            sqrtprod = sqrt(eles(i)%beta_a*eles(j)%beta_a)
            srdt_sums%h22000 = srdt_sums%h22000 + sgn*prod*(eles(i)%beta_a*eles(j)%beta_a)**(3./2.) * ( eles(i)%e3a/eles(j)%e3a + 3.0*eles(i)%ea/eles(j)%ea )
            srdt_sums%h00220 = srdt_sums%h00220 + sgn*prod*(eles(i)%beta_a*eles(j)%beta_a)**(1./2.) * (eles(i)%beta_b*eles(j)%beta_b) * ( &
                                        eles(i)%ea/eles(j)%ea*eles(i)%e2b/eles(j)%e2b + 4.0*eles(i)%ea/eles(j)%ea - eles(j)%ea/eles(i)%ea*eles(i)%e2b/eles(j)%e2b )
            srdt_sums%h11110 = srdt_sums%h11110 + sgn*prod*sqrt(eles(i)%beta_a*eles(j)%beta_a)*eles(i)%beta_b * ( &
                                                  eles(j)%beta_a*( eles(j)%ea/eles(i)%ea - eles(i)%ea/eles(j)%ea ) + &
                                                  eles(j)%beta_b*( eles(i)%ea/eles(j)%ea*eles(i)%e2b/eles(j)%e2b + eles(j)%ea/eles(i)%ea*eles(i)%e2b/eles(j)%e2b ) )

            srdt_sums%h31000 = srdt_sums%h31000 + sgn*prod*(eles(i)%beta_a*eles(j)%beta_a)**(3./2.) * eles(i)%e3a/eles(j)%ea !2Qx
            srdt_sums%h40000 = srdt_sums%h40000 + sgn*prod*(eles(i)%beta_a*eles(j)%beta_a)**(3./2.) * eles(i)%e3a*eles(j)%ea !4Qx
            srdt_sums%h20110 = srdt_sums%h20110 + sgn*prod*sqrtprod*eles(i)%beta_b * &
                                                  ( eles(j)%beta_a*( eles(j)%e3a/eles(i)%ea - eles(i)%ea*eles(j)%ea ) + &
                                                    2.0*eles(j)%beta_b*eles(i)%ea*eles(j)%ea*eles(i)%e2b/eles(j)%e2b ) 
            srdt_sums%h11200 = srdt_sums%h11200 + sgn*prod*sqrtprod*eles(i)%beta_b * &
                              ( eles(j)%beta_a*( eles(j)%ea/eles(i)%ea*eles(i)%e2b - eles(i)%ea/eles(j)%ea*eles(i)%e2b ) + &
                            2.0*eles(j)%beta_b*( eles(i)%ea/eles(j)%ea*eles(i)%e2b + eles(j)%ea/eles(i)%ea*eles(i)%e2b ) )
            srdt_sums%h20020 = srdt_sums%h20020 + sgn*prod*sqrtprod*eles(i)%beta_b * &
                                       ( eles(j)%beta_a*eles(j)%e3a/eles(i)%ea/eles(i)%e2b - (eles(j)%beta_a+4.0*eles(j)%beta_b)*eles(i)%ea*eles(j)%ea/eles(i)%e2b )
            srdt_sums%h20200 = srdt_sums%h20200 + sgn*prod*sqrtprod*eles(i)%beta_b * &
                                       ( eles(j)%beta_a*eles(j)%e3a/eles(i)%ea*eles(i)%e2b - (eles(j)%beta_a-4.0*eles(j)%beta_b)*eles(i)%ea*eles(j)%ea*eles(i)%e2b )
            srdt_sums%h00310 = srdt_sums%h00310 + sgn*prod*sqrtprod*eles(i)%beta_b*eles(j)%beta_b * &
                                       ( eles(i)%ea/eles(j)%ea*eles(i)%e2b - eles(j)%ea/eles(i)%ea*eles(i)%e2b )
            srdt_sums%h00400 = srdt_sums%h00400 + sgn*prod*sqrtprod*eles(i)%beta_b*eles(j)%beta_b * eles(i)%ea/eles(j)%ea*eles(i)%e2b*eles(j)%e2b
          endif
        endif
      enddo
    endif
  enddo
  srdt_sums%nux_Jx =  srdt_sums%nux_Jx / (-16.0) / pi
  srdt_sums%nux_Jy =  srdt_sums%nux_Jy / 8.0 / pi
  srdt_sums%nuy_Jy =  srdt_sums%nuy_Jy / (-16.0) / pi
  srdt_sums%h22000 =  srdt_sums%h22000 * i_imag / 64.0
  srdt_sums%h11110 =  srdt_sums%h11110 * i_imag / 16.0
  srdt_sums%h00220 =  srdt_sums%h00220 * i_imag / 64.0
  srdt_sums%h31000 =  srdt_sums%h31000 * i_imag / 32.0
  srdt_sums%h40000 =  srdt_sums%h40000 * i_imag / 64.0
  srdt_sums%h20110 =  srdt_sums%h20110 * i_imag / 32.0
  srdt_sums%h11200 =  srdt_sums%h11200 * i_imag / 32.0
  srdt_sums%h20020 =  srdt_sums%h20020 * i_imag / 64.0
  srdt_sums%h20200 =  srdt_sums%h20200 * i_imag / 64.0
  srdt_sums%h00310 =  srdt_sums%h00310 * i_imag / 32.0
  srdt_sums%h00400 =  srdt_sums%h00400 * i_imag / 64.0
endif

if(present(per_ele_out)) then
  allocate(per_ele_out(lat%n_ele_track))
  do i=1,lat%n_ele_track
    per_ele_out(i)%h20001 = sum(eles(:)%srdt%h20001, mask=eles(:)%ix==i)
    per_ele_out(i)%h00201 = sum(eles(:)%srdt%h00201, mask=eles(:)%ix==i)
    per_ele_out(i)%h10002 = sum(eles(:)%srdt%h10002, mask=eles(:)%ix==i)
    per_ele_out(i)%h21000 = sum(eles(:)%srdt%h21000, mask=eles(:)%ix==i)
    per_ele_out(i)%h30000 = sum(eles(:)%srdt%h30000, mask=eles(:)%ix==i)
    per_ele_out(i)%h10110 = sum(eles(:)%srdt%h10110, mask=eles(:)%ix==i)
    per_ele_out(i)%h10020 = sum(eles(:)%srdt%h10020, mask=eles(:)%ix==i)
    per_ele_out(i)%h10200 = sum(eles(:)%srdt%h10200, mask=eles(:)%ix==i)

    per_ele_out(i)%h31000 = sum(eles(:)%srdt%h31000, mask=eles(:)%ix==i)
    per_ele_out(i)%h40000 = sum(eles(:)%srdt%h40000, mask=eles(:)%ix==i)
    per_ele_out(i)%h20110 = sum(eles(:)%srdt%h20110, mask=eles(:)%ix==i)
    per_ele_out(i)%h11200 = sum(eles(:)%srdt%h11200, mask=eles(:)%ix==i)
    per_ele_out(i)%h20020 = sum(eles(:)%srdt%h20020, mask=eles(:)%ix==i)
    per_ele_out(i)%h20200 = sum(eles(:)%srdt%h20200, mask=eles(:)%ix==i)
    per_ele_out(i)%h00310 = sum(eles(:)%srdt%h00310, mask=eles(:)%ix==i)
    per_ele_out(i)%h00400 = sum(eles(:)%srdt%h00400, mask=eles(:)%ix==i)
    per_ele_out(i)%h22000 = sum(eles(:)%srdt%h22000, mask=eles(:)%ix==i)
    per_ele_out(i)%h00220 = sum(eles(:)%srdt%h00220, mask=eles(:)%ix==i)
    per_ele_out(i)%h11110 = sum(eles(:)%srdt%h11110, mask=eles(:)%ix==i)
  enddo
endif
  contains
  function res(i,j,nx,ny) result(r)
    implicit none

    real(rp) r
    integer i, j, nx, ny
    real(rp) arg_px, arg_py, arg_ux, arg_uy

    arg_px = nx*(eles(i)%phi_a-eles(j)%phi_a)
    arg_py = ny*(eles(i)%phi_b-eles(j)%phi_b)
    arg_ux = nx*pinux
    arg_uy = ny*pinuy

    r = cos(abs(arg_px+arg_py)-arg_ux-arg_uy)/sin(arg_ux+arg_uy)
  end function
end subroutine

!+
! Subroutine make_slices(lat,eles,n_slices_gen,n_slices_sext)
!
! Private subroutine in srdt_mod.f90.  Given lat, returns eles, which
! contains lattice functions needed for driving term calculations on
! sub-element spacing.  Also pre-computed some quantities to speed
! up driving term calculations.
!
! Input:
!   lat                -- lat_struct: lattice with Twiss parameters calculated.
!   n_slices_gen_opt   -- integer: number of times to slice elements other than sextupoles.
!   n_slices_sext_opt  -- integer: nubmer of times to slice sextupoles.
! Output:
!   eles               -- sliced_eles_struct: lattice parameters needed for driving term calculations.
!-
subroutine make_slices(lat,eles,n_slices_gen,n_slices_sext)

  implicit none

  type(lat_struct) lat
  type(sliced_eles_struct), allocatable :: eles(:)
  integer n_slices_gen, n_slices_sext

  type(coord_struct), allocatable :: co(:)
  type(ele_struct) :: elei

  integer i, j, w
  integer ns, pass

  logical good_ele, good_k2

  real(rp) k2l, k1, k2, slice_len, sj

  do pass=1,2
    w = 0
    do i=1,lat%n_ele_track
      if( lat%ele(i)%key == wiggler$ ) cycle
      if( lat%ele(i)%key == multipole$ ) then
        if(attribute_index(lat%ele(i),'K2L') /= 0) then
          if(abs(value_of_attribute(lat%ele(i),'K2L')) > 1e-8 .and. value_of_attribute(lat%ele(i), 'K2L') /= real_garbage$) then
            w = w + 1
            if(pass == 2) then
              k2l = value_of_attribute(lat%ele(i),'K2L')
              k2l = k2l / 2.0  ! Convention shown in Eqn. 8 of Bengsston paper: moments divided by n.
              eles(w)%ix = i
              eles(w)%good_k2 = .true.
              eles(w)%k1l = 0.0d0
              eles(w)%k2l = k2l
              eles(w)%l = -1
              eles(w)%s = lat%ele(i)%s
              eles(w)%eta_a  = lat%ele(i)%a%eta
              eles(w)%beta_a = lat%ele(i)%a%beta
              eles(w)%beta_b = lat%ele(i)%b%beta
              eles(w)%phi_a  = lat%ele(i)%a%phi
              eles(w)%phi_b  = lat%ele(i)%b%phi
              eles(w)%ea = exp(i_imag*lat%ele(i)%a%phi)
              eles(w)%eb = exp(i_imag*lat%ele(i)%b%phi)
              eles(w)%e2a = exp(i_imag*2.0d0*lat%ele(i)%a%phi)
              eles(w)%e2b = exp(i_imag*2.0d0*lat%ele(i)%b%phi)
              eles(w)%e3a = exp(i_imag*3.0d0*lat%ele(i)%a%phi)
            endif
          endif
        endif
      elseif(value_of_attribute(lat%ele(i), 'l')  >  1e-6) then
        good_ele = .false.
        good_k2 = .false.
        if(attribute_index(lat%ele(i),'K1')/=0 .or. attribute_index(lat%ele(i),'K2')/=0) then
          k1 = 0.0d0
          k2 = 0.0d0
          if(abs(value_of_attribute(lat%ele(i),'K1')) > 1e-8 .and. value_of_attribute(lat%ele(i), 'K1') /= real_garbage$) then
            k1 = value_of_attribute(lat%ele(i),'K1')
            good_ele = .true.
          endif
          if(abs(value_of_attribute(lat%ele(i),'K2')) > 1d-8 .and. value_of_attribute(lat%ele(i), 'K2') /= real_garbage$) then
            k2 = value_of_attribute(lat%ele(i),'K2')
            good_ele = .true.
            good_k2 = .true.
          endif
          if(good_ele) then
            k2 = k2 / 2.0  ! Convention shown in Eqn. 8 of Bengsston paper: moments divided by n.
            if (lat%ele(i)%key == sextupole$) then
              ns = n_slices_sext
            else
              ns = n_slices_gen
            endif
            slice_len = lat%ele(i)%value(l$) / ns
            do j = 1, ns
              w = w + 1
              if(pass == 2) then
                if(i > 1) then
                  sj = lat%ele(i-1)%s + slice_len*j
                else
                  sj = slice_len*j
                endif
                call twiss_and_track_at_s(lat,sj,elei,co)
                eles(w)%ix = i
                eles(w)%good_k2 = good_k2
                eles(w)%k1l = k1 * slice_len
                eles(w)%k2l = k2 * slice_len
                eles(w)%l = slice_len
                eles(w)%s = sj
                eles(w)%eta_a = elei%a%eta
                eles(w)%beta_a =elei%a%beta
                eles(w)%beta_b =elei%b%beta
                eles(w)%phi_a = elei%a%phi
                eles(w)%phi_b = elei%b%phi
                eles(w)%ea = exp(i_imag*elei%a%phi)
                eles(w)%eb = exp(i_imag*elei%b%phi)
                eles(w)%e2a = exp(i_imag*2.0d0*elei%a%phi)
                eles(w)%e2b = exp(i_imag*2.0d0*elei%b%phi)
                eles(w)%e3a = exp(i_imag*3.0d0*elei%a%phi)
              endif
            enddo
          endif
        endif
      endif
    enddo
    if(pass==1) allocate(eles(w))
  enddo
end subroutine

!+
! Subroutine srdt_lsq_solution(lat, ls_soln, n_slices_sext_opt, n_slices_gen_opt)
!
! Given lat, finds K2 moments that set the chromaticity and zeros-out the real
! and complex parts of the first order driving terms, that minimizes the sum of the squares
! of the K2 moments.  i.e. the weakest sextupole scheme that sets chromaticity
! and zeros out the first order terms.
!
! This subroutine assumes that Nsext > 18.
!
! Input:
!   lat                -- lat_struct: lattice with Twiss parameters calculated.
!   n_slices_gen_opt   -- integer, optional: number of times to slice elements other than sextupoles.  Default is 10.
!   n_slices_sext_opt  -- integer, optional: nubmer of times to slice sextupoles.  Default is 20.
! Output:
!   ls_soln(1:n_ele_track)  -- real(rp): contains K2 for 
!-
subroutine srdt_lsq_solution(lat, ls_soln, n_slices_sext_opt, n_slices_gen_opt)

  implicit none

  type(lat_struct) lat
  real(rp) ls_soln(:)
  integer, optional :: n_slices_sext_opt
  integer, optional :: n_slices_gen_opt

  type(sliced_eles_struct), allocatable :: eles(:)
  type(sliced_eles_struct), allocatable :: K2eles(:)

  integer i,j,w,nK2,nmag
  integer, allocatable :: k2mask(:), mags(:)

  real(rp), allocatable :: A(:,:), Ap(:,:)
  real(rp) B(18)
  real(rp), allocatable :: ls_soln_sliced(:)

  logical, allocatable :: mask(:)

  integer n_slices_sext
  integer n_slices_gen
 
  n_slices_gen = integer_option(10, n_slices_gen_opt)
  n_slices_sext = integer_option(20, n_slices_sext_opt)

  call make_slices(lat, eles, n_slices_gen, n_slices_sext)
  w = size(eles)

  nK2 = count(eles(:)%good_k2 .eqv. .true.)
  allocate(k2eles(nK2))
  !k2eles is a subset of eles containing only those slices with a valid sextupole moment.
  !Assumption here is that all slices with valid sextupole moment are in fact sextupoles which
  !can be adjusted to manipulate nonlinerities.
  j=0
  do i=1,w
    if(eles(i)%good_k2) then
      j = j + 1
      k2eles(j) = eles(i)
    endif
  enddo
  !Count the number of sextupoles represented in k2eles.  Note that accurate calculation of
  !RDTs requires sextupole slicing.  However, the K2 of each slice cannot be manipulated independently
  !to minimize the RDTs.
  nmag=0
  do i=1,size(ls_soln)
    if(count(k2eles(:)%ix==i) .gt. 0) nmag = nmag + 1
  enddo
  allocate(A(18,nmag))
  allocate(Ap(nmag,18))
  allocate(ls_soln_sliced(nmag))
  allocate(mask(nmag))
  allocate(mags(nmag))
  j=0
  do i=1,size(ls_soln)
    mask = k2eles(:)%ix==i
    if(count(mask) .gt. 0) then
      j = j + 1
      mags(j) = i
      A(1 ,j) = sum(k2eles(:)%l*(-2.0)*k2eles(:)%eta_a*k2eles(:)%beta_a,mask=mask)   ! chrom_x
      A(2 ,j) = sum(k2eles(:)%l*(-2.0)*k2eles(:)%eta_a*k2eles(:)%beta_b,mask=mask)   ! chrom_y
      A(3 ,j) = sum(k2eles(:)%l* real(-2.0*k2eles(:)%eta_a*k2eles(:)%beta_a * k2eles(:)%e2a),mask=mask)  !h20001
      A(4 ,j) = sum(k2eles(:)%l*aimag(-2.0*k2eles(:)%eta_a*k2eles(:)%beta_a * k2eles(:)%e2a),mask=mask)
      A(5 ,j) = sum(k2eles(:)%l* real( 2.0*k2eles(:)%eta_a*k2eles(:)%beta_b * k2eles(:)%e2b),mask=mask)  !h00201
      A(6 ,j) = sum(k2eles(:)%l*aimag( 2.0*k2eles(:)%eta_a*k2eles(:)%beta_b * k2eles(:)%e2b),mask=mask)
      A(7 ,j) = sum(k2eles(:)%l* real(k2eles(:)%eta_a*k2eles(:)%eta_a*sqrt(k2eles(:)%beta_a) * k2eles(:)%ea),mask=mask) !h10002
      A(8 ,j) = sum(k2eles(:)%l*aimag(k2eles(:)%eta_a*k2eles(:)%eta_a*sqrt(k2eles(:)%beta_a) * k2eles(:)%ea),mask=mask)
      A(9 ,j) = sum(k2eles(:)%l* real(k2eles(:)%beta_a**(3./2.) * k2eles(:)%ea),mask=mask)  !h21000
      A(10,j) = sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(3./2.) * k2eles(:)%ea),mask=mask)
      A(11,j) = sum(k2eles(:)%l* real(k2eles(:)%beta_a**(3./2.) * k2eles(:)%e3a),mask=mask)  !h30000
      A(12,j) = sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(3./2.) * k2eles(:)%e3a),mask=mask)
      A(13,j) = sum(k2eles(:)%l* real(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea),mask=mask)  !h10110
      A(14,j) = sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea),mask=mask)
      A(15,j) = sum(k2eles(:)%l* real(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea/k2eles(:)%e2b),mask=mask) !h10020
      A(16,j) = sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea/k2eles(:)%e2b),mask=mask)
      A(17,j) = sum(k2eles(:)%l* real(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea*k2eles(:)%e2b),mask=mask) !h10200
      A(18,j) = sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea*k2eles(:)%e2b),mask=mask)
    endif
  enddo
  call make_pseudoinverse(A,Ap)
  B(:) = 0.0d0
  B(1) = -4.0+sum(eles(:)%k1l*eles(:)%beta_a)
  B(2) = 4.0-sum(-eles(:)%k1l*eles(:)%beta_b)
  B(3) = sum(real(eles(:)%k1l*eles(:)%beta_a*eles(:)%e2a))
  B(4) = sum(aimag(eles(:)%k1l*eles(:)%beta_a*eles(:)%e2a))
  B(5) = sum(real(-eles(:)%k1l*eles(:)%beta_b*eles(:)%e2b))
  B(6) = sum(aimag(-eles(:)%k1l*eles(:)%beta_b*eles(:)%e2b))
  B(7) = sum(real(-eles(:)%k1l*eles(:)%eta_a*sqrt(eles(:)%beta_a) * eles(:)%ea))
  B(8) = sum(aimag(-eles(:)%k1l*eles(:)%eta_a*sqrt(eles(:)%beta_a) * eles(:)%ea))
  ls_soln_sliced = matmul(Ap,-1.0d0*B) 
  ls_soln = 0
  do i=1,nmag
    ls_soln(mags(i)) = 2.0*ls_soln_sliced(i)
  enddo
end subroutine

end module
