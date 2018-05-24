!+
! Implements Bengtsson's summation RDT formula's for first and second order resonant driving terms.
!-

module srdt_mod

use twiss_and_track_mod

implicit none

type summation_rdt_struct
  complex(rp) h20001
  complex(rp) h00201
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

character(6), parameter :: srdt_first(7) = [ 'h20001', 'h00201', 'h21000', 'h30000', 'h10110', 'h10020', 'h10200' ]
character(6), parameter :: srdt_second(14) = [ 'h31000', 'h40000', 'h20110', 'h11200', 'h20020', 'h20200', 'h00310', \
                                               'h00400', 'h22000', 'h00220', 'h11110', 'nux_Jx', 'nuy_Jy', 'nux_Jy' ]
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

subroutine srdt_calc (lat, srdt, order, n_slices_gen_opt, n_slices_sext_opt)

implicit none

type ele_special_struct
  real(rp) k1l,k2l,s,l
  real(rp) eta_a, beta_a, beta_b, phi_a, phi_b
  logical good_k2
end type

type(lat_struct) lat
type(summation_rdt_struct) srdt
integer order
integer, optional :: n_slices_sext_opt
integer, optional :: n_slices_gen_opt

type(ele_special_struct), allocatable :: eles(:)
type(ele_struct) :: elei
type(coord_struct), allocatable :: co(:)

real(rp) pinux, pinuy
real(rp) k2l, k1, k2
real(rp) slice_len
real(rp) ns, sj, sl
real(rp) dmux, dmuy, prod, sqrtprod, sgn

integer pass, w, i, j
integer n_slices_sext
integer n_slices_gen

logical good_ele, good_k2

!

n_slices_gen = integer_option(10, n_slices_gen_opt)
n_slices_sext = integer_option(20, n_slices_sext_opt)

pinux = lat%ele(lat%n_ele_track)%a%phi/2.0d0
pinuy = lat%ele(lat%n_ele_track)%b%phi/2.0d0

!pre-cache k1, k2 moments to speed up term summations in the loop that follows.
do pass=1,2
  w = 0
  do i=1,lat%n_ele_track
    if( lat%ele(i)%key == wiggler$ ) cycle
    if( lat%ele(i)%key == multipole$ ) then
      if(attribute_index(lat%ele(i),'K2L').ne.0) then
        if(abs(value_of_attribute(lat%ele(i),'K2L')).gt.1e-8 .and. value_of_attribute(lat%ele(i), 'K2L') .ne. real_garbage$) then
          w = w + 1
          if(pass == 2) then
            k2l = value_of_attribute(lat%ele(i),'K2L')
            k2l = k2l / 2.0  ! Convention shown in Eqn. 8 of Bengsston paper: moments divided by n.
            eles(w)%good_k2 = .true.
            eles(w)%k1l = 0.0d0
            eles(w)%k2l = k2l
            eles(w)%s = lat%ele(i)%s
            eles(w)%eta_a = lat%ele(i)%a%eta
            eles(w)%beta_a = lat%ele(i)%a%beta
            eles(w)%beta_b = lat%ele(i)%b%beta
            eles(w)%phi_a = lat%ele(i)%a%phi
            eles(w)%phi_b = lat%ele(i)%b%phi
          endif
        endif
      endif
    elseif(value_of_attribute(lat%ele(i), 'l') .gt. 1e-6) then
      good_ele = .false.
      good_k2 = .false.
      if(attribute_index(lat%ele(i),'K1').ne.0 .or. attribute_index(lat%ele(i),'K2').ne.0) then
        k1 = 0.0d0
        k2 = 0.0d0
        if(abs(value_of_attribute(lat%ele(i),'K1')).gt.1e-8 .and. value_of_attribute(lat%ele(i), 'K1') .ne. real_garbage$) then
          k1 = value_of_attribute(lat%ele(i),'K1')
          good_ele = .true.
        endif
        if(abs(value_of_attribute(lat%ele(i),'K2')).gt.1d-8 .and. value_of_attribute(lat%ele(i), 'K2') .ne. real_garbage$) then
          k2 = value_of_attribute(lat%ele(i),'K2')
          good_ele = .true.
          good_k2 = .true.
        endif
        if(good_ele) then
          k2 = k2 / 2.0  ! Convention shown in Eqn. 8 of Bengsston paper: moments divided by n.
          if( lat%ele(i)%key == sextupole$ ) then
            ns = n_slices_sext
          else
            ns = n_slices_gen
          endif
          slice_len = lat%ele(i)%value(l$) / ns
          do j=1,ns
            w = w + 1
            if(pass == 2) then
              if(i .gt. 1) then
                sj = lat%ele(i-1)%s + slice_len*j
              else
                sj = slice_len*j
              endif
              call twiss_and_track_at_s(lat,sj,elei,co)
              eles(w)%good_k2 = good_k2
              eles(w)%k1l = k1 * slice_len
              eles(w)%k2l = k2 * slice_len
              eles(w)%s = sj
              eles(w)%eta_a = elei%a%eta
              eles(w)%beta_a = elei%a%beta
              eles(w)%beta_b = elei%b%beta
              eles(w)%phi_a = elei%a%phi
              eles(w)%phi_b = elei%b%phi
            endif
          enddo
        endif
      endif
    endif
  enddo
  if(pass==1) allocate(eles(w))
enddo

! Calculate first order terms.
srdt%h20001 = sum((eles(:)%k1l-2.0*eles(:)%k2l*eles(:)%eta_a)*eles(:)%beta_a * exp( i_imag*2.0d0*eles(:)%phi_a ))
srdt%h00201 = sum((eles(:)%k1l-2.0*eles(:)%k2l*eles(:)%eta_a)*eles(:)%beta_b * exp( i_imag*2.0d0*eles(:)%phi_b ))
srdt%h21000 = sum(eles(:)%k2l*eles(:)%beta_a**(3./2.) * exp( i_imag*1.0d0*eles(:)%phi_a ))
srdt%h30000 = sum(eles(:)%k2l*eles(:)%beta_a**(3./2.) * exp( i_imag*3.0d0*eles(:)%phi_a ))
srdt%h10110 = sum(eles(:)%k2l*eles(:)%beta_a**(1./2.)*eles(:)%beta_b * exp( i_imag*1.0d0*eles(:)%phi_a ))
srdt%h10020 = sum(eles(:)%k2l*eles(:)%beta_a**(1./2.)*eles(:)%beta_b * exp( i_imag*(eles(:)%phi_a-2.0*eles(:)%phi_b) ))
srdt%h10200 = sum(eles(:)%k2l*eles(:)%beta_a**(1./2.)*eles(:)%beta_b * exp( i_imag*(eles(:)%phi_a+2.0*eles(:)%phi_b) ))
srdt%h21000 = srdt%h21000 / -8.0
srdt%h30000 = srdt%h30000 / -24.0
srdt%h10110 = srdt%h10110 /  4.0
srdt%h10020 = srdt%h10020 /  8.0
srdt%h10200 = srdt%h10200 /  8.0
srdt%h20001 = srdt%h20001 /  8.0
srdt%h00201 = srdt%h00201 / -8.0


!Calculate second order terms
srdt%h22000 = 0
srdt%h00220 = 0
srdt%nux_Jy = 0
srdt%nuy_Jy = 0
srdt%nux_Jx = 0
srdt%h31000 = 0
srdt%h40000 = 0
srdt%h20110 = 0
srdt%h11200 = 0
srdt%h20020 = 0
srdt%h20200 = 0
srdt%h00310 = 0
srdt%h00400 = 0
if(order .ge. 2) then
  do i=1,w
    if( eles(i)%good_k2 ) then
      do j=1,w
        if( eles(j)%good_k2 ) then
          dmux = eles(i)%phi_a-eles(j)%phi_a
          dmuy = eles(i)%phi_b-eles(j)%phi_b
          prod = eles(i)%k2l*eles(j)%k2l

          srdt%nux_Jx = srdt%nux_Jx + prod * (eles(i)%beta_a*eles(j)%beta_a)**(3./2.) * &
                     ( 3*cos(abs(dmux)-pinux)/sin(pinux) + cos(3*abs(dmux)-3*pinux)/sin(3*pinux) )
          srdt%nuy_Jy = srdt%nuy_Jy + prod * sqrt(eles(i)%beta_a*eles(j)%beta_a)*eles(i)%beta_b*eles(j)%beta_b * &
                     ( 4*cos(abs(dmux)-pinux)/sin(pinux) + &
                         cos(abs(dmux+2*dmuy)-(pinux+2*pinuy))/sin((pinux+2*pinuy)) + &
                         cos(abs(dmux-2*dmuy)-(pinux-2*pinuy))/sin((pinux-2*pinuy)) )
          srdt%nux_Jy = srdt%nux_Jy + prod * sqrt(eles(i)%beta_a*eles(j)%beta_a)*eles(i)%beta_b * &
                     ( 2*eles(j)%beta_a*cos(abs(dmux)-pinux)/sin(pinux) - &
                         eles(j)%beta_b*cos(abs(dmux+2*dmuy)-(pinux+2*pinuy))/sin((pinux+2*pinuy)) + &
                         eles(j)%beta_b*cos(abs(dmux-2*dmuy)-(pinux-2*pinuy))/sin((pinux-2*pinuy)) )

          if(i .ne. j) then
            if( i < j ) then
              sgn = 1
            elseif( i > j ) then
              sgn = -1
            endif
            sqrtprod = sqrt(eles(i)%beta_a*eles(j)%beta_a)
            srdt%h31000 = srdt%h31000 + sgn*prod*(eles(i)%beta_a*eles(j)%beta_a)**(3./2.) * exp( i_imag*(3.0*eles(i)%phi_a-eles(j)%phi_a) ) !2Qx
            srdt%h40000 = srdt%h40000 + sgn*prod*(eles(i)%beta_a*eles(j)%beta_a)**(3./2.) * exp( i_imag*(3.0*eles(i)%phi_a+eles(j)%phi_a) ) !4Qx
            srdt%h20110 = srdt%h20110 + sgn*prod*sqrtprod*eles(i)%beta_b * &
                              ( eles(j)%beta_a*(exp( -i_imag*(eles(i)%phi_a-3.0*eles(j)%phi_a) ) - exp( i_imag*(eles(i)%phi_a+eles(j)%phi_a) ) ) + &
                            2.0*eles(j)%beta_b*exp( i_imag*(eles(i)%phi_a+eles(j)%phi_a+2.0*eles(i)%phi_b-2.0*eles(j)%phi_b) ) )  !2Qx
            srdt%h11200 = srdt%h11200 + sgn*prod*sqrtprod*eles(i)%beta_b * &
                              ( eles(j)%beta_a*( exp(-i_imag*(eles(i)%phi_a-eles(j)%phi_a-2.0*eles(i)%phi_b) ) - &
                                                 exp( i_imag*(eles(i)%phi_a-eles(j)%phi_a+2.0*eles(i)%phi_b) ) ) + &
                            2.0*eles(j)%beta_b*( exp( i_imag*(eles(i)%phi_a-eles(j)%phi_a+2.0*eles(i)%phi_b) ) + &
                                                 exp(-i_imag*(eles(i)%phi_a-eles(j)%phi_a-2.0*eles(i)%phi_b) ) ) ) !2Qy
            srdt%h20020 = srdt%h20020 + sgn*prod*sqrtprod*eles(i)%beta_b * &
                                       ( eles(j)%beta_a*exp( -i_imag*(eles(i)%phi_a-3.0*eles(j)%phi_a+2.0*eles(i)%phi_b) ) - &
                                       (eles(j)%beta_a+4.0*eles(j)%beta_b)*exp(i_imag*(eles(i)%phi_a+eles(j)%phi_a-2.0*eles(i)%phi_b)) )! 2Qx+2Qy
            srdt%h20200 = srdt%h20200 + sgn*prod*sqrtprod*eles(i)%beta_b * &
                                       ( eles(j)%beta_a*exp( -i_imag*(eles(i)%phi_a-3.0*eles(j)%phi_a-2.0*eles(i)%phi_b) ) - &
                                       (eles(j)%beta_a-4.0*eles(j)%beta_b)*exp(i_imag*(eles(i)%phi_a+eles(j)%phi_a+2.0*eles(i)%phi_b)) )! 2Qx-2Qy
            srdt%h00310 = srdt%h00310 + sgn*prod*sqrtprod*eles(i)%beta_b*eles(j)%beta_b * &
                                       ( exp( i_imag*(eles(i)%phi_a-eles(j)%phi_a+2.0*eles(i)%phi_b) ) - &
                                         exp(-i_imag*(eles(i)%phi_a-eles(j)%phi_a-2.0*eles(i)%phi_b) ) )! 2Qy
            srdt%h00400 = srdt%h00400 + sgn*prod*sqrtprod*eles(i)%beta_b*eles(j)%beta_b * &
                                        exp(i_imag*(eles(i)%phi_a-eles(j)%phi_a+2.0*eles(i)%phi_b+2.0*eles(j)%phi_b) )  !4Qy
          endif
        endif
      enddo
    endif
  enddo
  srdt%nux_Jx = srdt%nux_Jx / -16.0 / pi
  srdt%nuy_Jy = srdt%nuy_Jy / 8.0 / pi
  srdt%nux_Jy = srdt%nux_Jy / -16.0 / pi
  srdt%h31000 = srdt%h31000 * i_imag / 32.0
  srdt%h40000 = srdt%h40000 * i_imag / 64.0
  srdt%h20110 = srdt%h20110 * i_imag / 32.0
  srdt%h11200 = srdt%h11200 * i_imag / 32.0
  srdt%h20020 = srdt%h20020 * i_imag / 64.0
  srdt%h20200 = srdt%h20200 * i_imag / 64.0
  srdt%h00310 = srdt%h00310 * i_imag / 32.0
  srdt%h00400 = srdt%h00400 * i_imag / 64.0
endif

end subroutine

end module
