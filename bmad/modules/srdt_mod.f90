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
end type

! srdt_first is first order terms list, srdt_second is second order terms list.

character(6), parameter :: srdt_first(10) = ['h20001', 'h00201', 'h10002', 'h21000', 'h30000', 'h10110', 'h10020', &
                                             'h10200', 'h11001', 'h00111']
character(6), parameter :: srdt_second(11) = ['h31000', 'h40000', 'h20110', 'h11200', 'h20020', 'h20200', 'h00310', &
                                              'h00400', 'h22000', 'h00220', 'h11110']

type sliced_eles_struct
  integer ix
  real(rp) k1l, k2l, s, l
  real(rp) eta_a, beta_a, beta_b, phi_a, phi_b
  logical good_k2
  type(summation_rdt_struct) srdt
  complex(rp) ea, eb, e2a, e2b, e3a
end type

private make_slices

contains

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!+ 
! Subroutine srdt_calc(lat, srdt_sums, order, n_slices_gen_opt, n_slices_sxt_opt)
!
! Calculate summation RDT terms up to order=1 or order=2 while slicing sextupoles 
! n_slices_sxt_opt times and all other elements n_slices_gen_opt times.
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
!   n_slices_sxt_opt   -- integer, optional: nubmer of times to slice sextupoles.  Default is 20.
!
! Output:
!   srdt_sums           -- summation_rdt_struct: contains complex RDT strengths.
!-

subroutine srdt_calc (lat, srdt_sums, order, n_slices_gen_opt, n_slices_sxt_opt, per_ele_out)

implicit none

type(lat_struct) lat
type(summation_rdt_struct) srdt_sums
integer order
integer, optional :: n_slices_gen_opt
integer, optional :: n_slices_sxt_opt

type(sliced_eles_struct), allocatable :: eles(:)
type(summation_rdt_struct), allocatable, optional :: per_ele_out(:)

real(rp) dmux, dmuy, prod, sqrtprod, sgn

integer pass, w, i, j, ns
integer n_slices_sext
integer n_slices_gen
integer nK2

!

n_slices_gen = integer_option(10, n_slices_gen_opt)
n_slices_sext = integer_option(20, n_slices_sxt_opt)

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
srdt_sums%h31000 = 0
srdt_sums%h40000 = 0
srdt_sums%h20110 = 0
srdt_sums%h11200 = 0
srdt_sums%h20020 = 0
srdt_sums%h20200 = 0
srdt_sums%h00310 = 0
srdt_sums%h00400 = 0
if (order .ge. 2) then
  do i=1, w
    if ( eles(i)%good_k2 ) then
      do j=1, w
        if ( eles(j)%good_k2 ) then
          prod = eles(i)%k2l*eles(j)%k2l

          if (i /= j) then
            if ( i < j ) then
              sgn = 1
            elseif ( i > j ) then
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

if (present(per_ele_out)) then
  allocate(per_ele_out(lat%n_ele_track))
  do i=1, lat%n_ele_track
    per_ele_out(i)%h20001 = sum(eles(:)%srdt%h20001, mask=eles(:)%ix==i) / 8.0
    per_ele_out(i)%h00201 = sum(eles(:)%srdt%h00201, mask=eles(:)%ix==i) / 8.0
    per_ele_out(i)%h10002 = sum(eles(:)%srdt%h10002, mask=eles(:)%ix==i) / 8.0
    per_ele_out(i)%h21000 = sum(eles(:)%srdt%h21000, mask=eles(:)%ix==i) / 8.0
    per_ele_out(i)%h30000 = sum(eles(:)%srdt%h30000, mask=eles(:)%ix==i) / 24.0
    per_ele_out(i)%h10110 = sum(eles(:)%srdt%h10110, mask=eles(:)%ix==i) / 4.0
    per_ele_out(i)%h10020 = sum(eles(:)%srdt%h10020, mask=eles(:)%ix==i) / 8.0
    per_ele_out(i)%h10200 = sum(eles(:)%srdt%h10200, mask=eles(:)%ix==i) / 8.0

    per_ele_out(i)%h31000 = sum(eles(:)%srdt%h31000, mask=eles(:)%ix==i) / 32.0
    per_ele_out(i)%h40000 = sum(eles(:)%srdt%h40000, mask=eles(:)%ix==i) / 64.0
    per_ele_out(i)%h20110 = sum(eles(:)%srdt%h20110, mask=eles(:)%ix==i) / 32.0
    per_ele_out(i)%h11200 = sum(eles(:)%srdt%h11200, mask=eles(:)%ix==i) / 32.0
    per_ele_out(i)%h20020 = sum(eles(:)%srdt%h20020, mask=eles(:)%ix==i) / 64.0
    per_ele_out(i)%h20200 = sum(eles(:)%srdt%h20200, mask=eles(:)%ix==i) / 64.0
    per_ele_out(i)%h00310 = sum(eles(:)%srdt%h00310, mask=eles(:)%ix==i) / 32.0
    per_ele_out(i)%h00400 = sum(eles(:)%srdt%h00400, mask=eles(:)%ix==i) / 64.0
    per_ele_out(i)%h22000 = sum(eles(:)%srdt%h22000, mask=eles(:)%ix==i) / 64.0
    per_ele_out(i)%h00220 = sum(eles(:)%srdt%h00220, mask=eles(:)%ix==i) / 64.0
    per_ele_out(i)%h11110 = sum(eles(:)%srdt%h11110, mask=eles(:)%ix==i) / 16.0
  enddo
endif

end subroutine

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!+
!-

subroutine srdt_calc_with_cache (lat, srdt_sums, order, n_slices_gen_opt, n_slices_sxt_opt, cache, per_ele_out)

implicit none

type(lat_struct) lat
type(summation_rdt_struct) srdt_sums
integer order
integer, optional :: n_slices_gen_opt
integer, optional :: n_slices_sxt_opt
complex(rp), allocatable :: cache(:,:,:)
type(summation_rdt_struct), allocatable, optional :: per_ele_out(:)

type(sliced_eles_struct), allocatable :: eles(:)
integer, allocatable :: ixK2(:)
real(rp) sgn
real(rp) prod

integer pass, w, i, j, ns
integer xi, xj, wK2
integer n_slices_gen, n_slices_sxt

!

n_slices_gen = integer_option(10, n_slices_gen_opt)
n_slices_sxt = integer_option(20, n_slices_sxt_opt)

call make_srdt_cache (lat, n_slices_gen, n_slices_sxt, eles, cache)

if (.not. allocated(cache)) then
  call srdt_calc(lat, srdt_sums, order, n_slices_gen_opt, n_slices_sxt_opt, per_ele_out)
  return
endif

!

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
srdt_sums%h31000 = 0
srdt_sums%h40000 = 0
srdt_sums%h20110 = 0
srdt_sums%h11200 = 0
srdt_sums%h20020 = 0
srdt_sums%h20200 = 0
srdt_sums%h00310 = 0
srdt_sums%h00400 = 0

if (order >= 2) then
  wK2 = count(eles(:)%good_k2)
  allocate(ixK2(wK2))
  ixK2 = pack([(i, i=1, w)], eles(:)%good_k2)
  do xj=1, wK2
    do xi=1, wK2
      prod = eles(ixK2(xi))%k2l*eles(ixK2(xj))%k2l
      srdt_sums%h22000 = srdt_sums%h22000 + prod*cache(1, xi, xj)
      srdt_sums%h00220 = srdt_sums%h00220 + prod*cache(2, xi, xj)
      srdt_sums%h11110 = srdt_sums%h11110 + prod*cache(3, xi, xj)
      srdt_sums%h31000 = srdt_sums%h31000 + prod*cache(4, xi, xj)
      srdt_sums%h40000 = srdt_sums%h40000 + prod*cache(5, xi, xj)
      srdt_sums%h20110 = srdt_sums%h20110 + prod*cache(6, xi, xj)
      srdt_sums%h11200 = srdt_sums%h11200 + prod*cache(7, xi, xj)
      srdt_sums%h20020 = srdt_sums%h20020 + prod*cache(8, xi, xj)
      srdt_sums%h20200 = srdt_sums%h20200 + prod*cache(9, xi, xj)
      srdt_sums%h00310 = srdt_sums%h00310 + prod*cache(10, xi, xj)
      srdt_sums%h00400 = srdt_sums%h00400 + prod*cache(11, xi, xj)
    enddo
  enddo
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

if (present(per_ele_out)) then
  allocate(per_ele_out(lat%n_ele_track))
  do i=1, lat%n_ele_track
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

end subroutine srdt_calc_with_cache 

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!+
!-

subroutine make_srdt_cache(lat, n_slices_gen, n_slices_sxt, eles, cache)

implicit none

type(lat_struct) lat
integer n_slices_gen, n_slices_sxt
complex(rp), allocatable :: cache(:, :, :)

type(sliced_eles_struct), allocatable :: eles(:)
real(rp) sqrtprod
integer ierr
integer i, j, sgn
integer xi, xj
integer w, wK2

character(*), parameter :: r_name = 'make_srdt_cache'

call make_slices(lat, eles, n_slices_gen, n_slices_sxt)

w = size(eles)
wK2 = count(eles(:)%good_K2)

if (.not. allocated(cache)) then
  write(*, '(a, f15.3, a)') "Allocating approximately ", wK2*wK2*11.*2.*rp/(1024.**3), " GB for SRDT cache."
  allocate(cache(11, wK2, wK2), stat=ierr)
  if (ierr /= 0) then
    call out_io (s_warn$, r_name, 'Insufficient memory available to make SRDT cache for specified n_slices_gen and n_slices_sxt.', &
                                  'Will proceed without cache.  Calculation may be much slower.', &
                                  'Consider reducing n_slices_gen and n_slices_sxt.')
    return
  endif 
elseif (size(cache, 2) < wK2) then
  deallocate(cache)
  allocate(cache(11, wK2, wK2))
endif

cache = (0.0d0, 0.0d0)
xi = 0
do i=1, size(eles)
  if (eles(i)%good_k2) then
    xi = xi + 1
    xj = 0
    do j=1, size(eles)
      if ( eles(j)%good_k2 ) then
        xj = xj + 1
        if (i /= j) then
          if ( i < j ) then
            sgn = 1
          elseif ( i > j ) then
            sgn = -1
          endif
          ! cache(1 , :, :)  h22000
          ! cache(2 , :, :)  h00220
          ! cache(3 , :, :)  h11110
          ! cache(4 , :, :)  h31000
          ! cache(5 , :, :)  h40000
          ! cache(6 , :, :)  h20110
          ! cache(7 , :, :)  h11200
          ! cache(8 , :, :)  h20020
          ! cache(9 , :, :)  h20200
          ! cache(10, :, :)  h00310
          ! cache(11, :, :)  h00400

          sqrtprod = sqrt(eles(i)%beta_a*eles(j)%beta_a)
          cache(1, xi, xj)  = sgn*(eles(i)%beta_a*eles(j)%beta_a)**(3./2.) * ( eles(i)%e3a/eles(j)%e3a + 3.0*eles(i)%ea/eles(j)%ea )
          cache(2, xi, xj)  = sgn*(eles(i)%beta_a*eles(j)%beta_a)**(1./2.) * (eles(i)%beta_b*eles(j)%beta_b) * ( &
                          eles(i)%ea/eles(j)%ea*eles(i)%e2b/eles(j)%e2b + 4.0*eles(i)%ea/eles(j)%ea - eles(j)%ea/eles(i)%ea*eles(i)%e2b/eles(j)%e2b )
          cache(3, xi, xj)  = sgn*sqrt(eles(i)%beta_a*eles(j)%beta_a)*eles(i)%beta_b * ( &
                          eles(j)%beta_a*( eles(j)%ea/eles(i)%ea - eles(i)%ea/eles(j)%ea ) + &
                          eles(j)%beta_b*( eles(i)%ea/eles(j)%ea*eles(i)%e2b/eles(j)%e2b + eles(j)%ea/eles(i)%ea*eles(i)%e2b/eles(j)%e2b ) )
          cache(4, xi, xj)  = sgn*(eles(i)%beta_a*eles(j)%beta_a)**(3./2.) * eles(i)%e3a/eles(j)%ea !2Qx
          cache(5, xi, xj)  = sgn*(eles(i)%beta_a*eles(j)%beta_a)**(3./2.) * eles(i)%e3a*eles(j)%ea !4Qx
          cache(6, xi, xj)  = sgn*sqrtprod*eles(i)%beta_b * &
                          ( eles(j)%beta_a*( eles(j)%e3a/eles(i)%ea - eles(i)%ea*eles(j)%ea ) + &
                          2.0*eles(j)%beta_b*eles(i)%ea*eles(j)%ea*eles(i)%e2b/eles(j)%e2b ) 
          cache(7, xi, xj)  = sgn*sqrtprod*eles(i)%beta_b * &
                          ( eles(j)%beta_a*( eles(j)%ea/eles(i)%ea*eles(i)%e2b - eles(i)%ea/eles(j)%ea*eles(i)%e2b ) + &
                          2.0*eles(j)%beta_b*( eles(i)%ea/eles(j)%ea*eles(i)%e2b + eles(j)%ea/eles(i)%ea*eles(i)%e2b ) )
          cache(8, xi, xj)  = sgn*sqrtprod*eles(i)%beta_b * &
                          ( eles(j)%beta_a*eles(j)%e3a/eles(i)%ea/eles(i)%e2b - (eles(j)%beta_a+4.0*eles(j)%beta_b)*eles(i)%ea*eles(j)%ea/eles(i)%e2b )
          cache(9, xi, xj)  = sgn*sqrtprod*eles(i)%beta_b * &
                          ( eles(j)%beta_a*eles(j)%e3a/eles(i)%ea*eles(i)%e2b - (eles(j)%beta_a-4.0*eles(j)%beta_b)*eles(i)%ea*eles(j)%ea*eles(i)%e2b )
          cache(10, xi, xj) = sgn*sqrtprod*eles(i)%beta_b*eles(j)%beta_b * &
                          ( eles(i)%ea/eles(j)%ea*eles(i)%e2b - eles(j)%ea/eles(i)%ea*eles(i)%e2b )
          cache(11, xi, xj) = sgn*sqrtprod*eles(i)%beta_b*eles(j)%beta_b * eles(i)%ea/eles(j)%ea*eles(i)%e2b*eles(j)%e2b
        endif
      endif
    enddo
  endif
enddo

end subroutine make_srdt_cache

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!+
! Subroutine make_slices(lat, eles, n_slices_gen, n_slices_sext, var_indexes)
!
! Private subroutine in srdt_mod.f90.  Given lat, returns eles, which
! contains lattice functions needed for driving term calculations on
! sub-element spacing.  Also pre-computed some quantities to speed
! up driving term calculations.
!
! All elements with a K1 and/or K2 are sliced according to n_slices_gen and 
! n_slices_sext put into eles.
!
! Input:
!   lat                -- lat_struct: lattice with Twiss parameters calculated.
!   n_slices_gen_opt   -- integer: number of times to slice elements other than sextupoles.
!   n_slices_sxt_opt   -- integer: nubmer of times to slice sextupoles.
!   var_indexes(:)     -- integer: make sure that good_K2 slices are made for these elements.
!
! Output:
!   eles               -- sliced_eles_struct: lattice parameters needed for driving term calculations.
!-
subroutine make_slices(lat, eles, n_slices_gen, n_slices_sext, var_indexes)

implicit none

type(lat_struct) lat
type(sliced_eles_struct), allocatable :: eles(:)
integer n_slices_gen, n_slices_sext
integer, optional :: var_indexes(:)

type(coord_struct), allocatable :: co(:)
type(ele_struct) :: elei

integer i, j, w
integer ns, pass
integer ix_pole_max
integer, allocatable :: var_indexes_use(:)

logical good_k2

real(rp) k2l, k1l, k1, k2, slice_len, sj
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

if (present(var_indexes)) then
  allocate(var_indexes_use(size(var_indexes)))
  var_indexes_use = var_indexes
else
  allocate(var_indexes_use(1))
  var_indexes_use = -1
endif

do pass=1, 2
  w = 0
  do i=1, lat%n_ele_track
    knl = 0.0
    tilt = 0.0
    call multipole_ele_to_kt(lat%ele(i), .true., ix_pole_max, knl, tilt, magnetic$, include_kicks$)
    if ( (ix_pole_max .ge. 1) .or. any(i==var_indexes_use) ) then
      k1l = knl(1)
      if ( ((ix_pole_max .ge. 2) .and. (abs(knl(2)) .gt. 1.0d-8)) .or. any(i==var_indexes_use) ) then
        k2l = knl(2) / 2.0  ! Convention shown in Eqn. 8 of Bengsston paper: moments divided by n.
        good_k2 = .true.
      else
        k2l = 0.0
        good_k2 = .false.
      endif

      if ( lat%ele(i)%key == multipole$ ) then
        w = w + 1
        if (pass == 2) then
          eles(w)%k2l = k2l
          eles(w)%k1l = k1l
          eles(w)%ix = i
          eles(w)%good_k2 = good_k2
          eles(w)%l = -1  !zero length
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
      else
        k1 = k1l / lat%ele(i)%value(l$)
        k2 = k2l / lat%ele(i)%value(l$)
        if (lat%ele(i)%key == sextupole$) then
          ns = n_slices_sext
        else
          ns = n_slices_gen
        endif
        slice_len = lat%ele(i)%value(l$) / ns
        do j = 1, ns
          w = w + 1
          if (pass == 2) then
            if (i > 1) then
              sj = lat%ele(i-1)%s + slice_len*j
            else
              sj = slice_len*j
            endif
            call twiss_and_track_at_s(lat, sj, elei, co)
            eles(w)%k2l = k2 * slice_len
            eles(w)%k1l = k1 * slice_len
            eles(w)%ix = i
            eles(w)%good_k2 = good_k2
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
  enddo
  if (pass==1) allocate(eles(w))
enddo

end subroutine make_slices

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!+
! Subroutine srdt_lsq_solution(lat, var_indexes, ls_soln, n_slices_gen_opt, n_slices_sxt_opt, 
!                                                     chrom_set_x_opt, chrom_set_y_opt, weight_in)
!
! Given lat, finds K2 moments that set the chromaticity and zeros-out the real
! and complex parts of the first order driving terms, that minimizes the sum of the squares
! of the K2 moments.  i.e. the weakest sextupole scheme that sets chromaticity
! and zeros out the first order terms.
!
! Note:  This subroutine does not, in its present form, work well with knobs, overlays, or in lattices where
!        multiple elements have the same name.
!
! This subroutine assumes that Nsext > 18.
!
! Input:
!   lat                -- lat_struct: lattice with Twiss parameters calculated.
!   var_indexes(:)     -- integer: indexes in lat%ele that are K2 variables.  Must be sorted smallest index to largest index.
!   n_slices_gen_opt   -- integer, optional: number of times to slice elements other than sextupoles.  Default is 10.
!   n_slices_sxt_opt   -- integer, optional: nubmer of times to slice sextupoles.  Default is 20.
!   chrom_set_x_opt    -- real(rp), optional: what to set x chromaticity to.  Default zero.
!   chrom_set_y_opt    -- real(rp), optional: what to set y chromaticity to.  Default zero.
!   weight_in(10)      -- real(rp), optional: moment weights. Terms are:
!                              [wgt_chrom_x, wgt_chrom_y, wgt_h20001, wgt_h00201, wgt_h10002, 
!                              wgt_h21000, wgt_h30000, wgt_h10110, wgt_h10020, wgt_h10200, 
!                          If present, any terms equal to zero are given default values which is 
!                          1.0e4 for wgt_chrom_x and wgt_chrom_y and is 1.0 for everything else.
!
! Output:
!   ls_soln(1:size(var_indexes))  -- real(rp): contains K2 for the indexes in var_indexes
!-

subroutine srdt_lsq_solution(lat, var_indexes, ls_soln, n_slices_gen_opt, n_slices_sxt_opt, &
                                                                 chrom_set_x_opt, chrom_set_y_opt, weight_in)

implicit none

type(lat_struct) lat
integer var_indexes(:)
real(rp), allocatable :: ls_soln(:)
integer, optional :: n_slices_sxt_opt
integer, optional :: n_slices_gen_opt
real(rp), optional :: chrom_set_x_opt, chrom_set_y_opt, weight_in(10)

type(sliced_eles_struct), allocatable :: eles(:)
type(sliced_eles_struct), allocatable :: K2eles(:)

integer i, j, w, nK2, nvar
integer, allocatable :: k2mask(:), mags(:)
integer n_slices_sext, n_slices_gen

real(rp), allocatable :: A(:, :), Ap(:, :)
real(rp), allocatable :: As(:, :), Asp(:, :)
real(rp) B(18), C(18), V(18), Weights(18)
real(rp), allocatable :: ls_soln_sliced(:)
real(rp) chrom_set_x, chrom_set_y

logical, allocatable :: mask(:)

character(8) err_str
character(17) :: r_name = 'srdt_lsq_solution'
 
!

n_slices_gen = integer_option(10, n_slices_gen_opt)
n_slices_sext = integer_option(20, n_slices_sxt_opt)
chrom_set_x = real_option(0.0d0, chrom_set_x_opt)
chrom_set_y = real_option(0.0d0, chrom_set_y_opt)

call make_slices(lat, eles, n_slices_gen, n_slices_sext, var_indexes)
w = size(eles)

!k2eles is a subset of eles containing only those slices with a valid sextupole moment.
nK2 = count(eles(:)%good_k2 .eqv. .true.)
allocate(k2eles(nK2))
k2eles(:) = pack(eles, eles(:)%good_k2)

nvar=size(var_indexes)
allocate(A(18, nvar))
allocate(Ap(nvar, 18))
allocate(ls_soln_sliced(nvar))
allocate(mask(nvar))
allocate(mags(nvar))


weights(1:2) = 1e4_rp
weights(3:) = 1.0_rp

if (present(weight_in)) then
  do i = 1, 10
    if (weight_in(i) == 0) cycle
    select case (i)
    case (1);  Weights(1) = weight_in(i)
    case (2);  Weights(2) = weight_in(i)
    case default
      Weights(2*i-3) = weight_in(i)
      Weights(2*i-2) = weight_in(i)
    end select
  enddo
endif

! A is a matrix of the coefficients of the sextupole moments that are valid variables s.t. A.K2vec = RDTs
! The columns of A are condensed such that each column corresponds to one physical sextupole. That is, 
! the contributions from each slice are added together.

do i=1, nvar
  mask = k2eles(:)%ix==var_indexes(i)
  if (count(mask) .gt. 0) then
    A(1 , i) = Weights(1)*sum(k2eles(:)%l*(-2.0)*k2eles(:)%eta_a*k2eles(:)%beta_a, mask=mask)/4.0   ! chrom_x
    A(2 , i) = Weights(2)*sum(k2eles(:)%l*(-2.0)*k2eles(:)%eta_a*k2eles(:)%beta_b, mask=mask)/(-4.0)   ! chrom_y
    A(3 , i) = Weights(3)*sum(k2eles(:)%l* real(-2.0*k2eles(:)%eta_a*k2eles(:)%beta_a * k2eles(:)%e2a), mask=mask)/8.0  !h20001
    A(4 , i) = Weights(4)*sum(k2eles(:)%l*aimag(-2.0*k2eles(:)%eta_a*k2eles(:)%beta_a * k2eles(:)%e2a), mask=mask)/8.0
    A(5 , i) = Weights(5)*sum(k2eles(:)%l* real( 2.0*k2eles(:)%eta_a*k2eles(:)%beta_b * k2eles(:)%e2b), mask=mask)/(-8.0)  !h00201
    A(6 , i) = Weights(6)*sum(k2eles(:)%l*aimag( 2.0*k2eles(:)%eta_a*k2eles(:)%beta_b * k2eles(:)%e2b), mask=mask)/(-8.0)
    A(7 , i) = Weights(7)*sum(k2eles(:)%l* real(k2eles(:)%eta_a*k2eles(:)%eta_a*sqrt(k2eles(:)%beta_a) * k2eles(:)%ea), mask=mask)/2.0 !h10002
    A(8 , i) = Weights(8)*sum(k2eles(:)%l*aimag(k2eles(:)%eta_a*k2eles(:)%eta_a*sqrt(k2eles(:)%beta_a) * k2eles(:)%ea), mask=mask)/2.0
    A(9 , i) = Weights(9)*sum(k2eles(:)%l* real(k2eles(:)%beta_a**(3./2.) * k2eles(:)%ea), mask=mask)/(-8.0)  !h21000
    A(10, i) = Weights(10)*sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(3./2.) * k2eles(:)%ea), mask=mask)/(-8.0)
    A(11, i) = Weights(11)*sum(k2eles(:)%l* real(k2eles(:)%beta_a**(3./2.) * k2eles(:)%e3a), mask=mask)/(-24.0)  !h30000
    A(12, i) = Weights(12)*sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(3./2.) * k2eles(:)%e3a), mask=mask)/(-24.0)
    A(13, i) = Weights(13)*sum(k2eles(:)%l* real(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea), mask=mask)/4.0  !h10110
    A(14, i) = Weights(14)*sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea), mask=mask)/4.0
    A(15, i) = Weights(15)*sum(k2eles(:)%l* real(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea/k2eles(:)%e2b), mask=mask)/8.0 !h10020
    A(16, i) = Weights(16)*sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea/k2eles(:)%e2b), mask=mask)/8.0
    A(17, i) = Weights(17)*sum(k2eles(:)%l* real(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea*k2eles(:)%e2b), mask=mask)/8.0 !h10200
    A(18, i) = Weights(18)*sum(k2eles(:)%l*aimag(k2eles(:)%beta_a**(1./2.)*k2eles(:)%beta_b * k2eles(:)%ea*k2eles(:)%e2b), mask=mask)/8.0
  else
    write(*, *) var_indexes(i) 
    write(err_str, '(i8)') var_indexes(i)
    call out_io (s_error$, r_name, 'element of var_indexes not identified as having valid K2', 'element ix = '//err_str, 'returning ls_soln=0')
    ls_soln=0
    return
  endif
enddo

! condense A

allocate(As(18, nvar))
allocate(Asp(nvar, 18))

do i=1, nvar
  do j=1, 18
    As(j, i) = sum(A(j, :)) !FOO /count(mask)
  enddo
enddo

call mat_pseudoinverse(As, Asp)

! C is a vector of the contributions to the DTs from those elements with valid K2 that are not variables. 
C=0.0d0
do i=1, nK2
  if (.not. any(k2eles(i)%ix==var_indexes(:))) then
    C(1 ) = C(1 ) + k2eles(i)%k2l*(-2.0)*k2eles(i)%eta_a*k2eles(i)%beta_a/4.0  ! chrom_x
    C(2 ) = C(2 ) + k2eles(i)%k2l*(-2.0)*k2eles(i)%eta_a*k2eles(i)%beta_b/(-4.0)   ! chrom_y
    C(3 ) = C(3 ) + k2eles(i)%k2l* real(-2.0*k2eles(i)%eta_a*k2eles(i)%beta_a * k2eles(i)%e2a)/8.0  !h20001
    C(4 ) = C(4 ) + k2eles(i)%k2l*aimag(-2.0*k2eles(i)%eta_a*k2eles(i)%beta_a * k2eles(i)%e2a)/8.0
    C(5 ) = C(5 ) + k2eles(i)%k2l* real( 2.0*k2eles(i)%eta_a*k2eles(i)%beta_b * k2eles(i)%e2b)/(-8.0)  !h00201
    C(6 ) = C(6 ) + k2eles(i)%k2l*aimag( 2.0*k2eles(i)%eta_a*k2eles(i)%beta_b * k2eles(i)%e2b)/(-8.0)
    C(7 ) = C(7 ) + k2eles(i)%k2l* real(k2eles(i)%eta_a*k2eles(i)%eta_a*sqrt(k2eles(i)%beta_a) * k2eles(i)%ea)/2.0 !h10002
    C(8 ) = C(8 ) + k2eles(i)%k2l*aimag(k2eles(i)%eta_a*k2eles(i)%eta_a*sqrt(k2eles(i)%beta_a) * k2eles(i)%ea)/2.0
    C(9 ) = C(9 ) + k2eles(i)%k2l* real(k2eles(i)%beta_a**(3./2.) * k2eles(i)%ea)/(-8.0)  !h21000
    C(10) = C(10) + k2eles(i)%k2l*aimag(k2eles(i)%beta_a**(3./2.) * k2eles(i)%ea)/(-8.0)
    C(11) = C(11) + k2eles(i)%k2l* real(k2eles(i)%beta_a**(3./2.) * k2eles(i)%e3a)/(-24.0)  !h30000
    C(12) = C(12) + k2eles(i)%k2l*aimag(k2eles(i)%beta_a**(3./2.) * k2eles(i)%e3a)/(-24.0)
    C(13) = C(13) + k2eles(i)%k2l* real(k2eles(i)%beta_a**(1./2.)*k2eles(i)%beta_b * k2eles(i)%ea)/4.0  !h10110
    C(14) = C(14) + k2eles(i)%k2l*aimag(k2eles(i)%beta_a**(1./2.)*k2eles(i)%beta_b * k2eles(i)%ea)/4.0
    C(15) = C(15) + k2eles(i)%k2l* real(k2eles(i)%beta_a**(1./2.)*k2eles(i)%beta_b * k2eles(i)%ea/k2eles(i)%e2b)/8.0 !h10020
    C(16) = C(16) + k2eles(i)%k2l*aimag(k2eles(i)%beta_a**(1./2.)*k2eles(i)%beta_b * k2eles(i)%ea/k2eles(i)%e2b)/8.0
    C(17) = C(17) + k2eles(i)%k2l* real(k2eles(i)%beta_a**(1./2.)*k2eles(i)%beta_b * k2eles(i)%ea*k2eles(i)%e2b)/8.0 !h10200
    C(18) = C(18) + k2eles(i)%k2l*aimag(k2eles(i)%beta_a**(1./2.)*k2eles(i)%beta_b * k2eles(i)%ea*k2eles(i)%e2b)/8.0
  endif
enddo

! B is a vector of the K1 contributions to the DTs.
B(1:18) = 0.0d0
B(1) = -chrom_set_x+sum(eles(:)%k1l*eles(:)%beta_a)/4.0
B(2) = chrom_set_y+sum(eles(:)%k1l*eles(:)%beta_b)/(-4.0)
B(3) = sum(real(eles(:)%k1l*eles(:)%beta_a*eles(:)%e2a))/8.0
B(4) = sum(aimag(eles(:)%k1l*eles(:)%beta_a*eles(:)%e2a))/8.0
B(5) = sum(real(-eles(:)%k1l*eles(:)%beta_b*eles(:)%e2b))/(-8.0)
B(6) = sum(aimag(-eles(:)%k1l*eles(:)%beta_b*eles(:)%e2b))/(-8.0)
B(7) = sum(real(-eles(:)%k1l*eles(:)%eta_a*sqrt(eles(:)%beta_a) * eles(:)%ea))/2.0
B(8) = sum(aimag(-eles(:)%k1l*eles(:)%eta_a*sqrt(eles(:)%beta_a) * eles(:)%ea))/2.0

do i=1, 18
  V(i) = Weights(i)*(-B(i)-C(i))
enddo

!ls_soln = 2.0*matmul(Ap, -B-C) 
allocate(ls_soln(nvar))
ls_soln = 2.0*matmul(Asp, V) 

end subroutine srdt_lsq_solution

end module

