!+
! Subroutine TWISS_AT_ELEMENT (ring, ix_ele, start, end, average)
! 
! Subroutine to return the twiss parameters at the beginning, 
! end or the average of an element
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring    -- Ring_struct: Ring holding the lattice.
!   ix_ele  -- Integer: Index of element
!
! Output:
!   start   -- Ele_struct: [Optional] Twiss and s at start of element.
!   end     -- Ele_struct: [Optional] Twiss and s at end of element.
!   average -- Ele_struct: [Optional] Average Twiss and s of element.
!
! Note: This Subroutine uses INDEXX from Numerical Recipes (F90 edition)
!
! Warning: This routine just takes the average to be the average of both ends.
! this does not do a good job of calculating the average for some elements.
!-

!$Id$
!$Log$
!Revision 1.5  2003/05/02 15:44:03  dcs
!F90 standard conforming changes.
!
!Revision 1.4  2003/01/27 14:40:46  dcs
!bmad_version = 56
!
!Revision 1.3  2002/02/23 20:32:28  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:32:00  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine twiss_at_element (ring, ix_ele, start, end, average)

  use bmad_struct
  use bmad_interface
  use nr
  
  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), optional :: start, end, average
  type (ele_struct), pointer :: s_ele, e_ele, ele

  integer ix_ele, ix1, ix2, n, ix_(100), indx(100)
  integer i, j, ix1_2, ix2_2, i_2, j_2

  real(rdef) rr, rr_2, coef_tot, coef_tot_2

! start and end

  if (ix_ele == 0) then
    print *, 'ERROR IN TWISS_AT_ELEMENT: IX_ELE = 0'
    call err_exit
  endif        

  ix1 = ring%ele_(ix_ele)%ix1_slave
  ix2 = ring%ele_(ix_ele)%ix2_slave

  if (ix_ele <= ring%n_ele_ring) then  ! in regular part of the ring.
    if (present(start)) start = ring%ele_(ix_ele - 1)
    if (present(end)) end = ring%ele_(ix_ele)
  elseif (ring%ele_(ix_ele)%control_type == super_lord$) then
    if (present(start)) start = ring%ele_(ring%control_(ix1)%ix_slave - 1)
    if (present(end)) end = ring%ele_(ring%control_(ix2)%ix_slave)
  else  ! overlay_lord$ or group_lord$
    if (present(start) .or. present(end)) then
      n = ring%ele_(ix_ele)%n_slave
      ix_(1:n) = ring%control_(ix1:ix2)%ix_slave
      ix_(2:n) = ix_(2:n) + ring%n_ele_ring * &
                             nint(float(ix_(1) - ix_(2:n)) / ring%n_ele_ring)
      call indexx (ix_(1:n), indx(1:n))
      if (present(start)) start = ring%ele_(ix_(indx(1)) - 1)
      if (present(end)) end = ring%ele_(ix_(indx(n)))
    endif
  endif

! average

  if (.not. present(average)) return

  call zero_ave (average)

  if (ix_ele <= ring%n_ele_ring) then
    call twiss_ave (average, ring%ele_(ix_ele - 1), &
                                             ring%ele_(ix_ele), 0.5_rdef)
    average%value = ring%ele_(ix_ele)%value

  else

    coef_tot = sum (ring%control_(ix1:ix2)%coef)

    do i = ix1, ix2
      j = ring%control_(i)%ix_slave
      average%value(l$) = average%value(l$) + ring%ele_(j)%value(l$)
      rr = ring%control_(i)%coef / (2 * coef_tot)
      if (ring%ele_(j)%n_slave == 0) then
        call twiss_ave (average, ring%ele_(j - 1), ring%ele_(j), rr)
      else
        ix1_2 = ring%ele_(j)%ix1_slave
        ix2_2 = ring%ele_(j)%ix2_slave
        coef_tot_2 = sum (ring%control_(ix1_2:ix2_2)%coef)
        do i_2 = ix1_2, ix2_2
          j_2 = ring%control_(i_2)%ix_slave
          rr_2 = rr * ring%control_(i_2)%coef / (2 * coef_tot_2)
          call twiss_ave (average, ring%ele_(j_2 - 1), ring%ele_(j_2), rr_2)
        enddo
      endif
    enddo
  endif

!--------------------------------------------------------------------------
contains

subroutine zero_ave (ave)

  type (ele_struct) ave

  ave%s         = 0
  ave%x%phi     = 0;   ave%y%phi     = 0
  ave%x%alpha   = 0;   ave%y%alpha   = 0
  ave%x%beta    = 0;   ave%y%beta    = 0
  ave%x%gamma   = 0;   ave%y%gamma   = 0
  ave%x%eta     = 0;   ave%y%eta     = 0
  ave%x%etap    = 0;   ave%y%etap    = 0
  ave%x%sigma   = 0;   ave%y%sigma   = 0
  ave%x%mobius_beta = 0;   ave%y%mobius_beta = 0
  ave%x%mobius_eta  = 0;   ave%y%mobius_eta  = 0
  ave%c_mat     = 0
  ave%gamma_c   = 0                                            
  ave%value(l$) = 0

end subroutine

!--------------------------------------------------------------------------

subroutine twiss_ave (ave, e1, e2, r)

  type (ele_struct) ave, e1, e2

  real(rdef) r

!

  ave%s         = ave%s         + r * e1%s         + r * e2%s
  ave%c_mat     = ave%c_mat     + r * e1%c_mat     + r * e2%c_mat
  ave%gamma_c   = ave%gamma_c   + r * e1%gamma_c   + r * e2%gamma_c
  ave%x%phi     = ave%x%phi     + r * e1%x%phi     + r * e2%x%phi
  ave%x%alpha   = ave%x%alpha   + r * e1%x%alpha   + r * e2%x%alpha
  ave%x%beta    = ave%x%beta    + r * e1%x%beta    + r * e2%x%beta
  ave%x%gamma   = ave%x%gamma   + r * e1%x%gamma   + r * e2%x%gamma
  ave%x%eta     = ave%x%eta     + r * e1%x%eta     + r * e2%x%eta
  ave%x%etap    = ave%x%etap    + r * e1%x%etap    + r * e2%x%etap
  ave%x%sigma   = ave%x%sigma   + r * e1%x%sigma   + r * e2%x%sigma
  ave%x%mobius_beta = ave%x%mobius_beta + r * e1%x%mobius_beta + &
                                                      r * e2%x%mobius_beta
  ave%x%mobius_eta  = ave%x%mobius_eta  + r * e1%x%mobius_eta  + &
                                                      r * e2%x%mobius_eta
  ave%y%phi     = ave%y%phi     + r * e1%y%phi     + r * e2%y%phi
  ave%y%alpha   = ave%y%alpha   + r * e1%y%alpha   + r * e2%y%alpha
  ave%y%beta    = ave%y%beta    + r * e1%y%beta    + r * e2%y%beta
  ave%y%gamma   = ave%y%gamma   + r * e1%y%gamma   + r * e2%y%gamma
  ave%y%eta     = ave%y%eta     + r * e1%y%eta     + r * e2%y%eta
  ave%y%etap    = ave%y%etap    + r * e1%y%etap    + r * e2%y%etap
  ave%y%sigma   = ave%y%sigma   + r * e1%y%sigma   + r * e2%y%sigma
  ave%y%mobius_beta = ave%y%mobius_beta + r * e1%y%mobius_beta + &
                                                      r * e2%y%mobius_beta
  ave%y%mobius_eta  = ave%y%mobius_eta  + r * e1%y%mobius_eta  + &
                                                      r * e2%y%mobius_eta

end subroutine

end subroutine
