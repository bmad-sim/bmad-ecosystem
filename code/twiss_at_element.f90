!+
! Subroutine twiss_at_element (ring, ix_ele, start, end, average)
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
! Warning: This routine just takes the average to be the average of both ends.
! this does not do a good job of calculating the average for some elements.
!-

#include "CESR_platform.inc"

subroutine twiss_at_element (ring, ix_ele, start, end, average)

  use bmad_struct
  use bmad_interface, except => twiss_at_element
  use nr
  
  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), optional :: start, end, average

  integer ix_ele, ix1, ix2
  integer i, n_slave
  integer, allocatable :: slave_list(:)

  real(rp) rr, l_tot, l_now

! Element 0 is easy.

  if (ix_ele == 0) then
    if (present(start))   start   = ring%ele_(0)
    if (present(end))     end     = ring%ele_(0)
    if (present(average)) average = ring%ele_(0)
    return
  endif        

! Regular elements are also easy.

  if (ix_ele <= ring%n_ele_use) then
    if (present(start)) start = ring%ele_(ix_ele - 1)
    if (present(end))   end = ring%ele_(ix_ele)
    if (present(average)) then
      call zero_ave (average)
      call twiss_ave (average, ring%ele_(ix_ele - 1), &
                                             ring%ele_(ix_ele), 0.5_rp)
      average%value = ring%ele_(ix_ele)%value
    endif
    return
  endif

! calculation for the lords...

  call get_element_slave_list (ring, ix_ele, slave_list, n_slave)

  if (ring%ele_(ix_ele)%control_type == super_lord$) then
    ix1 = ring%ele_(ix_ele)%ix1_slave
    ix2 = ring%ele_(ix_ele)%ix2_slave
    if (present(start)) start = ring%ele_(ring%control_(ix1)%ix_slave - 1)
    if (present(end)) end = ring%ele_(ring%control_(ix2)%ix_slave)

  else  ! overlay_lord$ or group_lord$
    ix1 = minval (slave_list(1:n_slave))
    ix2 = maxval (slave_list(1:n_slave))

    if (ix2-ix1 < ring%n_ele_use/2) then
      if (present(start)) start = ring%ele_(ix1-1)
      if (present(end))   end   = ring%ele_(ix2)
    else
      if (present(start)) start = ring%ele_(ix2-1)
      if (present(end))   end   = ring%ele_(ix1)
    endif

  endif

! Average calc for a lord.

  if (.not. present(average)) return
  call zero_ave (average)

  l_tot = 0
  do i = 1, n_slave
    l_tot = l_tot + ring%ele_(slave_list(i))%value(l$)
  enddo
  average%value(l$) = l_tot 

  do i = 1, n_slave
    if (l_tot == 0) then
      rr = 1.0 / n_slave
    else
      rr = ring%ele_(i)%value(l$) / l_tot
    endif
    call twiss_ave (average, ring%ele_(i - 1), ring%ele_(i), rr)
  enddo

!--------------------------------------------------------------------------
contains

subroutine zero_ave (ave)

  type (ele_struct) ave

  ave%s          = 0
  ave%x%phi      = 0;   ave%y%phi      = 0
  ave%x%alpha    = 0;   ave%y%alpha    = 0
  ave%x%beta     = 0;   ave%y%beta     = 0
  ave%x%gamma    = 0;   ave%y%gamma    = 0
  ave%x%eta      = 0;   ave%y%eta      = 0
  ave%x%etap     = 0;   ave%y%etap     = 0
  ave%x%sigma    = 0;   ave%y%sigma    = 0
  ave%x%eta_lab  = 0;   ave%y%eta_lab  = 0
  ave%x%etap_lab = 0;   ave%y%etap_lab = 0
  ave%c_mat      = 0
  ave%gamma_c    = 0                                            
  ave%value(l$)  = 0

end subroutine

!--------------------------------------------------------------------------

subroutine twiss_ave (ave, e1, e2, r)

  type (ele_struct) ave, e1, e2

  real(rp) r

!

  ave%s          = ave%s          + r * e1%s          + r * e2%s
  ave%c_mat      = ave%c_mat      + r * e1%c_mat      + r * e2%c_mat
  ave%gamma_c    = ave%gamma_c    + r * e1%gamma_c    + r * e2%gamma_c
  ave%x%phi      = ave%x%phi      + r * e1%x%phi      + r * e2%x%phi
  ave%x%alpha    = ave%x%alpha    + r * e1%x%alpha    + r * e2%x%alpha
  ave%x%beta     = ave%x%beta     + r * e1%x%beta     + r * e2%x%beta
  ave%x%gamma    = ave%x%gamma    + r * e1%x%gamma    + r * e2%x%gamma
  ave%x%eta      = ave%x%eta      + r * e1%x%eta      + r * e2%x%eta
  ave%x%etap     = ave%x%etap     + r * e1%x%etap     + r * e2%x%etap
  ave%x%sigma    = ave%x%sigma    + r * e1%x%sigma    + r * e2%x%sigma
  ave%x%eta_lab  = ave%x%eta_lab  + r * e1%x%eta_lab  + r * e2%x%eta_lab
  ave%x%etap_lab = ave%x%etap_lab + r * e1%x%etap_lab + r * e2%x%etap_lab
  ave%y%phi      = ave%y%phi      + r * e1%y%phi      + r * e2%y%phi
  ave%y%alpha    = ave%y%alpha    + r * e1%y%alpha    + r * e2%y%alpha
  ave%y%beta     = ave%y%beta     + r * e1%y%beta     + r * e2%y%beta
  ave%y%gamma    = ave%y%gamma    + r * e1%y%gamma    + r * e2%y%gamma
  ave%y%eta      = ave%y%eta      + r * e1%y%eta      + r * e2%y%eta
  ave%y%etap     = ave%y%etap     + r * e1%y%etap     + r * e2%y%etap
  ave%y%sigma    = ave%y%sigma    + r * e1%y%sigma    + r * e2%y%sigma
  ave%y%eta_lab  = ave%y%eta_lab  + r * e1%y%eta_lab  + r * e2%y%eta_lab
  ave%y%etap_lab = ave%y%etap_lab + r * e1%y%etap_lab + r * e2%y%etap_lab

end subroutine

end subroutine
