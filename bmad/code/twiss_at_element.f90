!+
! Subroutine twiss_at_element (ele, start, end, average)
! 
! Subroutine to return the twiss parameters at the beginning, 
! end or the average of an element.
!
! Note: Ele must be an element in a lattice.
!
! Warning: This routine just takes the average to be the average of both ends of the element.
! this does not do a good job of calculating the average for some elements.
!
! Input:
!   ele    -- ele_struct: Element to be averaged
!
! Output:
!   start   -- Ele_struct, optional: Twiss and s at start of element.
!   end     -- Ele_struct, optional: Twiss and s at end of element.
!   average -- Ele_struct, optional: Average Twiss and s of element.
!     %value(l$) -- "Effective" length which for groups and overlays
!                      are weighted by the control coefficient.
!-

recursive subroutine twiss_at_element (ele, start, end, average)

use expression_mod
use equal_mod, dummy => twiss_at_element

implicit none

type (branch_struct), pointer :: branch
type (ele_struct), optional :: start, end, average
type (ele_struct), target :: ele
type (ele_struct), pointer :: slave, ele2
type (ele_struct) slave_ave
type (ele_pointer_struct), allocatable :: slaves(:)
type (control_struct), pointer :: ctl

integer ix_ele, ix1, ix2
integer i, ix, n_slave, key

real(rp) rr, tot, l_now

logical err_flag

! Element 0 is easy.

branch => ele%branch
ix_ele = ele%ix_ele

if (ix_ele == 0) then
  if (present(start))   start   = ele
  if (present(end))     end     = ele
  if (present(average)) average = ele
  return
endif        

! Elements in the tracking part of the lattice are also easy.

if (ix_ele <= branch%n_ele_track) then
  if (present(start)) start = branch%ele(ix_ele - 1)
  if (present(end))   end = ele
  if (present(average)) then
    call zero_ave (average)
    call twiss_ave (average, branch%ele(ix_ele - 1), 0.5_rp)
    call twiss_ave (average, ele, 0.5_rp)
    average%value = ele%value
  endif
  return
endif

! Start and end calculation for the lord elements

select case (ele%lord_status)
case (multipass_lord$, super_lord$, girder_lord$)
  if (present(start)) then
    ele2 => pointer_to_slave(ele, 1)
    start = ele2
  endif
  if (present(end)) then
    ele2 => pointer_to_slave(ele, ele%n_slave)
    end = ele2
  endif

case default   ! overlay$ or group$

  call get_slave_list (ele, slaves, n_slave)
  ix1 = branch%n_ele_track;  ix2 = 0
  do i = 1, n_slave
    ix1 = min (ix1, slaves(i)%ele%ix_ele)
    ix2 = max (ix2, slaves(i)%ele%ix_ele)
  enddo

  if (ix2-ix1 < branch%n_ele_track/2) then
    if (present(start)) start = branch%ele(ix1-1)
    if (present(end))   end   = branch%ele(ix2)
  else
    if (present(start)) start = branch%ele(ix2-1)
    if (present(end))   end   = branch%ele(ix1)
  endif

end select

! Average calc for a lord.
! We need the "length" of the element.
! The "length" is tricky for group and overlay lords: Essentially it is the 
! length weighted by the coefficient.

if (.not. present(average)) return
call zero_ave (average)
if (ele%n_slave == 0) return
key = ele%key

tot = 0
do i = 1, ele%n_slave
  slave => pointer_to_slave(ele, i, ctl)
  if (key == group$ .or. key == overlay$) then
    tot = tot + abs(linear_coef(ctl%stack, err_flag)) * slave%value(l$)
  else
    tot = tot + slave%value(l$)
  endif
enddo

do i = 1, ele%n_slave
  slave => pointer_to_slave(ele, i, ctl)
  call twiss_at_element (slave, average = slave_ave)
  if (tot == 0) then
    rr = 1.0 / ele%n_slave
  elseif (key == group$ .or. key == overlay$) then
    rr = abs(linear_coef(ctl%stack, err_flag)) * slave_ave%value(l$) / tot
  else
    rr = slave_ave%value(l$) / tot
  endif
  call twiss_ave (average, slave_ave, rr)
enddo

!--------------------------------------------------------------------------
contains

subroutine zero_ave (ave)

type (ele_struct) ave

ave%s         = 0;   ave%s_start   = 0
ave%a%phi     = 0;   ave%b%phi     = 0
ave%a%alpha   = 0;   ave%b%alpha   = 0
ave%a%beta    = 0;   ave%b%beta    = 0
ave%a%gamma   = 0;   ave%b%gamma   = 0
ave%a%eta     = 0;   ave%b%eta     = 0
ave%a%etap    = 0;   ave%b%etap    = 0
ave%a%deta_ds = 0;   ave%b%deta_ds = 0
ave%a%sigma   = 0;   ave%b%sigma   = 0
ave%x%eta     = 0;   ave%y%eta     = 0
ave%x%etap    = 0;   ave%y%etap    = 0
ave%x%deta_ds = 0;   ave%y%deta_ds = 0
ave%c_mat     = 0
ave%gamma_c   = 0                                            
ave%value(l$) = 0

end subroutine

!--------------------------------------------------------------------------

subroutine twiss_ave (ave, ele1, r)

type (ele_struct) ave, ele1

real(rp) r

!

ave%s         = ave%s         + r * ele1%s          
ave%s_start   = ave%s_start   + r * ele1%s_start
ave%c_mat     = ave%c_mat     + r * ele1%c_mat      
ave%gamma_c   = ave%gamma_c   + r * ele1%gamma_c    
ave%a%phi     = ave%a%phi     + r * ele1%a%phi      
ave%a%alpha   = ave%a%alpha   + r * ele1%a%alpha    
ave%a%beta    = ave%a%beta    + r * ele1%a%beta     
ave%a%gamma   = ave%a%gamma   + r * ele1%a%gamma    
ave%a%eta     = ave%a%eta     + r * ele1%a%eta      
ave%a%etap    = ave%a%etap    + r * ele1%a%etap     
ave%a%deta_ds = ave%a%deta_ds + r * ele1%a%deta_ds     
ave%a%sigma   = ave%a%sigma   + r * ele1%a%sigma    
ave%x%eta     = ave%x%eta     + r * ele1%x%eta  
ave%x%etap    = ave%x%etap    + r * ele1%x%etap 
ave%x%deta_ds = ave%x%deta_ds + r * ele1%x%deta_ds 
ave%b%phi     = ave%b%phi     + r * ele1%b%phi      
ave%b%alpha   = ave%b%alpha   + r * ele1%b%alpha    
ave%b%beta    = ave%b%beta    + r * ele1%b%beta     
ave%b%gamma   = ave%b%gamma   + r * ele1%b%gamma    
ave%b%eta     = ave%b%eta     + r * ele1%b%eta      
ave%b%etap    = ave%b%etap    + r * ele1%b%etap     
ave%b%deta_ds = ave%b%deta_ds + r * ele1%b%deta_ds     
ave%b%sigma   = ave%b%sigma   + r * ele1%b%sigma    
ave%y%eta     = ave%y%eta     + r * ele1%y%eta  
ave%y%etap    = ave%y%etap    + r * ele1%y%etap 
ave%y%deta_ds = ave%y%deta_ds + r * ele1%y%deta_ds 
ave%value(l$) = ave%value(l$) + ele1%value(l$) 

end subroutine

end subroutine
