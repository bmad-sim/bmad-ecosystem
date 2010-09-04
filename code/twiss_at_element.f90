!+
! Subroutine twiss_at_element (lat, ix_ele, start, end, average)
! 
! Subroutine to return the twiss parameters at the beginning, 
! end or the average of an element
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat    -- lat_struct: Lat holding the lattice.
!   ix_ele  -- Integer: Index of element
!
! Output:
!   start   -- Ele_struct, optional: Twiss and s at start of element.
!   end     -- Ele_struct, optional: Twiss and s at end of element.
!   average -- Ele_struct, optional: Average Twiss and s of element.
!     %value(l$) -- "Effective" length which for groups and overlays
!                      are weighted by the control coefficient.
!
! Warning: This routine just takes the average to be the average of both ends.
! this does not do a good job of calculating the average for some elements.
!-

recursive subroutine twiss_at_element (lat, ix_ele, start, end, average)

  
  use lat_ele_loc_mod, except_dummy => twiss_at_element
  use nr
  
  implicit none

  type (lat_struct), target :: lat
  type (ele_struct), optional :: start, end, average
  type (ele_struct), pointer :: ele
  type (ele_struct) slave_ave
  type (ele_pointer_struct), allocatable, save :: slaves(:)

  integer ix_ele, ix1, ix2
  integer i, ix, n_slave, ct

  real(rp) rr, tot, l_now

!

  ele => lat%ele(ix_ele)

! Element 0 is easy.

  if (ix_ele == 0) then
    if (present(start))   start   = ele
    if (present(end))     end     = ele
    if (present(average)) average = ele
    return
  endif        

! Elements in the tracking part of the lattice are also easy.

  if (ix_ele <= lat%n_ele_track) then
    if (present(start)) start = lat%ele(ix_ele - 1)
    if (present(end))   end = ele
    if (present(average)) then
      call zero_ave (average)
      call twiss_ave (average, lat%ele(ix_ele - 1), 0.5_rp)
      call twiss_ave (average, ele, 0.5_rp)
      average%value = ele%value
    endif
    return
  endif

! Start and end calculation for the lord elements

  select case (ele%lord_status)
  case (multipass_lord$, super_lord$, girder_lord$)
    if (present(start)) then
      ele => pointer_to_slave(lat, ele, 1)
      start = ele
    endif
    if (present(end)) then
      ele => pointer_to_slave(lat, ele, ele%n_slave)
      end = ele
    endif

  case default   ! overlay_lord$ or group_lord$

    call get_element_slave_list (lat, ele, slaves, n_slave)
    ix1 = lat%n_ele_track;  ix2 = 0
    do i = 1, n_slave
      ix1 = min (ix1, slaves(i)%ele%ix_ele)
      ix2 = max (ix2, slaves(i)%ele%ix_ele)
    enddo

    if (ix2-ix1 < lat%n_ele_track/2) then
      if (present(start)) start = lat%ele(ix1-1)
      if (present(end))   end   = lat%ele(ix2)
    else
      if (present(start)) start = lat%ele(ix2-1)
      if (present(end))   end   = lat%ele(ix1)
    endif

  end select

! Average calc for a lord.
! %control(:)%coef is proportional to the length of the save for super_lords
! The "length" is tricky for group and overlay lords: Essentially it is the 
! length weighted by the coefficient.

  if (.not. present(average)) return
  call zero_ave (average)
  if (ele%n_slave == 0) return
  ct = ele%lord_status

  tot = 0
  do i = ele%ix1_slave, ele%ix2_slave
    ix = lat%control(i)%ix_slave
    if (ct == group_lord$ .or. ct == overlay_lord$) then
      tot = tot + abs(lat%control(i)%coef) * lat%ele(ix)%value(l$)
    else
      tot = tot + lat%ele(ix)%value(l$)
    endif
  enddo
  average%value(l$) = tot

  do i = ele%ix1_slave, ele%ix2_slave
    ix = lat%control(i)%ix_slave
    call twiss_at_element (lat, ix, average = slave_ave)
    if (tot == 0) then
      rr = 1.0 / ele%n_slave
    elseif (ct == group_lord$ .or. ct == overlay_lord$) then
      rr = abs(lat%control(i)%coef) * slave_ave%value(l$) / tot
    else
      rr = slave_ave%value(l$) / tot
    endif
    call twiss_ave (average, slave_ave, rr)
  enddo

!--------------------------------------------------------------------------
contains

subroutine zero_ave (ave)

  type (ele_struct) ave

  ave%s          = 0
  ave%a%phi      = 0;   ave%b%phi      = 0
  ave%a%alpha    = 0;   ave%b%alpha    = 0
  ave%a%beta     = 0;   ave%b%beta     = 0
  ave%a%gamma    = 0;   ave%b%gamma    = 0
  ave%a%eta      = 0;   ave%b%eta      = 0
  ave%a%etap     = 0;   ave%b%etap     = 0
  ave%a%sigma    = 0;   ave%b%sigma    = 0
  ave%x%eta  = 0;   ave%y%eta  = 0
  ave%x%etap = 0;   ave%y%etap = 0
  ave%c_mat      = 0
  ave%gamma_c    = 0                                            
  ave%value(l$)  = 0

end subroutine

!--------------------------------------------------------------------------

subroutine twiss_ave (ave, e1, r)

  type (ele_struct) ave, e1

  real(rp) r

!

  ave%s          = ave%s          + r * e1%s          
  ave%c_mat      = ave%c_mat      + r * e1%c_mat      
  ave%gamma_c    = ave%gamma_c    + r * e1%gamma_c    
  ave%a%phi      = ave%a%phi      + r * e1%a%phi      
  ave%a%alpha    = ave%a%alpha    + r * e1%a%alpha    
  ave%a%beta     = ave%a%beta     + r * e1%a%beta     
  ave%a%gamma    = ave%a%gamma    + r * e1%a%gamma    
  ave%a%eta      = ave%a%eta      + r * e1%a%eta      
  ave%a%etap     = ave%a%etap     + r * e1%a%etap     
  ave%a%sigma    = ave%a%sigma    + r * e1%a%sigma    
  ave%x%eta  = ave%x%eta  + r * e1%x%eta  
  ave%x%etap = ave%x%etap + r * e1%x%etap 
  ave%b%phi      = ave%b%phi      + r * e1%b%phi      
  ave%b%alpha    = ave%b%alpha    + r * e1%b%alpha    
  ave%b%beta     = ave%b%beta     + r * e1%b%beta     
  ave%b%gamma    = ave%b%gamma    + r * e1%b%gamma    
  ave%b%eta      = ave%b%eta      + r * e1%b%eta      
  ave%b%etap     = ave%b%etap     + r * e1%b%etap     
  ave%b%sigma    = ave%b%sigma    + r * e1%b%sigma    
  ave%y%eta  = ave%y%eta  + r * e1%y%eta  
  ave%y%etap = ave%y%etap + r * e1%y%etap 

end subroutine

end subroutine
