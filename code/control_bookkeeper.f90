!+
! Subroutine control_bookkeeper (ring, ix_ele)
!
! Subroutine to transfer attibute information from lord to slave elements.
! It is assumend that the element with index ix_ele is the only element whose
! attribute values have been changed.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   RING   -- Ring_struct: Ring to be used
!   IX_ELE -- Integer: Index of element whose attribute values have been
!               changed.
!-

subroutine control_bookkeeper (ring, ix_ele)

  use bmad_struct

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  integer ix_ele, i, j, ix, ix1, ix2
  integer ix_eles(300)

! Make a list of elements to update.
! we do not need to update free elements of group lords.

  ix1 = 0   ! index for processed elements
  ix_eles(1) = ix_ele
  ix2 = 1   ! index for last element in list

  do
    ix1 = ix1 + 1
    ele => ring%ele_(ix_eles(ix1))
    do i = ele%ix1_slave, ele%ix2_slave
      ix = ring%control_(i)%ix_slave
      if (ele%key == group_lord$ .and. ring%ele_(ix)%key == free$) cycle
      if (ix == ix_eles(ix2)) cycle   ! do not use duplicates
      ix2 = ix2 + 1
      ix_eles(ix2) = ix
    enddo
    if (ix1 == ix2) exit
  enddo

! First: Makup Lords

  do j = 1, ix2

    ix = ix_eles(j)
    ele => ring%ele_(ix)

    if (ele%control_type == group_lord$) then
      call makeup_group_slaves (ring, ix)

    elseif (ele%control_type == super_lord$ .and. ele%n_lord > 0) then
      call makeup_overlay_slave (ring, ix)

    elseif (ele%control_type == overlay_lord$ .and. ele%n_lord > 0) then
      call makeup_overlay_slave (ring, ix)

    endif

  enddo

! Second: Makeup Slaves

  do j = 1, ix2

    ix = ix_eles(j)
    ele => ring%ele_(ix)       

    if (ele%control_type == container_slave$) then
      call makeup_container_slave (ring, ix)
      
    elseif (ele%control_type == super_slave$) then
      call makeup_super_slave (ring, ix)

    elseif (ele%control_type == overlay_slave$) then
      call makeup_overlay_slave (ring, ix)

    endif

  enddo

end subroutine
                               
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_group_slaves (ring, ix_slave)
!
! Subroutine to calculate the attributes of group slave elements
!-

Subroutine makeup_group_slaves (ring, ix_lord)   

  use bmad_struct
  use dcslib_interface

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave

  real delta, coef

  integer ix_lord, ix, iv, ict

  logical moved

!

  lord => ring%ele_(ix_lord)

  delta = lord%value(command$) - lord%value(old_command$)    ! change
  lord%value(old_command$) = lord%value(command$) ! save old

  moved = .false.   ! have we longitudinally moved an element?

  do i = lord%ix1_slave, lord%ix2_slave

    ix = ring%control_(i)%ix_slave
    iv = ring%control_(i)%ix_attrib
    ict = ring%ele_(ix)%control_type
    slave => ring%ele_(ix)

    if (iv == l$) then
      moved = .true.
      if (ict /= free$ .and. ict /= super_slave$) then
        type *, 'ERROR IN CONTROL_BOOKKEEPER: A GROUP: ', lord%name
        type *, '      CONTROLS THE LENGTH OF A LORD ELEMENT: ', slave%name
        call err_exit
      endif
    endif
    coef = ring%control_(i)%coef
    slave%value(iv) = slave%value(iv) + delta * coef
  enddo

  if (moved) call s_calc (ring)       ! recompute s distances


end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine MAKEUP_SUPER_SLAVE (RING, IX_SLAVE)
!
! Subroutine to calcualte the attributes of overlay slave elements
!-
         
subroutine makeup_super_slave (ring, ix_slave)

  use bmad_struct
  use dcslib_interface

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave

  integer i, j, ix, ix_slave

  real tilt, k1_x, k1_y, x_kick, y_kick, ks, k1, coef
  real x_offset, y_offset, x_pitch, y_pitch, s_slave, s_del
  real sin_1, cos_1, x_off, y_off, a(0:n_pole_maxx), b(0:n_pole_maxx)
  real knl(0:n_pole_maxx), t(0:n_pole_maxx), value(n_attrib_maxx)
  real a_tot(0:n_pole_maxx), b_tot(0:n_pole_maxx)
                                                      
! Super_slave:

  slave => ring%ele_(ix_slave)
                     
  if (slave%control_type /= super_slave$) then
    type *, 'ERROR IN MAKEUP_SUPER_SLAVE: ELEMENT IS NOT AN SUPER SLAVE:'
    type *, '      ', slave%name
    call err_exit
  endif

! If this slave is the last slave for some lord (so that then longitudinal 
! end of the slave matches the end of the lord) then the limits of the lord
! are transfered to the slave.

!-----------------------------------------------------------------------
! 1 controller: just transfer attributes except length

  if (slave%n_lord == 1) then

    i = ring%ic_(slave%ic1_lord)  
    ix = ring%control_(i)%ix_lord
    lord => ring%ele_(ix)
    coef = ring%control_(i)%coef

    value = lord%value
    value(l$) = slave%value(l$)                 ! do not change slave length
    value(hkick$) = lord%value(hkick$) * coef
    value(vkick$) = lord%value(vkick$) * coef
    if (slave%key == rfcavity$) value(volt$) = lord%value(volt$) * coef
    if (i /= lord%ix2_slave) then   ! if not at end of the lord domain
      value(x_limit$) = 0
      value(y_limit$) = 0
    endif

    if (value(x_pitch$) /= 0 .or. value(y_pitch$) /= 0) then
      s_del = modulo2 ((slave%s - slave%value(l$)/2) - &
                  (lord%s - lord%value(l$)/2), ring%param%total_length/2)
      value(x_offset$) = value(x_offset$) + s_del * value(x_pitch$)
      value(y_offset$) = value(y_offset$) + s_del * value(y_pitch$)
    endif
      
    slave%value = value
    slave%is_on = lord%is_on
    return

  endif

!-----------------------------------------------------------------------
! Multiple controllers: must be a solenoid/quadrupole combo
! combine the lord elements.
                                           
  ks = 0
  k1_x = 0
  k1_y = 0
  x_kick = 0
  y_kick = 0
  x_offset = 0
  y_offset = 0
  x_pitch = 0
  y_pitch = 0
  a_tot = 0
  b_tot = 0

  do j = 1, n_attrib_maxx
    if (j /= l$) slave%value(j) = 0
  enddo

  s_slave = slave%s - slave%value(l$)/2  ! center of slave

! find the pitch of the slave.
! this is the pitch of the first non-solenoid lord element (if there is one).

  do j = slave%ic1_lord, slave%ic2_lord
    i = ring%ic_(j)
    ix = ring%control_(i)%ix_lord
    lord => ring%ele_(ix)
    x_pitch = lord%value(x_pitch$)
    y_pitch = lord%value(y_pitch$)
    if (lord%key /= solenoid$) exit
  enddo

! sum over all lords...
! all lords must have the same pitch except for a lord that is a solenoid.
! If a solenoid has a different pitch then it looks like a solenoid with a
! kick

  do j = slave%ic1_lord, slave%ic2_lord
    i = ring%ic_(j)
    ix = ring%control_(i)%ix_lord
    coef = ring%control_(i)%coef
    lord => ring%ele_(ix)

    if (lord%control_type /= super_lord$) then
      type *, 'ERROR IN MAKEUP_SUPER_SLAVE: SUPER_SLAVE HAS A'
      type *, '      CONTROL ELEMENT THAT IS NOT SUPER_LORD'
      type *, '      SLAVE: ', slave%name, ix_slave
      type *, '      LORD:  ', lord%name, ix
      call err_exit
    endif

    if (lord%value(x_pitch$) /= x_pitch .or. &
                                      lord%value(y_pitch$) /= y_pitch) then
      if (lord%key == solenoid$) then
        x_kick = x_kick - lord%value(ks$) * &
                           (lord%value(y_pitch$) - y_pitch) * slave%value(l$)
        y_kick = y_kick + lord%value(ks$) * &
                           (lord%value(x_pitch$) - x_pitch) * slave%value(l$)
      else
        type *, 'WARNING FROM CONTROL_BOOKKEEPER: I DO NOT YET KNOW HOW'
        type *, '        TO COMBINE DIFFERENT PITCH''S FROM 2 LORD ELEMENTS!'
        type *, '        SLAVE: ', slave%name
        call err_exit
      endif
    endif

    if (i == lord%ix2_slave) then      ! longitudinal ends match
      slave%value(x_limit$) = lord%value(x_limit$)
      slave%value(y_limit$) = lord%value(y_limit$)
    endif

! Offsets are computed so that the quadrupole field gives the correct kick.

    if (lord%is_on) then
      ks = ks + lord%value(ks$)
      x_kick = x_kick + lord%value(hkick$) * coef
      y_kick = y_kick + lord%value(vkick$) * coef
      tilt = lord%value(tilt$)
      cos_1 = lord%value(k1$) * cos(2 * tilt)
      sin_1 = lord%value(k1$) * sin(2 * tilt)
      k1_x = k1_x + cos_1
      k1_y = k1_y + sin_1
      s_del = modulo2 (s_slave - (lord%s - lord%value(l$)/2), &
                                                 ring%param%total_length/2)
      x_off = lord%value(x_offset$) + s_del * x_pitch
      y_off = lord%value(y_offset$) + s_del * y_pitch
      x_offset = x_offset + cos_1 * x_off + sin_1 * y_off
      y_offset = y_offset + sin_1 * x_off - cos_1 * y_off
      if (any(lord%value(ix1_m$:ix2_m$) /= 0)) then
        call multipole_to_vecs (lord, +1, knl, t)
        call multipole_kt_to_ab (knl/lord%value(l$), t, a, b)
        a_tot = a_tot + a
        b_tot = b_tot + b
      endif
    endif
  enddo

! stuff sums into slave element

  slave%value(ks$) = ks
  slave%value(hkick$) = x_kick
  slave%value(vkick$) = y_kick
  slave%value(x_pitch$) = x_pitch
  slave%value(y_pitch$) = y_pitch
  if (k1_x == 0 .and. k1_y == 0) then
    slave%value(k1$) = 0
    slave%value(tilt$) = 0
    slave%value(ix1_m$:ix2_m$) = 0
  else
    k1 = sqrt(k1_x**2 + k1_y**2)
    tilt = atan2(k1_y, k1_x) / 2
    if (tilt > pi/4) then
      k1 = -k1
      tilt = tilt - pi/2
    elseif (tilt < -pi/4) then
      k1 = -k1
      tilt = tilt + pi/2
    endif
    slave%value(k1$) = k1
    slave%value(tilt$) = tilt
    cos_1 = k1_x / (k1_x**2 + k1_y**2)
    sin_1 = k1_y / (k1_x**2 + k1_y**2)
    slave%value(x_offset$) =  cos_1 * x_offset + sin_1 * y_offset
    slave%value(y_offset$) =  sin_1 * x_offset - cos_1 * y_offset
    if (any(a_tot /= 0) .or. any(b_tot /= 0)) then
      call multipole_ab_to_kt(a_tot, b_tot, knl, t)
      call multipole_kt_to_ab(knl/k1, t-tilt, a, b)
      slave%value(ix1_m$:ix2_m$-1:2) = a
      slave%value(ix1_m$+1:ix2_m$:2) = b
      slave%value(radius$) = 1
    endif
  endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine makeup_overlay_slave (ring, ix_ele)

  use bmad_struct
  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  real value(n_attrib_maxx), coef
  integer i, j, ix, iv, ix_ele, icom, ct
  logical used(n_attrib_maxx)

!
                               
  ele => ring%ele_(ix_ele)
  ct = ele%control_type

  if (ct /= super_lord$ .and. ct /= overlay_slave$ .and. &
                                               ct /= overlay_lord$) then
    type *, 'ERROR IN MAKEUP_OVERLAY_SLAVE: ELEMENT IS NOT OF PROPER TYPE.'
    type *, '      RING INDEX:', ix_ele
    call type_ele (ele, .true., 0, .false., .true., ring)
    call err_exit
  endif

  value = 0
  used = .false.

  do i = ele%ic1_lord, ele%ic2_lord
    j = ring%ic_(i)
    ix = ring%control_(j)%ix_lord
    if (ring%ele_(ix)%control_type /= overlay_lord$) then
      type *, 'ERROR IN MAKEUP_OVERLAY_SLAVE:',  &
                          ' THE LORD IS NOT AN OVERLAY_LORD', ix_ele
      call type_ele (ele, .true., 0, .false., .true., ring)
      call err_exit
    endif     
    coef = ring%control_(j)%coef
    iv = ring%control_(j)%ix_attrib
    icom = ring%ele_(ix)%ix_value
    value(iv) = value(iv) + ring%ele_(ix)%value(icom)*coef
    used(iv) = .true.
  enddo

  do i = 1, n_attrib_maxx
    if (used(i)) ele%value(i) = value(i)
  enddo
                                            
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine MAKEUP_CONTAINER_SLAVE (RING, IX_ELE)
!
! Subroutine to make up a container element.
!-

subroutine makeup_container_slave (ring, ix_ele)
              
  use bmad_struct

  implicit none

  type (ring_struct)  ring

  integer i, ix_ele, j, ix
  real s_len, x_lim, y_lim

  s_len = 0
  x_lim = 0
  y_lim = 0

  do j = ring%ele_(ix_ele)%ic1_lord, ring%ele_(ix_ele)%ic2_lord
    i = ring%ic_(j)
    ix = ring%control_(i)%ix_lord
    s_len = s_len + ring%ele_(ix)%value(l$)
    if (x_lim * ring%ele_(ix)%value(x_limit$) == 0.0) then
      x_lim = max(x_lim, ring%ele_(ix)%value(x_limit$))
    else
      x_lim = min(x_lim, ring%ele_(ix)%value(x_limit$))
    endif
    if (y_lim * ring%ele_(ix)%value(y_limit$) == 0.0) then
      y_lim = max(y_lim, ring%ele_(ix)%value(y_limit$))
    else
      y_lim = min(y_lim, ring%ele_(ix)%value(y_limit$))
    endif
  enddo

  ring%ele_(ix_ele)%value(x_limit$) = x_lim
  ring%ele_(ix_ele)%value(y_limit$) = y_lim

  if (s_len /= ring%ele_(ix_ele)%value(l$)) then
    ring%ele_(ix_ele)%value(l$) = s_len
    call s_calc(ring)
  endif

end subroutine
