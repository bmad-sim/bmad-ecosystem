!+
! Subroutine control_bookkeeper (ring, ix_ele)
!
! Subroutine to transfer attibute information from lord to slave elements.
! It is assumend that the element with index ix_ele is the only element whose
! attribute values have been changed.
!
! Modules needed:
!   use bmad
!
! Input:
!   RING   -- Ring_struct: Ring to be used
!   IX_ELE -- Integer: Index of element whose attribute values have been
!               changed.
!-

!$Id$
!$Log$
!Revision 1.10  2003/03/18 20:36:44  dcs
!%num_steps for a slave now proportional to length
!
!Revision 1.9  2003/03/14 21:19:17  dcs
!Split bend bug fixed
!
!Revision 1.8  2003/01/27 14:40:32  dcs
!bmad_version = 56
!
!Revision 1.7  2002/06/13 14:54:24  dcs
!Interfaced with FPP/PTC
!
!Revision 1.6  2002/02/23 20:32:13  dcs
!Double/Single Real toggle added
!
!Revision 1.5  2002/01/08 21:44:38  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.4  2001/11/29 19:39:52  helms
!Updates from DCS including (*) -> (:)
!
!Revision 1.3  2001/10/02 18:49:11  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.2  2001/09/27 18:31:49  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine control_bookkeeper (ring, ix_ele)

  use bmad_struct
  use bmad_interface
  use multipole_mod

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  integer ix_ele, i, j, ix, ix1, ix2
  integer ix_eles(300)

! attribute bookkeeping

  call attribute_bookkeeper (ring%ele_(ix_ele), ring%param)

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

! First: Makup lords and group slaves.
! If a super_lord has lords in turn then these lords must be overlay_lords. 
! therefore treat the super_lord as an overlay_slave.
! the same is true if an overlay lord has lords.

  do j = 1, ix2

    ix = ix_eles(j)
    ele => ring%ele_(ix)

    if (ele%control_type == group_lord$) then
      call makeup_group_slaves (ring, ix)

    elseif (ele%control_type == super_lord$ .and. ele%n_lord > 0) then
      call adjust_super_lord_s_position (ring, ix)
      call makeup_overlay_slave (ring, ix)

    elseif (ele%control_type == overlay_lord$ .and. ele%n_lord > 0) then
      call makeup_overlay_slave (ring, ix)

    endif

  enddo

! Second: Makeup all slaves but group slaves.

  do j = 1, ix2

    ix = ix_eles(j)
    ele => ring%ele_(ix)       

    if (ele%control_type == super_slave$) then
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
! Subroutine adjust_super_lord_s_position (ring, ix_lord)
!
! Subroutine to adjust the positions of the slaves of a super_lord due
! to changes in the lord's s_offset.
!-

Subroutine adjust_super_lord_s_position (ring, ix_lord)

  use bmad

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave 

  integer ix_lord, ix

  real(rdef) s_start, s_start2, s_end

!

  lord => ring%ele_(ix_lord)

  if (lord%control_type /= super_lord$) then
    print *, 'ERROR IN ADJUST_SUPER_LORD_S_POSITION: ELEMENT IS NOT A LORD!'
    call err_exit
  endif

! If a super lord is moved then we just need to adjust the start and end edges.
! Adjust end position.

  s_end = lord%s + lord%value(s_offset$)
  ix = ring%control_(lord%ix2_slave)%ix_slave
  slave => ring%ele_(ix)
  s_start = slave%s - slave%value(l$)
  slave%value(l$) = s_end - s_start
  slave%s = s_end

! Adjust start position

  s_start = s_end - lord%value(l$)
  if (s_start < 0) s_start = s_start + ring%param%total_length
  ix = ring%control_(lord%ix1_slave)%ix_slave
  slave => ring%ele_(ix)
  s_start2 = slave%s - slave%value(l$)
  if (s_start < 0) s_start = s_start + ring%param%total_length 
  slave%value(l$) = slave%value(l$) + s_start2 - s_start

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

  use bmad
  use dcslib_interface

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave

  real(rdef) delta, coef

  integer ix_lord, ix, iv, ict, i

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

  use bmad
  use dcslib_interface

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave
  type (ele_struct) :: sol_quad

  integer i, ix_con, j, ix, ix_slave

  real(rdef) tilt, k1_x, k1_y, x_kick, y_kick, ks, k1, coef
  real(rdef) x_o, y_o, x_p, y_p, s_slave, s_del
  real(rdef) sin_2, cos_2, x_off, y_off, a(0:n_pole_maxx), b(0:n_pole_maxx)
  real(rdef) knl(0:n_pole_maxx), t(0:n_pole_maxx), value(n_attrib_maxx)
  real(rdef) a_tot(0:n_pole_maxx), b_tot(0:n_pole_maxx)
  real(rdef) sum_1, sum_2, sum_3, sum_4, ks_sum, ks_xp_sum, ks_xo_sum
  real(rdef) ks_yp_sum, ks_yo_sum, l_slave, r_off(4)
  real(rdef) t_1(4), t_2(4), T_end(4,4), mat4(4,4), mat4_inv(4,4), beta(4)
  real(rdef) T_tot(4,4), x_o_sol, x_p_sol, y_o_sol, y_p_sol

  logical, save :: init_needed = .true.

! init

  if (init_needed) then
    call init_ele (sol_quad)
    sol_quad%key = sol_quad$
    init_needed = .false.
  endif

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
! 1 super_lord for this super_slave: just transfer attributes except length

  if (slave%n_lord == 1) then

    ix_con = ring%ic_(slave%ic1_lord)  
    ix = ring%control_(ix_con)%ix_lord
    lord => ring%ele_(ix)
    coef = ring%control_(ix_con)%coef  ! = len_slave / len_lord

    value = lord%value
    value(l$) = slave%value(l$)                 ! do not change slave length
    value(hkick$) = lord%value(hkick$) * coef
    value(vkick$) = lord%value(vkick$) * coef
    slave%num_steps = lord%num_steps * coef + 1
    if (slave%key == rfcavity$) value(volt$) = lord%value(volt$) * coef
    if (ix_con /= lord%ix2_slave) then   ! if not at end of the lord domain
      value(x_limit$) = 0
      value(y_limit$) = 0
    endif

! s_del is the distance between lord and slave centers

    if (value(x_pitch$) /= 0 .or. value(y_pitch$) /= 0) then
      s_del = (slave%s - slave%value(l$)/2) - &
                  (lord%s + lord%value(s_offset$) - lord%value(l$)/2)
      s_del = modulo2 (s_del, ring%param%total_length/2)
      value(x_offset$) = value(x_offset$) + s_del * value(x_pitch$)
      value(y_offset$) = value(y_offset$) + s_del * value(y_pitch$)
    endif
      
    slave%value = value
    slave%is_on = lord%is_on

! if a wiggler: 
! must keep track of where we are in terms of the unsplit wiggler.
! This is for anything which does not try to make a homogeneous approximation.
! l_original is the length of the unsplit original wiggler.
! l_start is the starting point with respect to the original wiggler.
! l_end is the ending point with respect to the original wiggler.

    if (slave%key == wiggler$) then
      slave%value(n_pole$) = lord%value(n_pole$) * coef
      slave%value(l_original$) = lord%value(l$)
      slave%value(l_start$)    = (slave%s - slave%value(l$)) - &
                                                   (lord%s - lord%value(l$))
      slave%value(l_end$)      = slave%value(l_start$) + slave%value(l$)

      if (associated(lord%wig_term)) then
        if (.not. associated (slave%wig_term) .or. &
                size(slave%wig_term) /= size(lord%wig_term)) then
          deallocate (slave%wig_term)
          allocate (slave%wig_term(size(lord%wig_term)))
        endif
        do i = 1, size(lord%wig_term)
          slave%wig_term(i)%coef = lord%wig_term(i)%coef
          slave%wig_term(i)%kx = lord%wig_term(i)%kx
          slave%wig_term(i)%ky = lord%wig_term(i)%ky
          slave%wig_term(i)%kz = lord%wig_term(i)%kz
          slave%wig_term(i)%phi_z = lord%wig_term(i)%phi_z + &
                               lord%wig_term(i)%kz * slave%value(l_start$)
        enddo
      endif

    endif

! if an sbend:
!     1) renormalize the angles
!     2) zero the face angles next to the split

    if (slave%key == sbend$) then
      if (ix_con == lord%ix1_slave) then   ! first slave bend
        slave%value(e2$) = 0
      elseif (ix_con == lord%ix2_slave) then 
        slave%value(e1$) = 0
      else
        slave%value(e1$) = 0
        slave%value(e2$) = 0
      endif
    endif                       

    return

  endif

!-----------------------------------------------------------------------
! Multiple super_lords for this super_slave: 
! must be a solenoid/quadrupole combo.
! combine the lord elements.
                                           
  k1_x = 0
  k1_y = 0
  x_kick = 0
  y_kick = 0
  a_tot = 0
  b_tot = 0
  sum_1 = 0
  sum_2 = 0
  sum_3 = 0
  sum_4 = 0
  ks_sum = 0
  ks_xp_sum = 0
  ks_xo_sum = 0
  ks_yp_sum = 0
  ks_yo_sum = 0

  do j = 1, n_attrib_maxx
    if (j /= l$) slave%value(j) = 0
  enddo

  s_slave = slave%s - slave%value(l$)/2  ! center of slave

! sum over all lords...

  do j = slave%ic1_lord, slave%ic2_lord

    ix_con = ring%ic_(j)
    ix = ring%control_(ix_con)%ix_lord
    coef = ring%control_(ix_con)%coef
    lord => ring%ele_(ix)

    if (lord%control_type /= super_lord$) then
      type *, 'ERROR IN MAKEUP_SUPER_SLAVE: SUPER_SLAVE HAS A'
      type *, '      CONTROL ELEMENT THAT IS NOT SUPER_LORD'
      type *, '      SLAVE: ', slave%name, ix_slave
      type *, '      LORD:  ', lord%name, ix
      call err_exit
    endif

    if (ix_con == lord%ix2_slave) then      ! longitudinal ends match
      slave%value(x_limit$) = lord%value(x_limit$)
      slave%value(y_limit$) = lord%value(y_limit$)
    endif

    if (lord%is_on) then

      x_p = lord%value(x_pitch$);  x_o = lord%value(x_offset$)
      y_p = lord%value(y_pitch$);  y_o = lord%value(y_offset$)

      s_del = s_slave - (lord%s + lord%value(s_offset$) - lord%value(l$)/2)
      s_del = modulo2 (s_del, ring%param%total_length/2)

      ks = lord%value(ks$)

      ks_sum = ks_sum + ks

      ks_xp_sum = ks_xp_sum + ks * x_p
      ks_yp_sum = ks_yp_sum + ks * y_p

      ks_xo_sum = ks_xo_sum + ks * (x_o + x_p * s_del)
      ks_yo_sum = ks_yo_sum + ks * (y_o + y_p * s_del)

      x_kick = x_kick + lord%value(hkick$) * coef
      y_kick = y_kick + lord%value(vkick$) * coef

      tilt = lord%value(tilt$)
      cos_2 = lord%value(k1$) * cos(2 * tilt)
      sin_2 = lord%value(k1$) * sin(2 * tilt)

      k1_x = k1_x + cos_2
      k1_y = k1_y + sin_2

      sum_1 = sum_1 + cos_2 * x_p + sin_2 * y_p
      sum_2 = sum_2 + sin_2 * x_p - cos_2 * y_p

      sum_3 = sum_3 + cos_2 * (x_o + x_p * s_del) + sin_2 * (y_o + y_p * s_del)
      sum_4 = sum_4 + sin_2 * (x_o + x_p * s_del) - cos_2 * (y_o + y_p * s_del)

      if (associated(lord%a)) then
        call multipole_ele_to_kt (lord, +1, knl, t, .true.)
        call multipole_kt_to_ab (knl/lord%value(l$), t, a, b)
        a_tot = a_tot + a
        b_tot = b_tot + b
      endif

    endif

  enddo

! stuff sums into slave element

  ks = ks_sum
  slave%value(ks$) = ks
  slave%value(hkick$) = x_kick
  slave%value(vkick$) = y_kick

  if (k1_x == 0 .and. k1_y == 0 .and. ks == 0) return

  if (ks /= 0) then
    x_o_sol = ks_xo_sum / ks
    x_p_sol = ks_xp_sum / ks
    y_o_sol = ks_yo_sum / ks
    y_p_sol = ks_yp_sum / ks
  endif

  if (k1_x == 0 .and. k1_y == 0) then  ! pure solenoid
    slave%value(k1$) = 0
    slave%value(tilt$) = 0
    deallocate (slave%a, slave%b, stat = ix)
    slave%value(x_offset$) = x_o_sol
    slave%value(y_offset$) = y_o_sol
    slave%value(x_pitch$)  = x_p_sol
    slave%value(y_pitch$)  = y_p_sol
    return
  endif   

! here if have quadrupole component

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

  cos_2 = k1_x / (k1_x**2 + k1_y**2)
  sin_2 = k1_y / (k1_x**2 + k1_y**2)

  slave%value(x_pitch$)  = cos_2 * sum_1 + sin_2 * sum_2
  slave%value(y_pitch$)  = sin_2 * sum_1 - cos_2 * sum_2
  slave%value(x_offset$) = cos_2 * sum_3 + sin_2 * sum_4
  slave%value(y_offset$) = sin_2 * sum_3 - cos_2 * sum_4

  deallocate (slave%a, slave%b, stat = ix)
  if (any(a_tot /= 0) .or. any(b_tot /= 0)) then
    allocate (slave%a(0:n_pole_maxx), slave%b(0:n_pole_maxx))
    call multipole_ab_to_kt(a_tot, b_tot, knl, t)
    call multipole_kt_to_ab(knl/k1, t-tilt, a, b)
    slave%a = a
    slave%b = b
    slave%value(radius$) = 1
  endif

! if ks /= 0 then we have to recalculate the offsets and pitches.

  if (ks == 0) return

  x_p = slave%value(x_pitch$) - x_p_sol; x_o = slave%value(x_offset$) - x_o_sol
  y_p = slave%value(y_pitch$) - y_p_sol; y_o = slave%value(y_offset$) - y_o_sol

  if (x_p == 0 .and. x_o == 0 .and. y_p == 0 .and. y_o == 0) return

  t_2 = (/ x_o, x_p, y_o, y_p /)
  call tilt_coords (tilt, t_2, .true.)  ! set

  l_slave = slave%value(l$)

  t_1 = (/ t_2(2), 0.0_rdef, t_2(4), 0.0_rdef /)
  t_2(1) = t_2(1) + ks * t_2(4) / k1 
  t_2(3) = t_2(3) + ks * t_2(2) / k1
             
  call mat_unit (T_end, 4, 4)
  T_end(4,1) =  ks / 2
  T_end(2,3) = -ks / 2

  sol_quad%value(ks$) = ks
  sol_quad%value(k1$) = k1
  sol_quad%value(l$)  = l_slave
  call make_mat6 (sol_quad, ring%param)
  T_tot = sol_quad%mat6(1:4,1:4)

  r_off = matmul (T_end, l_slave * t_1 / 2 - t_2) 
  r_off = matmul (T_tot, r_off) + matmul (T_end, l_slave * t_1 / 2 + t_2)

  call mat_unit (mat4, 4, 4)
  mat4(:,2) = mat4(:,2) + l_slave * T_tot(:,1) / 2
  mat4(:,4) = mat4(:,4) + l_slave * T_tot(:,3) / 2
  mat4(1,2) = mat4(1,2) + l_slave / 2
  mat4(3,4) = mat4(3,4) + l_slave / 2
  mat4 = mat4 - T_tot

  call mat_inverse (mat4, mat4_inv)
  beta = matmul (mat4_inv, r_off)

  call tilt_coords (tilt, beta, .false.)  ! unset

  slave%value(x_offset$) = beta(1) + x_o_sol
  slave%value(x_pitch$)  = beta(2) + x_p_sol
  slave%value(y_offset$) = beta(3) + y_o_sol
  slave%value(y_pitch$)  = beta(4) + y_p_sol

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine makeup_overlay_slave (ring, ix_ele)

  use bmad                        

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  real(rdef) value(n_attrib_maxx), coef
  integer i, j, ix, iv, ix_ele, icom, ct
  logical used(n_attrib_maxx)

!
                               
  ele => ring%ele_(ix_ele)
  ct = ele%control_type

  if (ct /= super_lord$ .and. ct /= overlay_slave$ .and. &
                                               ct /= overlay_lord$) then
    type *, 'ERROR IN MAKEUP_OVERLAY_SLAVE: ELEMENT IS NOT OF PROPER TYPE.'
    type *, '      RING INDEX:', ix_ele
    call type_ele (ele, .true., 0, .false., 0, .true., ring)
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
      call type_ele (ele, .true., 0, .false., 0, .true., ring)
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
