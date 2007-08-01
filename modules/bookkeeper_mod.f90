#include "CESR_platform.inc"

module bookkeeper_mod

  use bmad_interface
  use bmad_utils_mod
  use multipole_mod

  integer, parameter :: off$ = 1, on$ = 2
  integer, parameter :: save_state$ = 3, restore_state$ = 4

  private this_bookkeeper
        
contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine control_bookkeeper (lattice, ix_ele)
!
! Subroutine to transfer attibute information from lord to slave elements.
! This subroutine will call attribute_bookkeeper.
! Note: To do a complete bookkeeping job on a lattice use:
!   lattice_bookkeeper
!
! Modules needed:
!   use bmad
!
! Input:
!   lattice   -- lat_struct: lattice to be used
!   ix_ele -- Integer, optional: Index of element whose attribute values 
!               have been changed. If not present bookkeeping will be done 
!               for all elements.
!-

subroutine control_bookkeeper (lattice, ix_ele)

  implicit none

  type (lat_struct), target :: lattice

  integer, optional :: ix_ele
  integer ie

! If ix_ele is present we only do bookkeeping for this one element

  if (present(ix_ele)) then
    call this_bookkeeper (lattice, ix_ele)

! Else we need to make up all the lords. 
! The group lords must be done last since their slaves don't know about them.

  else
    do ie = lattice%n_ele_track+1, lattice%n_ele_max
      if (lattice%ele(ie)%control_type /= group_lord$) call this_bookkeeper (lattice, ie)
    enddo

    do ie = lattice%n_ele_track+1, lattice%n_ele_max
      if (lattice%ele(ie)%control_type == group_lord$) call this_bookkeeper (lattice, ie)
    enddo

  endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine this_bookkeeper (lattice, ix_ele)
!
! This subroutine is only to be called from control_bookkeeper and is
! not meant for general use.
!-

subroutine this_bookkeeper (lattice, ix_ele)

  type (lat_struct), target :: lattice
  type (ele_struct), pointer :: lord, slave

  integer ix_ele, j, k, ix, ix1, ix2, ix_lord
  integer, allocatable, save :: ix_slaves(:), ix_super(:)

! Init

  call re_allocate (ix_slaves, lattice%n_ele_max)
  call re_allocate (ix_super, lattice%n_ele_max)
  ix_super = 0

! Attribute bookkeeping for this element

  call attribute_bookkeeper (lattice%ele(ix_ele), lattice%param)
  if (lattice%ele(ix_ele)%n_slave == 0 .and. &
            lattice%ele(ix_ele)%n_lord == 0) return  ! nothing more to do

! Make a list of slave elements to update.
! we do not need to update free elements of group lords.

  ix1 = 0   ! index for processed elements
  ix_slaves(1) = ix_ele
  ix2 = 1   ! index for last element in list

  do
    ix1 = ix1 + 1
    ix_lord = ix_slaves(ix1)
    lord => lattice%ele(ix_lord)
    do j = lord%ix1_slave, lord%ix2_slave
      ix = lattice%control(j)%ix_slave
      if (lord%control_type == group_lord$ .and. lattice%ele(ix)%control_type == free$) cycle
      if (ix == ix_slaves(ix2)) cycle   ! do not use duplicates
      ix2 = ix2 + 1
      ix_slaves(ix2) = ix
      if (lord%control_type == super_lord$) ix_super(ix2) = ix_lord
    enddo
    if (ix1 == ix2) exit
  enddo

! First: Makup lords and group slaves.
! If an overlay_lord has lords above it then these lords must be overlay_lords.
! Therefore treat the overlay_lord as an overlay_slave.
! The same is true if a super_lord has lords except in this case the lord
! may be a multipass_lord.

  do j = 1, ix2

    ix = ix_slaves(j)
    slave => lattice%ele(ix)

    if (slave%control_type == group_lord$) then
      call makeup_group_slaves (lattice, ix)
      call attribute_bookkeeper (slave, lattice%param)

    elseif (slave%control_type == super_lord$) then

      if (slave%n_lord > 0) then
        k =  lattice%ic(slave%ic1_lord)
        lord => lattice%ele(lattice%control(k)%ix_lord)
        if (lord%control_type == multipass_lord$) then
          call makeup_multipass_slave (lattice, ix)
        else
          call adjust_super_lord_s_position (lattice, ix)
          call makeup_overlay_and_girder_slave (lattice, ix)
        endif
        call attribute_bookkeeper (slave, lattice%param)
      endif

    elseif (slave%control_type == multipass_lord$ .and. slave%n_lord > 0) then
      call makeup_overlay_and_girder_slave (lattice, ix)
      call attribute_bookkeeper (slave, lattice%param)

    elseif (slave%control_type == overlay_lord$ .and. slave%n_lord > 0) then
      call makeup_overlay_and_girder_slave (lattice, ix)
      call attribute_bookkeeper (slave, lattice%param)

    endif

  enddo

! Second: Makeup all slaves but group slaves.

  do j = 1, ix2

    ix = ix_slaves(j)
    slave => lattice%ele(ix)       

    if (slave%control_type == super_slave$) then
      call makeup_super_slave (lattice, ix)
      call attribute_bookkeeper (slave, lattice%param)

    elseif (slave%control_type == overlay_slave$) then
      call makeup_overlay_and_girder_slave (lattice, ix)
      call attribute_bookkeeper (slave, lattice%param)

    elseif (slave%control_type == multipass_slave$) then
      call makeup_multipass_slave (lattice, ix)
      call attribute_bookkeeper (slave, lattice%param)

    endif

  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lattice_bookkeeper (lattice)
!
! Subroutine to do a complete bookkeeping job on a lattice.
!
! This this routine does a complete job of bookking and could be unacceptably
! slow if used, for example, in the inner loop of an optimizer. In this case
! consider using only control_bookkeeper instead.
!
! Modules needed:
!   use bmad
!
! Input:
!   lattice   -- lat_struct: Lattice needing bookkeeping.
!
! Output:
!   lattice   -- lat_struct: Lattice with bookkeeping done.
!-

subroutine lattice_bookkeeper (lattice)

  implicit none

  type (lat_struct) lattice
  integer i

! Control bookkeeper is called twice to make sure that the z_patch for a 
! wiggler super_lord is computed.

  call s_calc (lattice)
  call lat_geometry (lattice)
  call control_bookkeeper (lattice)
  call compute_reference_energy (lattice)
  call control_bookkeeper (lattice)

  do i = 1, lattice%n_ele_track
    call attribute_bookkeeper (lattice%ele(i), lattice%param)
  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine adjust_super_lord_s_position (lattice, ix_lord)
!
! Subroutine to adjust the positions of the slaves of a super_lord due
! to changes in the lord's s_offset.
!-

Subroutine adjust_super_lord_s_position (lattice, ix_lord)

  implicit none

  type (lat_struct), target :: lattice
  type (ele_struct), pointer :: lord, slave 

  integer ix_lord, ix

  real(rp) s_start, s_start2, s_end

  character(40) :: r_name = 'adjust_super_lord_s_position'

!

  lord => lattice%ele(ix_lord)

  if (lord%control_type /= super_lord$) then
     call out_io (s_abort$, r_name, 'ELEMENT IS NOT A LORD! ' // lord%name)
     call err_exit
  endif

! If a super lord is moved then we just need to adjust the start and end edges.
! Adjust end position.

  s_end = lord%s + lord%value(s_offset$)
  ix = lattice%control(lord%ix2_slave)%ix_slave
  slave => lattice%ele(ix)
  s_start = slave%s - slave%value(l$)
  slave%value(l$) = s_end - s_start
  slave%s = s_end

! Adjust start position

  s_start = s_end - lord%value(l$)
  if (s_start < 0) s_start = s_start + lattice%param%total_length
  ix = lattice%control(lord%ix1_slave)%ix_slave
  slave => lattice%ele(ix)
  s_start2 = slave%s - slave%value(l$)
  if (s_start < 0) s_start = s_start + lattice%param%total_length 
  slave%value(l$) = slave%value(l$) + s_start2 - s_start

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_group_slaves (lattice, ix_slave)
!
! Subroutine to calculate the attributes of group slave elements
!-

Subroutine makeup_group_slaves (lattice, ix_lord)   

  implicit none

  type (lat_struct), target :: lattice
  type (ele_struct), pointer :: lord, slave

  real(rp) delta, coef

  integer ix_lord, ix, iv, ict, i

  logical moved

  character(20) :: r_name = 'makeup_group_slaves'

!

  lord => lattice%ele(ix_lord)

  delta = lord%value(command$) - lord%value(old_command$)    ! change
  lord%value(old_command$) = lord%value(command$) ! save old

  moved = .false.   ! have we longitudinally moved an element?

  do i = lord%ix1_slave, lord%ix2_slave

    ix = lattice%control(i)%ix_slave
    iv = lattice%control(i)%ix_attrib
    ict = lattice%ele(ix)%control_type
    slave => lattice%ele(ix)

    if (iv == l$) then
      moved = .true.
      if (ict /= free$ .and. ict /= super_slave$) then
        call out_io (s_abort$, r_name, "A GROUP: " // lord%name, &
                    "CONTROLS THE LENGTH OF A LORD ELEMENT: " // slave%name)
        call err_exit
      endif
    endif
    coef = lattice%control(i)%coef
    slave%value(iv) = slave%value(iv) + delta * coef
  enddo

  if (moved) call s_calc (lattice)       ! recompute s distances


end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_multipass_slave (lattice, ix_slave)
!
! Subroutine to calcualte the attributes of multipass slave elements.
! This routine is not meant for general use.
!-

subroutine makeup_multipass_slave (lattice, ix_slave)

  implicit none

  type (lat_struct), target :: lattice
  type (ele_struct), pointer :: lord, slave

  real(rp) s, val(n_attrib_maxx)
  integer j, ix_slave

!

  slave => lattice%ele(ix_slave)
  j =  lattice%ic(slave%ic1_lord)
  lord => lattice%ele(lattice%control(j)%ix_lord)

  val = slave%value

  slave%value = lord%value
  if (lord%key == lcavity$ .or. lord%key == rfcavity$) then
    slave%value(dphi0$)        = val(dphi0$)
    slave%value(E_tot_start$) = val(E_tot_start$)
    slave%value(p0c_start$)    = val(p0c_start$)
  endif

  slave%value(E_tot$) = val(E_tot$)
  slave%value(p0c$)         = val(p0c$)

  if (associated (slave%a_pole)) then
    slave%a_pole = lord%a_pole
    slave%b_pole = lord%b_pole
  endif

  if (associated (slave%r)) slave%r = lord%r
  if (associated (slave%const)) slave%const = lord%const
  if (associated (slave%wake)) then
    slave%wake%sr_table       = lord%wake%sr_table
    slave%wake%sr_mode_long  = lord%wake%sr_mode_long
    slave%wake%sr_mode_trans = lord%wake%sr_mode_trans
    slave%wake%lr        = lord%wake%lr
  endif

  slave%mat6_calc_method = lord%mat6_calc_method
  slave%tracking_method  = lord%tracking_method
  slave%num_steps        = lord%num_steps      
  slave%is_on            = lord%is_on
  slave%aperture_at      = lord%aperture_at
  slave%coupler_at       = lord%coupler_at

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave (lattice, ix_slave)
!
! Subroutine to calcualte the attributes of superposition slave elements.
! This routine is not meant for general use.
!-
         
subroutine makeup_super_slave (lattice, ix_slave)

  implicit none

  type (lat_struct), target :: lattice
  type (ele_struct), pointer :: lord, slave
  type (ele_struct), save :: sol_quad

  integer i, j, ix_con, ix, ix_slave

  real(rp) tilt, k_x, k_y, x_kick, y_kick, ks, k1, coef
  real(rp) x_o, y_o, x_p, y_p, s_slave, s_del, k2, k3, c, s
  real(rp) sin_n, cos_n, a(0:n_pole_maxx), b(0:n_pole_maxx)
  real(rp) knl(0:n_pole_maxx), t(0:n_pole_maxx), value(n_attrib_maxx)
  real(rp) a_tot(0:n_pole_maxx), b_tot(0:n_pole_maxx)
  real(rp) sum_1, sum_2, sum_3, sum_4, ks_sum, ks_xp_sum, ks_xo_sum
  real(rp) ks_yp_sum, ks_yo_sum, l_slave, r_off(4), leng
  real(rp) t_1(4), t_2(4), T_end(4,4), mat4(4,4), mat4_inv(4,4), beta(4)
  real(rp) T_tot(4,4), x_o_sol, x_p_sol, y_o_sol, y_p_sol

  logical, save :: init_needed = .true.

  character(20) :: r_name = 'makeup_super_slave'

! init

  if (init_needed) then
    call init_ele (sol_quad)
    sol_quad%key = sol_quad$
    init_needed = .false.
  endif

! Super_slave:

  slave => lattice%ele(ix_slave)

  if (slave%control_type /= super_slave$) then
     call out_io(s_abort$, r_name, "ELEMENT IS NOT AN SUPER SLAVE: " // slave%name)
    call err_exit
  endif

! If this slave is the last slave for some lord (so that then longitudinal 
! end of the slave matches the end of the lord) then the limits of the lord
! are transfered to the slave.

!-----------------------------------------------------------------------
! 1 super_lord for this super_slave: just transfer attributes except length

  if (slave%n_lord == 1) then

    ix_con = lattice%ic(slave%ic1_lord)  
    ix = lattice%control(ix_con)%ix_lord
    lord => lattice%ele(ix)
    coef = lattice%control(ix_con)%coef  ! = len_slave / len_lord

    ! If this is not the first slave: Transfer reference orbit from previous slave

    if (ix_con /= lord%ix1_slave) then
      slave%ref_orb_in = lattice%ele(ix_slave-1)%ref_orb_out
    endif

    !

    value = lord%value
    value(check_sum$) = slave%value(check_sum$) ! do not change the check_sum
    value(l$) = slave%value(l$)                 ! do not change slave length
    if (lord%key == wiggler$) then
      value(z_patch$) = slave%value(z_patch$)
    endif
    if (lord%key == hkicker$ .or. lord%key == vkicker$) then
      value(kick$) = lord%value(kick$) * coef
    else
      value(hkick$) = lord%value(hkick$) * coef
      value(vkick$) = lord%value(vkick$) * coef
    endif
    slave%num_steps = lord%num_steps * coef + 1
    if (slave%key == rfcavity$) value(voltage$) = lord%value(voltage$) * coef

    slave%aperture_at = no_end$
    call compute_slave_aperture (value, slave, lord, ix_con)

    if (slave%key == lcavity$) then
      slave%coupler_at = no_end$
      call compute_slave_coupler (value, slave, lord, ix_con)
    endif

    ! s_del is the distance between lord and slave centers

    s_del = (slave%s - slave%value(l$)/2) - &
                  (lord%s + lord%value(s_offset$) - lord%value(l$)/2)
    s_del = modulo2 (s_del, lattice%param%total_length/2)
    value(x_pitch$)  = value(x_pitch_tot$)
    value(y_pitch$)  = value(y_pitch_tot$)
    value(x_offset$) = value(x_offset_tot$) + s_del * value(x_pitch_tot$)
    value(y_offset$) = value(y_offset_tot$) + s_del * value(y_pitch_tot$)
    value(tilt$)     = value(tilt_tot$)

    slave%value = value
    slave%is_on = lord%is_on
    slave%mat6_calc_method = lord%mat6_calc_method
    slave%tracking_method  = lord%tracking_method

! If a wiggler: 
! must keep track of where we are in terms of the unsplit wiggler.
! This is for anything which does not try to make a homogeneous approximation.
! l_original is the length of the unsplit original wiggler.
! l_start is the starting point with respect to the original wiggler.
! l_end is the ending point with respect to the original wiggler.

    if (slave%key == wiggler$) then
      slave%value(n_pole$) = lord%value(n_pole$) * coef
      slave%value(l_original$) = lord%value(l$)

      leng = 0 ! length of all slaves before this one
      do i = lord%ix1_slave, ix_con-1
        j = lattice%control(i)%ix_slave
        leng = leng + lattice%ele(j)%value(l$)
      enddo
      slave%value(l_start$)    = leng
      slave%value(l_end$)      = slave%value(l_start$) + slave%value(l$)

      if (associated(lord%wig_term)) then
        if (.not. associated (slave%wig_term) .or. &
                size(slave%wig_term) /= size(lord%wig_term)) then
          if (associated (slave%wig_term)) deallocate (slave%wig_term)
          allocate (slave%wig_term(size(lord%wig_term)))
        endif
        do i = 1, size(lord%wig_term)
          slave%wig_term(i) = lord%wig_term(i)
          slave%wig_term(i)%phi_z = lord%wig_term(i)%phi_z + &
                               lord%wig_term(i)%kz * slave%value(l_start$)
        enddo
      else
        if (associated (slave%wig_term)) deallocate (slave%wig_term)
      endif

    endif

! If a custom element: 
! Must keep track of where we are in terms of the unsplit element.
! See wiggler above for more details.

    if (slave%key == custom$) then
      slave%value(l_original$) = lord%value(l$)
      slave%value(l_start$)    = (slave%s - slave%value(l$)) - &
                                                   (lord%s - lord%value(l$))
      slave%value(l_end$)      = slave%value(l_start$) + slave%value(l$)
    endif

! If an sbend:
!     1) renormalize the angles
!     2) zero the face angles next to the split

    if (slave%key == sbend$) then
      if (ix_con == lord%ix1_slave) then   ! first slave bend
        slave%value(e2$)    = 0
        slave%value(h2$)    = 0
        slave%value(fintx$) = 0
        slave%value(hgapx$) = 0
      elseif (ix_con == lord%ix2_slave) then 
        slave%value(e1$)    = 0
        slave%value(h1$)    = 0
        slave%value(fint$)  = 0
        slave%value(hgap$)  = 0
      else
        slave%value(e1$)    = 0
        slave%value(h1$)    = 0
        slave%value(fint$)  = 0
        slave%value(hgap$)  = 0
        slave%value(e2$)    = 0
        slave%value(h2$)    = 0
        slave%value(fintx$) = 0
        slave%value(hgapx$) = 0
      endif
    endif                       

    return

  endif

!-----------------------------------------------------------------------
! Multiple super_lords for this super_slave: 
! combine the lord elements.
                                           
  k_x = 0
  k_y = 0
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

  value = 0
  value(l$) = slave%value(l$)
  value(E_tot$) = slave%value(E_tot$)
  value(p0c$) = slave%value(p0c$)
  value(check_sum$) = slave%value(check_sum$) ! do not change the check_sum

  s_slave = slave%s - value(l$)/2  ! center of slave
  slave%is_on = .false.

! sum over all lords...

  do j = slave%ic1_lord, slave%ic2_lord

    ix_con = lattice%ic(j)
    ix = lattice%control(ix_con)%ix_lord
    coef = lattice%control(ix_con)%coef
    lord => lattice%ele(ix)

    if (lord%control_type /= super_lord$) then
      call out_io (s_abort$, r_name, &
            "SUPER_SLAVE HAS A CONTROL ELEMENT THAT IS NOT A SUPER_LORD", &
            'SLAVE: ' //  slave%name // '  \i\ ', &
            'LORD:  ' //  lord%name  // '  \i\ ', i_array = (/ ix_slave, ix /) )
      call err_exit
    endif

    ! If this is not the first slave: Transfer reference orbit from previous slave

    if (ix_con /= lord%ix1_slave) then
      slave%ref_orb_in = lattice%ele(ix_slave-1)%ref_orb_out
    endif

    !

    call compute_slave_aperture (value, slave, lord, ix_con)
    if (slave%key == lcavity$) call compute_slave_coupler (value, slave, lord, ix_con)

    if (j == slave%ic1_lord) then
      slave%mat6_calc_method = lord%mat6_calc_method
      slave%tracking_method  = lord%tracking_method
    else
      if (slave%mat6_calc_method /= lord%mat6_calc_method) then
        ix = lattice%control(lattice%ic(slave%ic1_lord))%ix_lord
        call out_io(s_abort$, r_name, 'MAT6_CALC_METHOD DOES NOT AGREE FOR DIFFERENT', &
             'SUPERPOSITION LORDS: ' // trim(lord%name) // ', ' // trim(lattice%ele(ix)%name))
        call err_exit
      endif
      if (slave%tracking_method /= lord%tracking_method) then
        ix = lattice%control(lattice%ic(slave%ic1_lord))%ix_lord
        call out_io(s_abort$, r_name, ' TRACKING_METHOD DOES NOT AGREE FOR DIFFERENT', &
             'SUPERPOSITION LORDS: ' // trim(lord%name) // ', ' // trim(lattice%ele(ix)%name))
        call err_exit
      endif
    endif

    if (.not. lord%is_on) cycle
    slave%is_on = .true.  ! on if at least one lord is on

    tilt = lord%value(tilt_tot$)

    if (lord%key == hkicker$) then
      x_kick = x_kick + lord%value(kick$) * cos(tilt) * coef
      y_kick = y_kick + lord%value(kick$) * sin(tilt) * coef
    elseif (lord%key == vkicker$) then
      x_kick = x_kick - lord%value(kick$) * sin(tilt) * coef
      y_kick = y_kick + lord%value(kick$) * cos(tilt) * coef
    elseif (lord%key == kicker$) then
      c = cos(tilt) * coef
      s = sin(tilt) * coef
      x_kick = x_kick + c * lord%value(hkick$) - s * lord%value(vkick$)
      y_kick = y_kick + s * lord%value(hkick$) + c * lord%value(vkick$)
    else
      x_kick = x_kick + lord%value(hkick$) * coef
      y_kick = y_kick + lord%value(vkick$) * coef
    endif

    if (associated(lord%a_pole)) then
      call multipole_ele_to_kt (lord, +1, knl, t, .true.)
      call multipole_kt_to_ab (knl/lord%value(l$), t, a, b)
      a_tot = a_tot + a
      b_tot = b_tot + b
    endif

!------

    select case (slave%key)


! sextupole

    case (sextupole$) 

      cos_n = lord%value(k2$) * cos(3 * tilt)
      sin_n = lord%value(k2$) * sin(3 * tilt)            
      
      k_x = k_x + cos_n
      k_y = k_y + sin_n
  
! octupole

    case (octupole$)

      cos_n = lord%value(k3$) * cos(4 * tilt)
      sin_n = lord%value(k3$) * sin(4 * tilt)        
      
      k_x = k_x + cos_n
      k_y = k_y + sin_n

! solenoid/quadrupole combo.

    case (solenoid$, sol_quad$, quadrupole$)

      x_p = lord%value(x_pitch_tot$);  x_o = lord%value(x_offset_tot$)
      y_p = lord%value(y_pitch_tot$);  y_o = lord%value(y_offset_tot$)

      s_del = s_slave - (lord%s + lord%value(s_offset_tot$) - lord%value(l$)/2)
      s_del = modulo2 (s_del, lattice%param%total_length/2)

      ks = lord%value(ks$)

      ks_sum = ks_sum + ks

      ks_xp_sum = ks_xp_sum + ks * x_p
      ks_yp_sum = ks_yp_sum + ks * y_p

      ks_xo_sum = ks_xo_sum + ks * (x_o + x_p * s_del)
      ks_yo_sum = ks_yo_sum + ks * (y_o + y_p * s_del)

      cos_n = lord%value(k1$) * cos(2 * tilt)
      sin_n = lord%value(k1$) * sin(2 * tilt)

      k_x = k_x + cos_n
      k_y = k_y + sin_n

      sum_1 = sum_1 + cos_n * x_p + sin_n * y_p
      sum_2 = sum_2 + sin_n * x_p - cos_n * y_p

      sum_3 = sum_3 + cos_n * (x_o + x_p * s_del) + sin_n * (y_o + y_p * s_del)
      sum_4 = sum_4 + sin_n * (x_o + x_p * s_del) - cos_n * (y_o + y_p * s_del)

! bend_sol_quad

    case (bend_sol_quad$)
      call out_io (s_abort$, r_name, &
                 'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
      call err_exit

! hkicker, vkicker, kicker

    case (hkicker$, vkicker$, kicker$)

! default

    case default
      call out_io (s_abort$, r_name, &
                 'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
      call err_exit

    end select

  enddo

!------------------------------
! stuff sums into slave element

  if (x_kick == 0 .and. y_kick == 0) then
    if (slave%key == hkicker$ .or. slave%key == vkicker$) then
      value(kick$) = 0
    else
      value(hkick$) = 0
      value(vkick$) = 0
    endif
  elseif (slave%key == hkicker$) then
    value(kick$) = sqrt(x_kick**2 + y_kick**2)
    value(tilt$) = atan2(y_kick, x_kick)
  elseif (slave%key == vkicker$) then
    value(kick$) = sqrt(x_kick**2 + y_kick**2)
    value(tilt$) = atan2(-x_kick, y_kick)
  elseif (slave%key == kicker$) then
    value(tilt$) = 0
    value(hkick$) = x_kick
    value(vkick$) = y_kick
  else
    value(hkick$) = x_kick
    value(vkick$) = y_kick
  endif

  slave%value = value

  if (any(a_tot /= 0) .or. any(b_tot /= 0)) then
    call multipole_init(slave)
    call multipole_ab_to_kt(a_tot, b_tot, knl, t)
    call multipole_kt_to_ab(knl*slave%value(l$), t-tilt, a, b)
    slave%a_pole = a
    slave%b_pole = b
    slave%value(radius$) = 1
  elseif (associated(slave%a_pole)) then
    deallocate (slave%a_pole, slave%b_pole)
  endif

!-----------------------------

  select case (slave%key)

  case (sextupole$) 

    if (k_x == 0 .and. k_y == 0) return

    k2 = sqrt(k_x**2 + k_y**2)
    tilt = atan2(k_y, k_x) / 3

    if (tilt > pi/6) then
      k2 = -k2
      tilt = tilt - pi/3
    elseif (tilt < -pi/6) then
      k2 = -k2
      tilt = tilt + pi/3
    endif

    slave%value(k2$) = k2
    slave%value(tilt$) = tilt

! octupole

  case (octupole$)

    if (k_x == 0 .and. k_y == 0 .and. ks == 0) return

    k3 = sqrt(k_x**2 + k_y**2)
    tilt = atan2(k_y, k_x) / 4

    if (tilt > pi/8) then
      k3 = -k3
      tilt = tilt - pi/4
    elseif (tilt < -pi/8) then
      k3 = -k3
      tilt = tilt + pi/4
    endif

    slave%value(k3$) = k3
    slave%value(tilt$) = tilt

! sol_quad, etc.

  case (solenoid$, sol_quad$, quadrupole$)

    ks = ks_sum
    slave%value(ks$) = ks

    if (k_x == 0 .and. k_y == 0 .and. ks == 0) return

    if (ks /= 0) then
      x_o_sol = ks_xo_sum / ks
      x_p_sol = ks_xp_sum / ks
      y_o_sol = ks_yo_sum / ks
      y_p_sol = ks_yp_sum / ks
    endif

    if (k_x == 0 .and. k_y == 0) then  ! pure solenoid
      slave%value(k1$) = 0
      slave%value(tilt$) = 0
     deallocate (slave%a_pole, slave%b_pole, stat = ix)
      slave%value(x_offset$) = x_o_sol
      slave%value(y_offset$) = y_o_sol
      slave%value(x_pitch$)  = x_p_sol
      slave%value(y_pitch$)  = y_p_sol
      return
    endif   

! here if have quadrupole component

    k1 = sqrt(k_x**2 + k_y**2)
    tilt = atan2(k_y, k_x) / 2

    if (tilt > pi/4) then
      k1 = -k1
      tilt = tilt - pi/2
    elseif (tilt < -pi/4) then
      k1 = -k1
      tilt = tilt + pi/2
    endif

    slave%value(k1$) = k1
    slave%value(tilt$) = tilt

    cos_n = k_x / (k_x**2 + k_y**2)
    sin_n = k_y / (k_x**2 + k_y**2)

    slave%value(x_pitch$)  = cos_n * sum_1 + sin_n * sum_2
    slave%value(y_pitch$)  = sin_n * sum_1 - cos_n * sum_2
    slave%value(x_offset$) = cos_n * sum_3 + sin_n * sum_4
    slave%value(y_offset$) = sin_n * sum_3 - cos_n * sum_4

! if ks /= 0 then we have to recalculate the offsets and pitches.

    if (ks == 0) return

    x_p = slave%value(x_pitch$) - x_p_sol; x_o = slave%value(x_offset$) - x_o_sol
    y_p = slave%value(y_pitch$) - y_p_sol; y_o = slave%value(y_offset$) - y_o_sol

    if (x_p == 0 .and. x_o == 0 .and. y_p == 0 .and. y_o == 0) return

    t_2 = (/ x_o, x_p, y_o, y_p /)
    call tilt_coords (tilt, t_2, .true.)  ! set

    l_slave = slave%value(l$)

    t_1 = (/ t_2(2), 0.0_rp, t_2(4), 0.0_rp /)
    t_2(1) = t_2(1) + ks * t_2(4) / k1 
    t_2(3) = t_2(3) + ks * t_2(2) / k1
             
    call mat_make_unit (T_end)
    T_end(4,1) =  ks / 2
    T_end(2,3) = -ks / 2

    sol_quad%value(ks$) = ks
    sol_quad%value(k1$) = k1
    sol_quad%value(l$)  = l_slave
    call make_mat6 (sol_quad, lattice%param)
    T_tot = sol_quad%mat6(1:4,1:4)

    r_off = matmul (T_end, l_slave * t_1 / 2 - t_2) 
    r_off = matmul (T_tot, r_off) + matmul (T_end, l_slave * t_1 / 2 + t_2)

    call mat_make_unit (mat4)
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

! bend_sol_quad

  case (bend_sol_quad$)
    call out_io (s_abort$, r_name, &
                 'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
    call err_exit

  end select

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine compute_slave_aperture (value, slave, lord, ix_con)
!
! This routine is not meant for general use.
!-

subroutine compute_slave_aperture (value, slave, lord, ix_con)

  use bmad_struct

  implicit none

  type (ele_struct) slave, lord
  real(rp) value(n_attrib_maxx)
  integer ix_con

!

  select case (lord%aperture_at)
  case (exit_end$) 
    if (ix_con == lord%ix2_slave) slave%aperture_at = exit_end$
  case (entrance_end$)
    if (ix_con == lord%ix1_slave) slave%aperture_at = entrance_end$
  case (both_ends$)
    if (ix_con == lord%ix1_slave .and. ix_con == lord%ix2_slave) then
      slave%aperture_at = both_ends$
    elseif (ix_con == lord%ix1_slave) then
      slave%aperture_at = entrance_end$
    elseif (ix_con == lord%ix2_slave) then 
      slave%aperture_at = exit_end$
    endif
  end select

  if (slave%aperture_at == no_end$) then
    value(x_limit$) = 0
    value(y_limit$) = 0
    value(aperture$) = 0
  endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine compute_slave_coupler (value, slave, lord, ix_con)
!
! This routine is not meant for general use.
!-

subroutine compute_slave_coupler (value, slave, lord, ix_con)

  use bmad_struct

  implicit none

  type (ele_struct) slave, lord
  real(rp) value(n_attrib_maxx)
  integer ix_con

!

  select case (lord%coupler_at)
  case (exit_end$) 
    if (ix_con == lord%ix2_slave) slave%coupler_at = exit_end$
  case (entrance_end$)
    if (ix_con == lord%ix1_slave) slave%coupler_at = entrance_end$
  case (both_ends$)
    if (ix_con == lord%ix1_slave .and. ix_con == lord%ix2_slave) then
      slave%coupler_at = both_ends$
    elseif (ix_con == lord%ix1_slave) then
      slave%coupler_at = entrance_end$
    elseif (ix_con == lord%ix2_slave) then 
      slave%coupler_at = exit_end$
    endif
  end select

  if (slave%coupler_at == no_end$) then
    value(coupler_strength$) = 0
  endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_overlay_and_girder_slave (lattice, ix_ele)
!
! This routine is not meant for general use.
!-

subroutine makeup_overlay_and_girder_slave (lattice, ix_ele)

  implicit none

  type (lat_struct), target :: lattice
  type (ele_struct), pointer :: slave, lord

  real(rp) value(n_attrib_maxx), coef, ds
  integer i, j, ix, iv, ix_ele, icom, ct
  logical used(n_attrib_maxx)

  character(40) :: r_name = 'makeup_overlay_and_girder_slave'

!
                               
  slave => lattice%ele(ix_ele)
  ct = slave%control_type

  if (ct /= super_lord$ .and. ct /= overlay_slave$ .and. &
           ct /= overlay_lord$ .and. ct /= multipass_lord$) then
    call out_io(s_abort$, r_name, 'ELEMENT IS NOT OF PROPER TYPE. lattice INDEX: \i\ ', ix_ele)
    call type_ele (slave, .true., 0, .false., 0, .true., lattice)
    call err_exit
  endif

  value = 0
  used = .false.
  slave%on_an_girder = .false.

  do i = slave%ic1_lord, slave%ic2_lord
    j = lattice%ic(i)
    ix = lattice%control(j)%ix_lord
    lord => lattice%ele(ix)

    if (lord%control_type == multipass_lord$) cycle

    if (lord%control_type == girder_lord$) then
      ds = (slave%s - slave%value(l$)/2) - lord%value(s_center$) 
      slave%value(x_offset_tot$) = slave%value(x_offset$) + &
                     ds * lord%value(x_pitch$) + lord%value(x_offset$)
      slave%value(y_offset_tot$) = slave%value(y_offset$) + &
                     ds * lord%value(y_pitch$) + lord%value(y_offset$)
      slave%value(s_offset_tot$) = slave%value(s_offset$) + lord%value(s_offset$)
      slave%value(x_pitch_tot$)  = slave%value(x_pitch$)  + lord%value(x_pitch$)
      slave%value(y_pitch_tot$)  = slave%value(y_pitch$)  + lord%value(y_pitch$)
      slave%value(tilt_tot$)     = slave%value(tilt$)     + lord%value(tilt$)
      slave%on_an_girder = .true.
      cycle
    endif

    if (lord%control_type /= overlay_lord$) then
      call out_io (s_abort$, r_name, 'THE LORD IS NOT AN OVERLAY_LORD \i\ ', ix_ele)
      call type_ele (slave, .true., 0, .false., 0, .true., lattice)
      call err_exit
    endif     

    coef = lattice%control(j)%coef
    iv = lattice%control(j)%ix_attrib
    icom = lord%ix_value
    value(iv) = value(iv) + lord%value(icom)*coef
    used(iv) = .true.
  enddo

  where (used) slave%value = value

! If no girder then simply transfer tilt to tilt_tot, etc.

  if (.not. slave%on_an_girder) then
    slave%value(tilt_tot$)     = slave%value(tilt$)
    slave%value(x_offset_tot$) = slave%value(x_offset$)
    slave%value(y_offset_tot$) = slave%value(y_offset$)
    slave%value(s_offset_tot$) = slave%value(s_offset$)
    slave%value(x_pitch_tot$)  = slave%value(x_pitch$)
    slave%value(y_pitch_tot$)  = slave%value(y_pitch$)
  endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine attribute_bookkeeper (ele, param)
!
! Subroutine to recalculate the dependent attributes of an element.
! If the attributes have changed then any Taylor Maps will be killed.
! Note: This routine does not do any other bookkeeping. Consider using
! control_bookkeeper or lattice_bookkeeper instead.
! 
! BEAMBEAM:   
!     bbi_const$ = param%n_part * m_electron * charge$ * r_e /
!                           (2 * pi * p0c$ * (sig_x$ + sig_y$)
!
! ELSEPARATOR:
!     e_field$ = sqrt(hkick$**2 + vkick$**2) * p0c$ / l$
!     voltage$ = e_field$ * gap$ 
!
! LCAVITY:    
!     delta_e$ = gradient$ * L$ 
!     E_tot$   = E_tot$ + gradient$ * l$ * cos(phase)
!     p0c$     = sqrt(E_tot$**2 - mc2^2)
! 
! RFCAVITY:   
!     rf_frequency$ = harmon$ * c_light / param%total_length (only if harmon$ /= 0)
!
! SBEND:      
!     angle$   = L$ * G$
!     l_chord$ = 2 * sin(Angle$/2) / G$
!     rho$     = 1 / G$
!     k2$      = 2 * B(2) / L$
!
! WIGGLER:    
!     k1$  = -0.5 * (c_light * b_max$ / p0c$)**2
!     rho$ = p0c$ / (c_light * b_max$)
!     n_pole$ = L$ / l_pole$
!     z_patch$
!     x_patch$ (for a periodic wiggler)
!
! Modules needed:
!   use bmad
!
! Input:
!   ele        -- Ele_struct: Element with attributes 
!   param      -- lat_param_struct: 
!
! Output:
!   ele            -- Ele_struct: Element with self-consistant attributes.
!     %ref_orb_out -- Reference orbit to be used for the next
!                         super_slave wiggler z_patch calculation. 
!                         This is to be only used by the control_bookkeeper routine.
!
! Programming Note: If the dependent attributes are changed then 
!       the attribute_free routine must be modified.
!-

subroutine attribute_bookkeeper (ele, param)

  use symp_lie_mod, only: symp_lie_bmad

  implicit none

  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) start, end

  real(rp) factor, check_sum, phase, E_tot
  
  character(20) ::  r_name = 'attribute_bookkeeper'

! Transfer tilt to tilt_tot, etc.

  if (.not. ele%on_an_girder .and. ele%key /= match$) then
    ele%value(tilt_tot$)     = ele%value(tilt$)
    ele%value(x_offset_tot$) = ele%value(x_offset$)
    ele%value(y_offset_tot$) = ele%value(y_offset$)
    ele%value(s_offset_tot$) = ele%value(s_offset$)
    ele%value(x_pitch_tot$)  = ele%value(x_pitch$)
    ele%value(y_pitch_tot$)  = ele%value(y_pitch$)
  endif

! field_master

  if (ele%field_master) then

    if (ele%value(p0c$) == 0) then
      factor = 0
    else
      factor = c_light / ele%value(p0c$)
    endif

    select case (ele%key)
    case (quadrupole$)
      ele%value(k1$) = factor * ele%value(B1_gradient$)
    case (sextupole$)
      ele%value(k2$) = factor * ele%value(B2_gradient$)
    case (octupole$)
      ele%value(k3$) = factor * ele%value(B3_gradient$)
    case (solenoid$)
      ele%value(ks$) = factor * ele%value(Bs_field$)
    case (sol_quad$)
      ele%value(ks$) = factor * ele%value(Bs_field$)
      ele%value(k1$) = factor * ele%value(B1_gradient$)
    case (sbend$)
      ele%value(g$)     = factor * ele%value(B_field$)
      ele%value(g_err$) = factor * ele%value(B_field_err$)
    case (hkicker$)
      ele%value(kick$) = factor * ele%value(BL_kick$)
    case (vkicker$)
      ele%value(kick$) = factor * ele%value(BL_kick$)
    case (lcavity$, drift$, monitor$, instrument$, &
          ecollimator$, rcollimator$, elseparator$)
    case default
       call out_io(s_abort$,r_name,' "FIELD_MASTER" NOT IMPLEMENTED FOR: ' // trim(ele%name))
      call err_exit
    end select

    ele%value(hkick$) = factor * ele%value(BL_hkick$)
    ele%value(vkick$) = factor * ele%value(BL_vkick$)

  else

    factor = ele%value(p0c$) / c_light

    select case (ele%key)
    case (quadrupole$)
      ele%value(B1_gradient$) = factor * ele%value(k1$)
    case (sextupole$)
      ele%value(B2_gradient$) = factor * ele%value(k2$)
    case (octupole$)
      ele%value(B3_gradient$) = factor * ele%value(k3$)
    case (solenoid$)
      ele%value(Bs_field$)    = factor * ele%value(ks$)
    case (sol_quad$)
      ele%value(Bs_field$ )   = factor * ele%value(ks$)
      ele%value(B1_gradient$) = factor * ele%value(k1$)
    case (sbend$)
      ele%value(B_field$)     = factor * ele%value(g$)
      ele%value(B_field_err$) = factor * ele%value(g_err$)
    case (hkicker$)
      ele%value(BL_kick$) = factor * ele%value(kick$)
    case (vkicker$) 
      ele%value(BL_kick$) = factor * ele%value(kick$)
    end select

    ele%value(BL_hkick$) = factor * ele%value(hkick$)
    ele%value(BL_vkick$) = factor * ele%value(vkick$)

  endif

! Dependent attribute bookkeeping.
! Note: If the dependent attributes are changed then attribute_free 
!       must be modified.

  select case (ele%key)

! Bends

  case (sbend$)

    ele%value(angle$) = ele%value(l$) * ele%value(g$)

    if (ele%value(l$) == 0 .or. ele%value(g$) == 0) then
      ele%value(l_chord$) = 0
    else
      ele%value(l_chord$) = 2 * sin(ele%value(angle$)/2) / ele%value(g$)
    endif

    if (ele%value(g$) == 0) then
      ele%value(rho$) = 0
    else
      ele%value(rho$) = 1 / ele%value(g$)
    endif

    if (ele%value(l$) /= 0 .and. associated(ele%a_pole)) then
      ele%value(k2$) = 2 * ele%b_pole(2) / ele%value(l$)
    else
      ele%value(k2$) = 0
    endif

! Lcavity
! Only do the calculation if the starting energy is not zero since 
! attribute_bookkeeper can be called before the attributes are set.

  case (lcavity$)
    if (ele%value(E_tot_start$) /= 0) then
      ele%value(delta_e$) = ele%value(gradient$) * ele%value(L$) 
      phase = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
      E_tot = ele%value(E_tot_start$) + ele%value(gradient$) * &
                                                ele%value(l$) * cos(phase)
      E_tot = E_tot - ele%value(e_loss$) * param%n_part * e_charge
      ele%value(E_tot$) = E_tot
      call convert_total_energy_to (E_tot, param%particle, pc = ele%value(p0c$))
    endif

! RFcavity

  case (rfcavity$)
    if (ele%value(harmon$) /= 0) ele%value(rf_frequency$) =  &
                              ele%value(harmon$) * c_light / param%total_length 

! BeamBeam

  case (beambeam$)

    if (ele%value(n_slice$) == 0) ele%value(n_slice$) = 1.0 ! revert to default

    if (ele%value(charge$) == 0 .or. param%n_part == 0) then
      ele%value(bbi_const$) = 0

    else

      if (ele%value(sig_x$) == 0 .or. ele%value(sig_y$) == 0) then
        call out_io(s_abort$, r_name, 'ZERO SIGMA IN BEAMBEAM ELEMENT!')
        call type_ele(ele, .true., 0, .false., 0, .false.)
        call err_exit
      endif

      ele%value(bbi_const$) = &
        -param%n_part * m_electron * ele%value(charge$) * r_e /  &
        (2 * pi * ele%value(p0c$) * (ele%value(sig_x$) + ele%value(sig_y$)))

    endif

! Elseparator

  case (elseparator$)

    if (ele%value(l$) == 0 .or. ele%value(gap$) == 0) then
      ele%value(e_field$) = 0
      ele%value(voltage$) = 0
    else
      ele%value(e_field$) = sqrt(ele%value(hkick$)**2 + ele%value(vkick$)**2) * &
                                               ele%value(p0c$) / ele%value(l$)
      ele%value(voltage$) = ele%value(e_field$) * ele%value(gap$) 
    endif


! Wiggler
! Periodic_type wigglers have a single %wig_term for use with tracking, etc.

  case (wiggler$) 

    if (ele%value(p0c$) == 0) then
      ele%value(k1$) = 0
    else
      ele%value(k1$) = -0.5 * &
                    (c_light * ele%value(b_max$) / ele%value(p0c$))**2
    endif

    if (ele%value(b_max$) == 0) then
      ele%value(rho$) = 0
    else
      ele%value(rho$) = ele%value(p0c$) / (c_light * ele%value(b_max$))
    endif

    if (ele%value(l_pole$) == 0) then
      ele%value(n_pole$) = 0
    else
      ele%value(n_pole$) = ele%value(l$) / ele%value(l_pole$)
    endif

    if (ele%sub_key == periodic_type$) then
      if (.not. associated(ele%wig_term)) allocate (ele%wig_term(1))

      if (ele%value(l_pole$) == 0) then
        ele%wig_term(1)%ky = 0
      else
        ele%wig_term(1)%ky = pi / ele%value(l_pole$)
      endif
      ele%wig_term(1)%coef   = ele%value(b_max$)
      ele%wig_term(1)%kx     = 0
      ele%wig_term(1)%kz     = ele%wig_term(1)%ky
      ele%wig_term(1)%phi_z  = (ele%value(l_pole$) - ele%value(l$)) / 2
      ele%wig_term(1)%type   = hyper_y$
    endif

  end select

! num_steps

  if (ele%value(ds_step$) /= 0) ele%num_steps = abs(nint(ele%value(l$) / ele%value(ds_step$)))
  if (ele%num_steps == 0) ele%num_steps = 1

! We need to kill the Taylor Map, etc. if things have changed.
! calculate a check sum to see if things have changed.
! ele%value(check_sum$) == 0 means that the check_sum has never been 
! computed so in this case do not kill the Taylor Map

  select case (ele%key)

  case (wiggler$)
    check_sum = ele%value(polarity$)
    check_sum = check_sum

  case (quadrupole$)
    check_sum = ele%value(k1$) 

  case (sol_quad$)
    check_sum = ele%value(ks$) + ele%value(k1$)

  case (solenoid$)
    check_sum = ele%value(ks$)

  case (sbend$)
    check_sum = ele%value(g$) + ele%value(g_err$) + ele%value(e1$) + ele%value(e2$)

  case (sextupole$)
    check_sum = ele%value(k2$)

  case (octupole$)
    check_sum = ele%value(k3$)

  case (rfcavity$)
    check_sum = ele%value(voltage$) + ele%value(phi0$) 

  case (elseparator$)
    check_sum = check_sum + ele%value(hkick$) + ele%value(vkick$)

  case (lcavity$)
    check_sum = ele%value(gradient$) + ele%value(phi0$) + ele%value(gradient_err$) + &
                ele%value(phi0_err$)
  case default
    return

  end select

  check_sum = check_sum + ele%value(l$) + ele%value(ds_step$)

  if (ele%map_with_offsets) check_sum = check_sum + ele%value(x_offset$) + &
        ele%value(y_offset$) + ele%value(x_pitch$) + ele%value(y_pitch$) + &
        ele%value(s_offset$) + ele%value(tilt$)

  ! For some very strange reason there can be round off error in the check sum.
  ! Hence we use a non-exact test.

  if (abs(ele%value(check_sum$) - check_sum) > &
                      1d-14 * (abs(check_sum) + abs(ele%value(check_sum$)))) then

    if (ele%value(check_sum$) /= 0) then
      if (associated(ele%taylor(1)%term)) call kill_taylor(ele%taylor)
      if (associated(ele%gen_field)) call kill_gen_field(ele%gen_field)
      if (ele%key == wiggler$) then
        ele%value(z_patch$) = 0
        ele%value(x_patch$) = 0
      endif
    endif

    ele%value(check_sum$) = check_sum

  endif

! compute the z_patch for a wiggler if needed.
! The starting reference orbit is in ele%value(ref_orb$:ref_orb$+3).
! This is normally zero except for split wiggler sections.

  if (ele%key == wiggler$ .and. &
      ele%value(z_patch$) == 0 .and. ele%value(p0c$) /= 0) then
    start%vec = 0
    start%vec = ele%ref_orb_in
    call symp_lie_bmad (ele, param, start, end, .false., offset_ele = .false.)
    ele%value(z_patch$) = end%vec(5)
    if (ele%value(z_patch$) == 0) ele%value(z_patch$) = 1e-30 ! something non-zero.
    if (ele%sub_key == periodic_type$) ele%value(x_patch$) = end%vec(1)
    end%vec(5) = end%vec(5) - ele%value(z_patch$)
    ele%ref_orb_out = end%vec             ! save for next super_slave
  endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine changed_attribute_bookkeeper (lat, ix_ele, a_ptr)
!
! Subroutine to do bookkeeping when a particular attribute has been altered.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat    -- lat_struct: Lattice with the changed attribute.
!   ix_ele -- Integer: Index of element if an element attribute is being modified.
!               Otherwise should be set to -1.
!   a_ptr  -- Real(rp), pointer: Pointer to the changed attribute.
!
! Output:
!   lat  -- lat_struct: Lattice with appropriate changes.
!-

subroutine changed_attribute_bookkeeper (lat, ix_ele, a_ptr)

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

real(rp), pointer :: a_ptr
real(rp) v_mat(4,4), v_inv_mat(4,4), eta_vec(4), eta_xy_vec(4)

integer ix_ele

logical coupling_change

!

if (ix_ele > -1) then

  ele => lat%ele(ix_ele)

  if (associated(ele%taylor(1)%term)) call kill_taylor(ele%taylor)

  if (ele%ix_ele == 0) then
    coupling_change = .false.

    if (associated(a_ptr, ele%a%beta) .or. associated(a_ptr, ele%a%alpha)) then
      if (ele%a%beta /= 0) ele%a%gamma = (1 + ele%a%alpha**2) / ele%a%beta
      return
    endif

    if (associated(a_ptr, ele%b%beta) .or. associated(a_ptr, ele%b%alpha)) then
      if (ele%b%beta /= 0) ele%b%gamma = (1 + ele%b%alpha**2) / ele%b%beta
      return
    endif

    if (associated(a_ptr, ele%c_mat(1,1)) .or. associated(a_ptr, ele%c_mat(1,2)) .or. & 
            associated(a_ptr, ele%c_mat(2,1)) .or. associated(a_ptr, ele%c_mat(2,2))) then
      ele%gamma_c = sqrt(1 - ele%c_mat(1,1)*ele%c_mat(2,2) + &
                                                  ele%c_mat(1,2)*ele%c_mat(2,1))
      coupling_change = .true.
    endif

    if (associated(a_ptr, ele%x%eta) .or. associated(a_ptr, ele%x%etap) .or. &
        associated(a_ptr, ele%y%eta) .or. associated(a_ptr, ele%y%etap) .or. &
        coupling_change) then 
      call make_v_mats (ele, v_mat, v_inv_mat)
      eta_xy_vec = (/ ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap /)
      eta_vec = matmul (v_inv_mat, eta_xy_vec)
      ele%a%eta  = eta_vec(1)
      ele%a%etap = eta_vec(2)
      ele%b%eta  = eta_vec(3)
      ele%b%etap = eta_vec(4)
      return
    endif

    if (associated(a_ptr, ele%a%eta) .or. associated(a_ptr, ele%a%etap) .or. &
        associated(a_ptr, ele%b%eta) .or. associated(a_ptr, ele%b%etap)) then 
      call make_v_mats (ele, v_mat, v_inv_mat)
      eta_vec = (/ ele%a%eta, ele%a%etap, ele%b%eta, ele%b%etap /)
      eta_xy_vec = matmul (v_mat, eta_vec)
      ele%x%eta  = eta_xy_vec(1)
      ele%x%etap = eta_xy_vec(2)
      ele%y%eta  = eta_xy_vec(3)
      ele%y%etap = eta_xy_vec(4)
      return
    endif

  endif

endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_lat_taylors (lattice_in, lattice_out, 
!                                              type_out, transfered_all)
!
! Subroutine to transfer the taylor maps from the elements of one lattice to
! the elements of another. The elements are matched between the lattices so 
! that the appropriate element in lattice_out will get the correct Taylor map
! even if the order of the elements is different in the 2 lattices.
!
! Note: The transfered Taylor map will be truncated to bmad_com%taylor_order.
! Note: If the taylor_order of an element in lattice_in is less than 
!   bmad_com%taylor_order then it will not be used.  
!
! Modules needed:
!   use bmad
!
! Input:
!   lattice_in   -- lat_struct: Input lattice with Taylor maps.
!   type_out  -- Logical: If True then print a message for each Taylor map
!                 transfered.
!
! Output:
!   lattice_out  -- lat_struct: lattice to receive the Taylor maps.
!   transfered_all -- Logical, optional: Set True if a Taylor map is found
!                 for all elements in lattice_out that need one. False otherwise.
!-

subroutine transfer_lat_taylors (lattice_in, lattice_out, type_out, transfered_all)

  implicit none

  type (lat_struct), target, intent(in) :: lattice_in
  type (lat_struct), target, intent(inout) :: lattice_out
  type (ele_struct), pointer :: ele_in, ele_out

  integer i, j
  integer n_in, ix_in(ubound(lattice_in%ele, 1))
 
  logical, intent(in)  :: type_out
  logical, optional :: transfered_all

  character(25) :: r_name = 'transfer_lat_taylors'

! check global parameters

  if (present(transfered_all)) transfered_all = .true.

  if (lattice_in%ele(0)%value(E_tot$) /= &
                              lattice_out%ele(0)%value(E_tot$)) then
    if (type_out) then
       call out_io (s_warn$, r_name, &
              'THE LATTICE ENERGIES ARE DIFFERENT. TAYLOR MAPS NOT TRANSFERED.')
    endif
    if (present(transfered_all)) transfered_all = .false.
    return
  endif

! Find the taylor series in the first lattice.

  n_in = 0
  do i = 1, lattice_in%n_ele_max
    if (associated(lattice_in%ele(i)%taylor(1)%term)) then
      if (bmad_com%taylor_order > lattice_in%ele(i)%taylor_order) cycle
      n_in = n_in + 1
      ix_in(n_in) = i
    endif
  enddo

! Go through lattice_out and match elements.
! If we have a match transfer the Taylor map.
! Call attribute_bookkeeper before transfering the taylor map to make sure
! the check_sum is correct. 

  out_loop: do i = 1, lattice_out%n_ele_max

    ele_out => lattice_out%ele(i)

    do j = 1, n_in

      ele_in => lattice_in%ele(ix_in(j))

      if (equivalent_eles (ele_in, ele_out)) then
        if (type_out) call out_io (s_info$, r_name, &
            ' Reusing Taylor from: ' // trim(ele_in%name) // '  to: ' //  ele_out%name)
        call attribute_bookkeeper (ele_out, lattice_out%param)
        call transfer_ele_taylor (ele_in, ele_out, bmad_com%taylor_order)
        cycle out_loop
      endif

    enddo

    if (ele_out%tracking_method == taylor$ .or. &
                    ele_out%mat6_calc_method == taylor$ .and. type_out) then
      call out_io (s_warn$, r_name, ' NO TAYLOR FOR: ' // ele_out%name)
      if (present(transfered_all)) transfered_all = .false.
    endif

  enddo out_loop

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine set_on_off (key, lat, switch, orb, use_ref_orb)
!
! Subroutine to turn on or off a set of elements (quadrupoles, rfcavities,
! etc.) in a lattice. An element that is turned off acts like a drift.
! lat_make_mat6 will be called to remake lat%ele()%mat6.
!
! Modules needed:
!   use bmad
!
! Input:
!   key          -- Integer: Key name of elements to be turned on or off.
!                      [Key = quadrupole$, etc.]
!   lat             -- lat_struct: lattice structure holding the elements
!   switch       -- Integer: 
!                     on$            => Turn elements on.  
!                     off$           => Turn elements off. 
!                     save_state$    => Save present on/off state. 
!                                         No turning on or off is done.
!                     restore_state$ => Restore saved on/off state.
!   orb(0:)     -- Coord_struct, optional: Needed for lat_make_mat6
!   use_ref_orb -- Logical, optional: If present and true then use the
!                    present ele%ref_orb. Default is false.
!
! Output:
!   lat -- lat_struct: Modified lattice.
!-

#include "CESR_platform.inc"
                                    
subroutine set_on_off (key, lat, switch, orb, use_ref_orb)

  use bmad_struct
  use bmad_interface

  implicit none

  type (lat_struct) lat
  type (coord_struct), optional :: orb(0:)
  type (coord_struct) ref_orb

  integer i, key               
  integer, intent(in) :: switch

  logical, optional :: use_ref_orb
  logical old_state

  character(20) :: r_name = 'set_on_off'

!

  do i = 1, lat%n_ele_max

    if (lat%ele(i)%key /= key) cycle

    old_state = lat%ele(i)%is_on

    select case (switch)
    case (on$) 
      lat%ele(i)%is_on = .true.
    case (off$)
      lat%ele(i)%is_on = .false.
    case (save_state$)
      lat%ele(i)%old_is_on = lat%ele(i)%is_on
      cycle
    case (restore_state$)
      lat%ele(i)%is_on = lat%ele(i)%old_is_on
    case default
      call out_io (s_abort$, r_name, 'BAD SWITCH: \i\ ', switch)
      call err_exit
    end select

    if (old_state .neqv. lat%ele(i)%is_on) then
      if (logic_option (.false., use_ref_orb)) then
        ref_orb%vec = lat%ele(i)%ref_orb_in
        call make_mat6(lat%ele(i), lat%param, ref_orb)
      else
        call lat_make_mat6(lat, i, orb)
      endif
    endif

  enddo

end subroutine

end module
