!+
! Subroutine check_ring_controls (ring, exit_on_error)
!
! Subroutine to check if the control links in a ring structure are valid
!
! Modules needed:
!   use bmad
!
! Input:
!   ring -- Ring_struct: Ring to check
!   exit_on_error -- Logical: Exit if an error detected?
!-

#include "CESR_platform.inc"

subroutine check_ring_controls (ring, exit_on_error)

  use bmad_struct
  use bmad_interface

  implicit none
       
  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele, ele2

  integer i_t, j, i_t2, ix, t_type, t2_type, n, cc(100), i
  integer n_count, ix1, ix2, ii, key

  logical exit_on_error, found_err, good_control(10,10)
  logical i_beam_here

! check energy

  if (any(ring%ele_(:)%key == lcavity$) .and. &
                          ring%param%lattice_type /= linac_lattice$) then
    print *, 'ERROR IN CHECK_RING_CONTROLS: THERE IS A LCAVITY BUT THE'
    print *, '      LATTICE_TYPE IS NOT SET TO LINAC_LATTICE!'
  endif

!

  good_control = .false.
  good_control(group_lord$, (/ group_lord$, overlay_lord$, super_lord$, &
                              i_beam_lord$, free$, overlay_slave$ /)) = .true.
  good_control(i_beam_lord$, (/ super_lord$, overlay_slave$ /)) = .true.
  good_control(overlay_lord$, (/ overlay_lord$, &
                       i_beam_lord$, overlay_slave$, super_lord$ /)) = .true.
  good_control(super_lord$, (/ super_slave$ /)) = .true.

  found_err = .false.
             
! loop over all elements

  do i_t = 1, ring%n_ele_max

    ele => ring%ele_(i_t)
    t_type = ele%control_type

! check that element is in correct part of the ele_(:) array

    if (ele%key == null_ele$ .and. i_t > ring%n_ele_use) cycle      

    if (i_t > ring%n_ele_use) then
      if (t_type == free$ .or. t_type == super_slave$ .or. &
          t_type == overlay_slave$) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: ELEMENT: ', ele%name
        print *, '      WHICH IS A: ', control_name(t_type)
        print *, '      IS *NOT* IN THE REGULAR PART OF RING LIST AT', i_t
        found_err = .true.
      endif                                             
    else                                                         
      if (t_type == super_lord$ .or. t_type == overlay_lord$ .or. &
          t_type == group_lord$ .or. t_type == i_beam_lord$) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: ELEMENT: ', ele%name
        print *, '      WHICH IS A: ', control_name(t_type)
        print *, '      IS IN THE REGULAR PART OF RING LIST AT', i_t
        found_err = .true.
      endif
    endif

    if (.not. any( (/ free$, super_slave$, overlay_slave$, i_beam_lord$, &
                    super_lord$, overlay_lord$, group_lord$ /) == t_type)) then
      print *, 'ERROR IN CHECK_RING_CONTROLS: ELEMENT: ', ele%name
      print *, '      HAS UNKNOWN CONTROL INDEX: ', t_type
      found_err = .true.
    endif

    if (ele%n_slave /= ele%ix2_slave - ele%ix1_slave + 1) then
      print *, 'ERROR IN CHECK_RING_CONTROLS: LORD: ', ele%name, i_t
      print *, '      HAS SLAVE NUMBER MISMATCH:', &
                                  ele%n_slave, ele%ix1_slave, ele%ix2_slave
      found_err = .true.
      cycle
    endif

    if (ele%n_lord /= ele%ic2_lord - ele%ic1_lord + 1) then
      print *, 'ERROR IN CHECK_RING_CONTROLS: SLAVE: ', ele%name, i_t
      print *, '      HAS LORD NUMBER MISMATCH:', &
                                  ele%n_lord, ele%ic1_lord, ele%ic2_lord
      found_err = .true.
      cycle
    endif

    if (t_type == overlay_slave$ .and. ele%n_lord == 0) then
      print *, 'ERROR IN CHECK_RING_CONTROLS: OVERLAY_SLAVE: ', ele%name, i_t
      print *, '      DOES HAS ZERO LORDS'
      found_err = .true.
    endif

    if (t_type == super_slave$ .and. ele%n_lord == 0) then
      print *, 'ERROR IN CHECK_RING_CONTROLS: OVERLAY_SLAVE: ', ele%name, i_t
      print *, '      DOES HAS ZERO LORDS'
      found_err = .true.
    endif

! check that super_lord elements have their slaves in the correct order

    if (t_type == super_lord$) then
      do i = ele%ix1_slave+1, ele%ix2_slave
        ix1 = ring%control_(i-1)%ix_slave
        ix2 = ring%control_(i)%ix_slave
        if (ix2 > ix1) then
          do ii = ix1+1, ix2-1
            if (ring%ele_(ii)%value(l$) /= 0) goto 9000   ! error
          enddo
        elseif (ix2 < ix1) then
          do ii = ix1+1, ring%n_ele_use
            if (ring%ele_(ii)%value(l$) /= 0) goto 9000   ! error
          enddo
          do ii = 1, ix2-1            
            if (ring%ele_(ii)%value(l$) /= 0) goto 9000   ! error
          enddo
        else
          print *, 'ERROR IN CHECK_RING_CONTROLS: DUPLICATE SUPER_SLAVES: ', &
                                                      ring%ele_(ix1)%name, ii
          print *, '      FOR SUPER_SLAVE: ', ele%name, i_t
          found_err = .true.
        endif
      enddo
    endif

! check slaves

    do j = ele%ix1_slave, ele%ix2_slave

      if (j < 1 .or. j > ring%n_control_max) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: LORD: ', ele%name, i_t
        print *, '      HAS IX_SLAVE INDEX OUT OF BOUNDS:', &
                                  ele%ix1_slave, ele%ix2_slave
        found_err = .true.
      endif

      if (ring%control_(j)%ix_lord /= i_t) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: LORD: ', ele%name, i_t
        print *, '      HAS A %IX_LORD POINTER MISMATCH:', &
                                                 ring%control_(j)%ix_lord
        print *, '      AT:', j
        found_err = .true.
      endif

      i_t2 = ring%control_(j)%ix_slave

      if (i_t2 < 1 .or. i_t2 > ring%n_ele_max) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: LORD: ', ele%name, i_t
        print *, '      HAS A SLAVE INDEX OUT OF RANGE:', i_t2
        print *, '      AT:', j
        found_err = .true.
        cycle
      endif

      ele2 => ring%ele_(i_t2)  
      t2_type = ele2%control_type      

      if (.not. good_control(t_type, t2_type) .and. &
                        ring%control_(j)%ix_attrib /= l$) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: LORD: ', ele%name, i_t
        print *, '      WHICH IS A: ', control_name(t_type)
        print *, '      HAS A SLAVE: ', ele2%name, i_t2
        print *, '      WHICH IS A: ', control_name(t2_type)
        found_err = .true.
      endif

      if (t_type /= group_lord$ .and. t_type /= i_beam_lord$) then
        n = ele2%ic2_lord - ele2%ic1_lord + 1
        cc(1:n) = (/ (ring%ic_(i), i = ele2%ic1_lord, ele2%ic2_lord) /)
        if (.not. any(ring%control_(cc(1:n))%ix_lord == i_t)) then
          print *, 'ERROR IN CHECK_RING_CONTROLS: SLAVE: ', ele2%name, i_t2
          print *, '      WHICH IS A: ', control_name(t2_type)
          print *, '      DOES NOT HAVE A POINTER TO ITS LORD: ', ele%name, i_t
          found_err = .true.
        endif
      endif

    enddo


    if (t_type == super_lord$) then
      do j = ele%ix1_slave, ele%ix2_slave
        i_t2 = ring%control_(j)%ix_slave
        ele2 => ring%ele_(i_t2)  
        if (ele2%key == sol_quad$) then
          if (ele%key == solenoid$) cycle
          if (ele%key == quadrupole$) cycle
        endif
        if (ele2%key == ele%key) cycle
        print *, 'ERROR IN CHECK_RING_CONTROLS: SUPER_LORD: ', ele%name
        print *, '      WHICH IS A: ', key_name(ele%key)
        print *, '      CANNOT HAVE A SUPER_SLAVE: ', ele2%name
        print *, '      WHICH IS A: ', key_name(ele2%key)
        found_err = .true.
      enddo
    endif    

! check lords

    i_beam_here = .false.

    do ix = ele%ic1_lord, ele%ic2_lord

      if (ix < 1 .or. ix > ring%n_control_max) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: SLAVE: ', ele%name, i_t
        print *, '      HAS IC_LORD INDEX OUT OF BOUNDS:', &
                                  ele%ic1_lord, ele%ic2_lord
        found_err = .true.
      endif

      j = ring%ic_(ix)

      if (j < 1 .or. j > ring%n_control_max) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: SLAVE: ', ele%name, i_t
        print *, '      HAS IC_ INDEX OUT OF BOUNDS:', ix, j
        found_err = .true.
      endif
          
      i_t2 = ring%control_(j)%ix_lord

      if (i_t2 < 1 .or. i_t2 > ring%n_ele_max) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: SLAVE: ', ele%name, i_t
        print *, '      HAS A LORD INDEX OUT OF RANGE:', ix, j, i_t2
        found_err = .true.
        cycle
      endif

      if (ring%control_(j)%ix_slave /= i_t) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: SLAVE: ', ele%name, i_t
        print *, '      HAS A %IX_SLAVE POINTER MISMATCH:', &
                                                 ring%control_(j)%ix_slave
          print *, '      AT:', ix, j
        found_err = .true.
      endif

      ele2 => ring%ele_(i_t2)
      t2_type = ele2%control_type

      if (.not. good_control(t2_type, t_type)) then
        print *, 'ERROR IN CHECK_RING_CONTROLS: SLAVE: ', ele%name, i_t
        print *, '      WHICH IS A: ', control_name(t_type)
        print *, '      HAS A LORD: ', ele2%name, i_t2
        print *, '      WHICH IS A: ', control_name(t2_type)
        found_err = .true.
      endif

      if (t2_type == i_beam_lord$) then
        if (i_beam_here) then
          print *, 'ERROR IN CHECK_RING_CONTROLS: SLAVE: ', ele%name, i_t
          print *, '      HAS MORE THAN ONE I_BEAM_LORD.'
          found_err = .true.
        endif
        i_beam_here = .true.
      endif


    enddo

  enddo

  if (found_err .and. exit_on_error) call err_exit
  return

!---------------------------------
! super_lord error

9000 continue

  print *, 'ERROR IN CHECK_RING_CONTROLS: SUPER_SLAVES: ', &
                                ring%ele_(ix1)%name, ring%ele_(ix2)%name
  print *, '      NOT IN CORRECT ORDER FOR SUPER_LORD: ', ele%name, i_t

  if (exit_on_error) call err_exit

end subroutine
