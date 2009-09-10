!+
! Subroutine check_lat_controls (lat, exit_on_error)
!
! Subroutine to check if the control links in a lat structure are valid
!
! Modules needed:
!   use bmad
!
! Input:
!   lat -- lat_struct: Lat to check
!   exit_on_error -- Logical: Exit if an error detected?
!-

#include "CESR_platform.inc"

subroutine check_lat_controls (lat, exit_on_error)

use bmad_struct
use bmad_interface, except_dummy => check_lat_controls

implicit none
     
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, slave, lord
type (branch_struct), pointer :: branch

integer i_t, j, i_t2, ix, s_stat, l_stat, t2_type, n, cc(100), i
integer ix1, ix2, ii, i_b, i_b2

character(24) :: r_name = 'check_lat_controls'

logical exit_on_error, found_err, good_control(12,12)
logical girder_here

! check energy

if (any(lat%ele(:)%key == lcavity$) .and. &
                        lat%param%lattice_type /= linear_lattice$) then
  call out_io (s_fatal$, r_name, &
            'THERE IS A LCAVITY BUT THE LATTICE_TYPE IS NOT SET TO LINEAR_LATTICE!')
endif

! good_control specifies what elements can control what other elements.

good_control = .false.
good_control(group_lord$, (/ free$, overlay_slave$, multipass_slave$ /)) = .true.
good_control(girder_lord$, (/ overlay_slave$, multipass_slave$ /)) = .true.
good_control(overlay_lord$, (/ overlay_slave$, multipass_slave$ /)) = .true.
good_control(super_lord$, (/ super_slave$ /)) = .true.
good_control(multipass_lord$, (/ multipass_slave$ /)) = .true.

found_err = .false.
           
! loop over all elements

do i_b = 0, ubound(lat%branch, 1)

  branch => lat%branch(i_b)

  if (branch%ix_branch /= i_b) then
    call out_io (s_fatal$, r_name, &
              'BRANCH: ' // branch%name, &
              'IS OUT OF ORDER IN THE BRANCH(:) ARRAY \i0\ ', &
              i_array = (/ branch%ix_branch /) )
    found_err = .true.
  endif


  if (i_b > 0) then
    ix = branch%ix_from_branch

    if (ix < 0 .or. ix == i_b .or. ix > ubound(lat%branch, 1)) then 
      call out_io (s_fatal$, r_name, &
                'BRANCH: ' // branch%name, &
                'HAS A IX_FROM_BRANCH INDEX OUT OF RANGE: \i0\ ', i_array = (/ ix /) )
      found_err = .true.
    endif

    if (branch%ix_from_ele < 1 .or. branch%ix_from_ele > lat%branch(ix)%n_ele_track) then
      call out_io (s_fatal$, r_name, &
                'BRANCH: ' // branch%name, &
                'HAS A IX_FROM_ELE INDEX OUT OF RANGE: \i0\ ', i_array = (/ branch%ix_from_ele /) )
      found_err = .true.
    endif

    slave => lat%branch(ix)%ele(branch%ix_from_ele)
    if (slave%key /= branch$ .and. slave%key /= photon_branch$) then
      call out_io (s_fatal$, r_name, &
            'BRANCH: ' // branch%name, &
            'HAS A FROM ELEMENT THAT IS NOT A BRANCH NOR A PHOTON_BRANCH ELEMENT: ' // slave%name)
      found_err = .true.
    endif

  endif

  !--------------------------------

  do i_t = 1, branch%n_ele_max

    ele => branch%ele(i_t)
    l_stat = ele%lord_status
    s_stat = ele%slave_status

    ! branch check

    if (ele%key == branch$ .or. ele%key == photon_branch$) then
      ix = nint(ele%value(ix_branch_to$))
      if (ix < 0 .or. ix > ubound(lat%branch, 1) .or. ix == i_b) then
        call out_io (s_fatal$, r_name, &
                  'ELEMENT: ' // ele%name, &
                  'WHICH IS A: BRANCH OR PHOTON_BRANCH ELEMENT', &
                  'HAS A IX_BRANCH_TO INDEX OUT OF RANGE: \i0\ ', i_array = (/ ix /) )
        found_err = .true.
      endif
    endif

    ! match elements with match_end set should only appear in linear_lattices

    if (ele%key == match$) then
      if (ele%value(match_end$) /= 0 .and. lat%param%lattice_type /= linear_lattice$) then
        call out_io (s_fatal$, r_name, &
                  'ELEMENT: ' // ele%name, &
                  'WHICH IS A: MATCH ELEMENT', &
                  'HAS THE MATCH_END ATTRIBUTE SET BUT THIS IS NOT A LINEAR LATTICE!')
        found_err = .true.
      endif
    endif

    ! sbend multipass lord must have non-zero ref_energy.
    ! Check both p0c and e_tot since if this routine is called by bmad_parser we
    ! can have one zero and the other non-zero.

    if (ele%key == sbend$ .and. l_stat == multipass_lord$) then
      if (ele%value(p0c$) == 0 .and. ele%value(e_tot$) == 0) then
        call out_io (s_fatal$, r_name, &
                  'BEND: ' // ele%name, &
                  'WHICH IS A: MULTIPASS_LORD', &
                  'DOES NOT HAVE A REFERENCE ENERGY DEFINED')
        found_err = .true.
      endif
      if (ele%ref_orbit == match_none$ .and. ele%value(n_ref_pass$) /= 1) then
        call out_io (s_fatal$, r_name, &
                  'BEND: ' // ele%name, &
                  'WHICH IS A: MULTIPASS_LORD', &
                  'HAS REF_ORBIT = MATCH_NONE BUT N_REF_PASS IS NOT SET TO 1!')
        found_err = .true.
      endif
    endif

    ! A multipass lord that is a magnetic or electric element must either:
    !   1) Have field_master = True or
    !   2) Have a defined reference energy.

    if (l_stat == multipass_lord$ .and. .not. ele%field_master .and. ele%value(p0c$) == 0) then
      select case (ele%key)
      case (quadrupole$, sextupole$, octupole$, solenoid$, sol_quad$, sbend$, &
            hkicker$, vkicker$, kicker$, elseparator$, bend_sol_quad$)
        call out_io (s_fatal$, r_name, &
              'FOR MULTIPASS LORD: ' // ele%name, &
              'N_REF_PASS, E_TOT, AND P0C ARE ALL ZERO AND FIELD_MASTER = FALSE!')
        found_err = .true.
      end select
    endif

    ! check that element is in correct part of the ele(:) array

    if (ele%key == null_ele$ .and. i_t > branch%n_ele_track) cycle      

    if (i_t > branch%n_ele_track) then
      if (s_stat == super_slave$) then
        call out_io (s_fatal$, r_name, &
                  'ELEMENT: ' // ele%name, &
                  'WITH SLAVE_STATUS: ' // control_name(s_stat), &
                  'IS *NOT* IN THE TRACKING PART OF LAT LIST AT: \i0\ ', &
                  i_array = (/ i_t /) )
        found_err = .true.
      endif                                             
    else                                                         
      if (l_stat == super_lord$ .or. l_stat == overlay_lord$ .or. &
          l_stat == group_lord$ .or. l_stat == girder_lord$ .or. &
          l_stat == multipass_lord$) then
        call out_io (s_fatal$, r_name, &
                  'ELEMENT: ' // ele%name, &
                  'WITH LORD_STATUS: ' // control_name(l_stat), &
                  'IS IN THE TRACKING PART OF LAT LIST AT: \i0\ ', i_array = (/ i_t /) )
        found_err = .true.
      endif
    endif

    if (.not. any( (/ free$, girder_lord$, super_lord$, overlay_lord$, group_lord$, &
                      multipass_lord$ /) == l_stat)) then
      call out_io (s_fatal$, r_name, &
                'ELEMENT: ' // trim(ele%name) // '  (\i0\)', &
                'HAS UNKNOWN LORD_STATUS INDEX: \i0\ ', i_array = (/ i_t, l_stat /) )
      found_err = .true.
    endif

    if (.not. any( (/ free$, super_slave$, overlay_slave$, multipass_slave$ /) == s_stat)) then
      call out_io (s_fatal$, r_name, &
                'ELEMENT: ' // trim(ele%name) // '  (\i0\)', &
                'HAS UNKNOWN SLAVE_STATUS INDEX: \i0\ ', i_array = (/ i_t, s_stat /) )
      found_err = .true.
    endif

    if (ele%n_slave /= ele%ix2_slave - ele%ix1_slave + 1) then
      call out_io (s_fatal$, r_name, &
                'LORD: ' // trim(ele%name) // '  (\i0\)',  &
                'HAS SLAVE NUMBER MISMATCH: \3i5\ ', &
                i_array = (/ i_t, ele%n_slave, ele%ix1_slave, ele%ix2_slave /) )
      found_err = .true.
      cycle
    endif

    if (ele%n_lord /= ele%ic2_lord - ele%ic1_lord + 1) then
      call out_io (s_fatal$, r_name, &
                'SLAVE: ' // trim(ele%name) // '  (\i0\)', &
                'HAS LORD NUMBER MISMATCH: \3i5\ ', &
                i_array = (/ i_t, ele%n_lord, ele%ic1_lord, ele%ic2_lord /) )
      found_err = .true.
      cycle
    endif

    if (s_stat == overlay_slave$ .and. ele%n_lord == 0) then
      call out_io (s_fatal$, r_name, &
                'OVERLAY_SLAVE: ' // trim(ele%name) // '  (\i0\)', &
                'HAS ZERO LORDS!', i_array = (/ i_t /) )
      found_err = .true.
    endif

    if (s_stat == super_slave$ .and. ele%n_lord == 0) then
      call out_io (s_fatal$, r_name, &
                'SUPER_SLAVE: ' // trim(ele%name) // '  (\i0\)', &
                'HAS ZERO LORDS!', i_array = (/ i_t /) )
      found_err = .true.
    endif

    ! check that super_lord elements have their slaves in the correct order

    if (l_stat == super_lord$) then
      do i = ele%ix1_slave+1, ele%ix2_slave
        ix1 = lat%control(i-1)%ix_slave
        ix2 = lat%control(i)%ix_slave
        if (ix2 > ix1) then
          do ii = ix1+1, ix2-1
            if (lat%ele(ii)%value(l$) /= 0) goto 9000   ! error
          enddo
        elseif (ix2 < ix1) then
          do ii = ix1+1, branch%n_ele_track
            if (lat%ele(ii)%value(l$) /= 0) goto 9000   ! error
          enddo
          do ii = 1, ix2-1            
            if (lat%ele(ii)%value(l$) /= 0) goto 9000   ! error
          enddo
        else
          call out_io (s_fatal$, r_name, &
                    'DUPLICATE SUPER_SLAVES: ', trim(lat%ele(ix1)%name) // '  (\i0)', &
                    'FOR SUPER_LORD: ' // trim(ele%name) // '  (\i0\)', &
                    i_array = (/ ii, i_t /) )
          found_err = .true.
        endif
      enddo
    endif

    ! The slaves of a multipass_lord cannot be controlled by anything else.

    if (l_stat == multipass_lord$) then
      do i = ele%ix1_slave, ele%ix2_slave
        ii = lat%control(i)%ix_slave
        if (lat%ele(ii)%n_lord /= 1) then
          call out_io (s_fatal$, r_name, &
                    'SLAVE OF A MULTIPASS_LORD: ' // trim(lat%ele(ii)%name) // '  (\i0\)', &
                    'HAS MORE THAN ONE LORD.', &
                    'FOR MULTIPASS_LORD: ' // trim(ele%name) // '  (\i0\)', &
                    i_array = (/ ii, i_t /) )
          found_err = .true.
        endif
      enddo
    endif

    ! check slaves

    do j = ele%ix1_slave, ele%ix2_slave

      if (j < 1 .or. j > lat%n_control_max) then
        call out_io (s_fatal$, r_name, &
                  'LORD: ' // trim(ele%name)  // '  (\i0\)', &
                  'HAS IX_SLAVE INDEX OUT OF BOUNDS: \2i5\ ', &
                  i_array = (/ i_t, ele%ix1_slave, ele%ix2_slave /) )
        found_err = .true.
      endif

      if (lat%control(j)%ix_lord /= i_t) then
        call out_io (s_fatal$, r_name, &
                  'LORD: ' // trim(ele%name) // '  (\i0\)', &
                  'HAS A %IX_LORD POINTER MISMATCH: \i0\ ', &
                  'AT: \i0\ ', &
                  i_array = (/ i_t, lat%control(j)%ix_lord, j /) )
        found_err = .true.
      endif

      i_t2 = lat%control(j)%ix_slave
      i_b2 = lat%control(j)%ix_branch

      if (i_t2 < 1 .or. i_t2 > lat%branch(i_b2)%n_ele_max) then
        call out_io (s_fatal$, r_name, &
                  'LORD: ' // trim(ele%name) // '  (\i0\)', &
                  'HAS A SLAVE INDEX OUT OF RANGE: \i0\ ', &
                  'AT: \i0\ ', &
                  i_array = (/ i_t, i_t2, j /) )
        found_err = .true.
        cycle
      endif

      slave => lat%branch(i_b2)%ele(i_t2)  
      t2_type = slave%slave_status      

      if (.not. good_control(l_stat, t2_type) .and. &
                        lat%control(j)%ix_attrib /= l$) then
        call out_io (s_fatal$, r_name, &
                  'LORD: ' // trim(ele%name) // '  (\i0\)',  &
                  'WITH LORD_STATUS: ' // control_name(l_stat), &
                  'HAS A SLAVE: ' // trim(slave%name) // '  (\i0\)', &
                  'WITH SLAVE_STATUS: ' // control_name(t2_type), &
                  i_array = (/ i_t, i_t2 /) )
        found_err = .true.
      endif

      if (l_stat /= group_lord$ .and. l_stat /= girder_lord$) then
        n = slave%ic2_lord - slave%ic1_lord + 1
        cc(1:n) = (/ (lat%ic(i), i = slave%ic1_lord, slave%ic2_lord) /)
        if (.not. any(lat%control(cc(1:n))%ix_lord == i_t)) then
          call out_io (s_fatal$, r_name, &
                    'SLAVE: ', trim(slave%name) // '  (\i0\)', &
                    'WITH SLAVE_STATUS: ' // control_name(t2_type), &
                    'DOES NOT HAVE A POINTER TO ITS LORD: ' // trim(ele%name) // '  (\i0\)', &
                    i_array = (/ i_t2, i_t /) )
          found_err = .true.
        endif
      endif

    enddo

    ! check lords

    girder_here = .false.

    do ix = ele%ic1_lord, ele%ic2_lord

      if (ix < 1 .or. ix > lat%n_control_max) then
        call out_io (s_fatal$, r_name, &
                  'SLAVE: ' // trim(ele%name) // '  (\i0\)', &
                  'HAS IC_LORD INDEX OUT OF BOUNDS: \2i5\ ', &
                  i_array = (/ i_t, ele%ic1_lord, ele%ic2_lord /) )
        found_err = .true.
      endif

      j = lat%ic(ix)

      if (j < 1 .or. j > lat%n_control_max) then
        call out_io (s_fatal$, r_name, &
                  'SLAVE: ' // trim(ele%name) // '  (\i0\)', &
                  'HAS IC INDEX OUT OF BOUNDS: \2i5\ ', & 
                  i_array = (/ i_t, ix, j /) )
        found_err = .true.
      endif
          
      i_t2 = lat%control(j)%ix_lord

      if (i_t2 < 1 .or. i_t2 > branch%n_ele_max) then
        call out_io (s_fatal$, r_name, &
                  'SLAVE: ' // trim(ele%name) // '  (\i0\)', &
                  'HAS A LORD INDEX OUT OF RANGE: \3i5\ ', &
                  i_array = (/ i_t, ix, j, i_t2 /) )
        found_err = .true.
        cycle
      endif

      if (lat%control(j)%ix_slave /= i_t) then
        call out_io (s_fatal$, r_name, &
                  'SLAVE: ' // trim(ele%name) // '  (\i0\)', &
                  'HAS A %IX_SLAVE POINTER MISMATCH: \i0\ ', &
                  'AT: \2i5\ ', &
                  i_array = (/ i_t, lat%control(j)%ix_slave, ix, j /) )
        found_err = .true.
      endif

      lord => lat%ele(i_t2)
      t2_type = lord%lord_status

      if (.not. good_control(t2_type, s_stat)) then
        call out_io (s_fatal$, r_name, &
                  'SLAVE: ' // trim(ele%name) // '  (\i0\)', &
                  'WITH SLAVE_STATUS: ' // control_name(s_stat), &
                  'HAS A LORD: ' // trim(lord%name) // '  (\i0\)', &
                  'WITH LORD_STATUS: ' // control_name(t2_type), &
                  i_array = (/ i_t, i_t2 /) )
        found_err = .true.
      endif

      if (t2_type == girder_lord$) then
        if (girder_here) then
          call out_io (s_fatal$, r_name, &
                    'SLAVE: ' // trim(ele%name) // '  (\i0\)', &
                    'HAS MORE THAN ONE GIRDER_LORD.', &
                    i_array = (/ i_t /) )
          found_err = .true.
        endif
        girder_here = .true.
      endif


    enddo

  enddo

enddo

if (found_err .and. exit_on_error) call err_exit
return

!---------------------------------
! super_lord error

9000 continue

call out_io (s_fatal$, r_name, &
            'SUPER_SLAVES: ' // trim(branch%ele(ix1)%name) // ', ' // branch%ele(ix2)%name, &
            'NOT IN CORRECT ORDER FOR SUPER_LORD: ' // trim(ele%name) // '  (\i0\)', &
            i_array = (/ i_t /) )

if (exit_on_error) call err_exit

end subroutine
