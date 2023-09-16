!+
! Module aperture_by_tracking_mod
!
! This module contains code for the aperture_by_tracking program.  The code was moved here
! from aperture_by_tracking.f90 to make pzap_by_tracking.f90 easier to read.
!-
MODULE aperture_by_tracking_mod

USE precision_def
USE bmad !needed for err_exit
USE slice_mod

IMPLICIT NONE

CONTAINS

!+
! Subroutine progress_indicator(i,max)
!
! Where i is the current iteration and max it the total number of iterations to do,
! this subroutine prints out the percent done in integer increments.
!
! Input:
!   i    -- INTEGER, INTENT(IN): current iteration
!   max  -- INTEGER, INTENT(IN): total number of iterations to do
! Output:
!   None
!-
SUBROUTINE progress_indicator(i,max)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: max
  INTEGER pct_done
  INTEGER, SAVE :: last_print = -1

  pct_done = FLOOR((1.0*i)/(1.0*max)*100.0)

  IF(pct_done .gt. last_print) THEN
    WRITE(*,'(A,I6,A)') "Progress: ", pct_done, "%"
    last_print = pct_done
  ENDIF
END SUBROUTINE progress_indicator

!+
! Subroutine prep_lat_ring(lat,co,Qx,Qy,Qz,master)
!-
SUBROUTINE prep_lat_ring(lat,co,Qx,Qy,Qz,master)

  IMPLICIT NONE

  TYPE(lat_struct) lat
  TYPE(coord_struct), ALLOCATABLE :: co(:)
  type (ele_pointer_struct), allocatable :: eles(:)
  REAL(rp), ALLOCATABLE :: dk1(:)
  REAL(rp) Qx, Qy, Qz
  LOGICAL master
  LOGICAL ok

  CALL set_on_off(rfcavity$, lat, on$)
  bmad_com%radiation_damping_on = .true.
  CALL closed_orbit_calc(lat,co,6)
  call lat_make_mat6(lat, -1, co)
  CALL twiss_at_start(lat)
  CALL twiss_propagate_all(lat)

  !Set lattice tunes
  IF( (Qx .gt. 0.0) .or. (Qy .gt. 0.0) ) THEN
    IF(master) THEN
      WRITE(*,*) "before tuning:"
      WRITE(*,'(A,F11.3)') "  a tune: ", lat%ele(lat%n_ele_track)%a%phi/twopi
      WRITE(*,'(A,F11.3)') "  b tune: ", lat%ele(lat%n_ele_track)%b%phi/twopi
    ENDIF
    CALL set_on_off(rfcavity$, lat, off$)
    ALLOCATE(dk1(lat%n_ele_max))
    CALL choose_quads_for_set_tune(lat%branch(0), dk1, eles)
    ok = set_tune(Qx*twopi, Qy*twopi, dk1, eles, lat%branch(0), co)
    DEALLOCATE(dk1)
    CALL set_on_off(rfcavity$, lat, on$)
    IF(master) THEN
      WRITE(*,*) "after tuning:"
      WRITE(*,'(A,F11.3)') "  a tune: ", lat%ele(lat%n_ele_track)%a%phi/twopi
      WRITE(*,'(A,F11.3)') "  b tune: ", lat%ele(lat%n_ele_track)%b%phi/twopi
    ENDIF
  ENDIF

  IF( Qz .gt. -998.0 ) THEN
    IF(master) THEN
      WRITE(*,*) "before tuning:"
      CALL calc_z_tune(lat%branch(0))
      WRITE(*,'(A,F11.3)') "  z tune: ", lat%z%tune/twopi
    ENDIF
    CALL set_z_tune(lat%branch(0), Qz*twopi, ok)
    IF(master) THEN
      WRITE(*,*) "after tuning:"
      WRITE(*,'(A,F11.3)') "  z tune: ", lat%z%tune/twopi
    ENDIF
  ENDIF
END SUBROUTINE prep_lat_ring

SUBROUTINE make_slices(lat,slice_method,end_s,start_s,slice_length,slices,n_slices)
  IMPLICIT NONE

  TYPE(lat_struct) lat
  type (ele_struct), pointer :: ele
  CHARACTER slice_method*6
  REAL(rp) s, end_s, start_s, slice_length
  REAL(rp), ALLOCATABLE :: temp_slices(:), slices(:)
  INTEGER n_slices, n_skipped, ix_ele

  INTEGER i, j, n
  REAL(rp) actual_slice_length

  IF(slice_method .eq. 'bystep') THEN
    !Slice by distance.  Each slice is the same length.
    !slices(:) is setup so that slice(n_slices)
    !is lat%param%total_length-actual_slice_length.
    n_slices=CEILING((end_s-start_s)/slice_length)
    ALLOCATE(temp_slices(n_slices))
    actual_slice_length = (end_s-start_s)/n_slices
    n_skipped=0
    s = start_s
    n = 0
    DO i=1,n_slices
      s = start_s + actual_slice_length*(i-1)
      ix_ele = element_at_s(lat, s, .false.)
      ! Skip unslicable elements
      if (lat%ele(ix_ele)%key == taylor$) cycle
      ! Add to slice list
      n = n + 1
      temp_slices(n) = s
    ENDDO
    n_slices = n
    allocate(slices(n_slices))
    slices(1:n) = temp_slices(1:n)
    ! Make sure the last slice is within the lattice
    slices(n_slices) = MIN(temp_slices(n_slices),lat%param%total_length)
    ! Cleanup
    deallocate(temp_slices)
  ELSEIF(slice_method .eq. 'byelem') THEN
    !Slices coincide with element locations.
    n_slices = 1
    !Count slices
    DO i=1,lat%n_ele_track
      if (lat%ele(i)%key == taylor$) cycle
      if( lat%ele(i)%s .gt. end_s) exit
      IF(lat%ele(i)%value(l$) .gt. 0) THEN
        n_slices=n_slices+1
      ENDIF
    ENDDO
    ALLOCATE(slices(1:n_slices))
    j=1
    slices(j) = 0.0

    !Populate slices
    DO i=1,lat%n_ele_track-1
      if (lat%ele(i)%key == taylor$) cycle
      if( lat%ele(i)%s .gt. end_s) exit
      IF(lat%ele(i)%value(l$) .gt. 0.) THEN
        j=j+1
        slices(j) = lat%ele(i)%s
      ENDIF
    ENDDO
    slices(n_slices) = end_s
  ELSE
    WRITE(*,*) "FATAL: unknown slice_method"
    STOP
  ENDIF
END SUBROUTINE make_slices

!+
!-
SUBROUTINE check_if_lost_ring(lat,start_s,vec_start,nturns, track_state)
  IMPLICIT NONE

  TYPE(lat_struct) lat
  TYPE(coord_struct), ALLOCATABLE, SAVE :: orb(:)
  INTEGER nturns
  INTEGER i
  REAL(rp) start_s
  real(rp) vec_start(6)
  TYPE(coord_struct) coord_start
  TYPE(coord_struct) coord_end

  TYPE(ele_pointer_struct), ALLOCATABLE :: eles(:)
  INTEGER n_loc
  LOGICAL err
  REAL(rp) freq, half_period
  INTEGER track_state

  CALL lat_ele_locator('rfcavity::*', lat, eles, n_loc, err)
  freq = eles(1)%ele%value(rf_frequency$)
  half_period = c_light/freq/2.0

  IF( .not.ALLOCATED(orb) ) ALLOCATE(orb(0:lat%n_ele_track))

  !First track from start_s to the last element of the ring
  call init_coord(coord_start, vec_start, lat%ele(element_at_s(lat,start_s,.true.)), element_end=upstream_end$)
  CALL track_from_s_to_s(lat, start_s, lat%param%total_length, coord_start, coord_end, track_state=track_state)
  IF(track_state .eq. moving_forward$) THEN
    !particle was not lost.
    !Track over many turns
    orb(0) = coord_end
    DO i=1,nturns
      CALL track_all(lat,orb,track_state=track_state)

      IF(track_state .ne. moving_forward$) THEN
        !particle was lost
        coord_end = orb(track_state)
        coord_end%state = lost$
        EXIT
      ENDIF

      IF(ABS(orb(lat%n_ele_track)%vec(5)) .gt. half_period) THEN 
        !particle is outside RF bucket 
        coord_end = orb(lat%n_ele_track) 
        coord_end%state = lost$ 
        track_state = lat%n_ele_track 
        EXIT 
      ENDIF 
      orb(0) = orb(lat%n_ele_track)
    ENDDO
  ENDIF
END SUBROUTINE check_if_lost_ring

!+
!-
SUBROUTINE check_if_lost_linac(lat,start_s,vec_start,track_state,track_till,halo,halo_aperture,halo_h_emittance)
  IMPLICIT NONE

  TYPE(lat_struct) lat
  REAL(rp) start_s
  real(rp) vec_start(6)
  TYPE(coord_struct) coord_start
  TYPE(coord_struct) coord_end
  INTEGER track_state
  REAL(rp) track_till
  LOGICAL halo 
  REAL(rp) halo_aperture
  REAL(rp) halo_h_emittance

  TYPE(ele_struct) ele_at_end
  REAL(rp) pt_x, pt_xp
  REAL(rp) psi, pt_r2, e_r2
  REAL(rp) betaF, alphaF, gammaF

  call init_coord(coord_start, vec_start, lat%ele(element_at_s(lat,start_s,.true.)), element_end=upstream_end$)

  CALL track_from_s_to_s(lat, start_s, track_till, coord_start, coord_end, track_state=track_state)
  IF(halo) THEN
    IF(coord_end%state .eq. alive$) THEN
      CALL twiss_and_track_at_s(lat,track_till,ele_at_end)
      alphaF = ele_at_end%a%alpha
      betaF =  ele_at_end%a%beta
      gammaF = ele_at_end%value(E_TOT$)/(0.511E6)
      !check if particle lays outside n-sigma of phase space at last ele, where
      !n = halo_aperture
      !see ehrlichman lab book 2 page 60
      pt_x = coord_end%vec(1)
      pt_xp = coord_end%vec(2)
      IF((betaF*pt_xp + alphaF*pt_x) .ne. 0.0_rp) THEN
        psi = ATAN(pt_x/(betaF*pt_xp + alphaF*pt_x))
        pt_r2 = pt_x*pt_x + pt_xp*pt_xp
        e_r2 = (halo_aperture**2)*halo_h_emittance/gammaF*(betaF*SIN(psi)**2 + 1./betaF*(COS(psi)-alphaF*SIN(psi))**2)
        IF(pt_r2 .gt. e_r2) THEN
          coord_end%state = lost$
        ENDIF
      ENDIF
    ENDIF
  ENDIF
END SUBROUTINE check_if_lost_linac

!+
!-
SUBROUTINE binary_search(lost,delta_m,accuracy,reset) !updates delta_m according to one iteration of binary search
  LOGICAL lost
  REAL(rp) delta_m
  REAL(rp) accuracy
  LOGICAL reset

  LOGICAL, SAVE :: high_found
  LOGICAL, SAVE :: low_found
  LOGICAL, SAVE :: bracket_found
  REAL(rp), SAVE :: last_high
  REAL(rp), SAVE :: last_low

  IF(reset) THEN
    high_found = .false.
    low_found = .false.
    bracket_found = .false.
    accuracy = 10.0
  ELSE
    ! If an upper and lower bound for the aperture have have not yet been found, double or half delta_m
    IF(.not. bracket_found) THEN
      IF(lost) THEN
        high_found = .true.
        last_high = delta_m
        IF(.not. low_found) THEN
          delta_m = delta_m / 2.0_rp
        ENDIF
      ELSE
        low_found = .true.
        last_low = delta_m
        IF(.not. high_found) THEN
          delta_m = delta_m * 2.0_rp
        ENDIF
      ENDIF
      bracket_found = (high_found .and. low_found)
      accuracy = 10.0
    ENDIF

    ! If and upper and lower bound have been found, move half the distance
    IF(bracket_found) THEN
      IF(lost) THEN
        ! Decrease delta_m
        last_high = delta_m
        delta_m = delta_m - ABS(delta_m-last_low)/2.0_rp
      ELSE
        ! Increase delta_m
        last_low = delta_m
        delta_m = delta_m + ABS(last_high-delta_m)/2.0_rp
      ENDIF
      accuracy = ABS(last_high-last_low)/MIN(ABS(last_high),ABS(last_low))
    ENDIF

    ! Check for extreme values.  If the aperture is larger then 100%, record 100% and move on.
    IF(ABS(delta_m) .gt. 1.0) THEN
      accuracy = 0.0
      !The next two lines result in 1.0 being recorded as the aperture for this location.
      last_high = 1.1_rp
      last_low = 0.9_rp
    ENDIF
  ENDIF
END SUBROUTINE binary_search

END MODULE aperture_by_tracking_mod






