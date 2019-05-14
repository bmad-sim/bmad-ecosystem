!+
! Module slice_mod
!
! Contains a routine to assist in partial element tracking.
!
! track_s_to_s is passed the slices array and two integers, slix_start and slix_stop.
! It tracks from slices(slix_start) to slices(slix_stop).  It returns the orbit
! at locations in the slices array, and also information if the particle is lost.
!-
MODULE slice_mod

USE bmad

PUBLIC track_s_to_s

CONTAINS

!+
! Subroutine track_s_to_s(lat, slices, slix_start, slix_stop, slorbit, slix_lost, lost_to_col)
!
! For description of this subroutine, see slice_mod description.
!
! Input:
!   lat                 -- TYPE(lat_struct), INTENT(INOUT): lattice
!   slices(:)           -- REAL(rp), INTENT(IN): array of lattice locations
!   slorbit(slix_start) -- TYPE(coord_struct): initial coordinates of particle
!   slix_start          -- INTEGER, INTENT(IN): tracking starts at slices(slix_start)
!   slix_stop           -- INTEGER, INTENT(IN): tracking stops at slices(slix_stop)
! Output:
!   slorbit(1:)  -- TYPE(coord_struct), INTENT(INOUT): orbit at locations corresponding to slices array
!   slix_lost    -- INTEGER, OPTIONAL, INTENT(OUT): if particle lost, set to first slix after where particle was lost
!   lost_to_col  -- INTEGER, OPTIONAL, INTENT(OUT): if particle lost to collimator, set to collimator element number
!-
SUBROUTINE track_s_to_s(lat, slices, slix_start, slix_stop, slorbit, slix_lost, lost_to_col, plane_lost_at)
  USE bmad

  IMPLICIT NONE
  
  TYPE(lat_struct), INTENT(INOUT) :: lat
  REAL(rp), INTENT(IN) :: slices(:)
  INTEGER, INTENT(IN) :: slix_start, slix_stop
  TYPE(coord_struct), INTENT(INOUT) :: slorbit(1:)
  INTEGER, OPTIONAL :: slix_lost
  INTEGER, OPTIONAL :: lost_to_col
  INTEGER, OPTIONAL :: plane_lost_at

  INTEGER i
  INTEGER track_state
  INTEGER last_slice
  REAL(rp) slice_b

  DO i=slix_start, slix_stop-1
    IF(slices(i+1) .ge. lat%param%total_length) THEN
      slice_b = lat%param%total_length - 0.0001_rp
    ELSE
      slice_b = slices(i+1)
    ENDIF
    CALL track_from_s_to_s(lat, slices(i), slice_b, slorbit(i), slorbit(i+1),track_state=track_state)
    last_slice = i+1
    IF(track_state .ne. moving_forward$) THEN
      EXIT
    ENDIF
  ENDDO

  IF(PRESENT(lost_to_col)) THEN
    lost_to_col = -1
    IF(track_state .ne. moving_forward$) THEN
      IF( (lat%ele(track_state)%key .eq. ecollimator$) .or. \
          (lat%ele(track_state)%key .eq. rcollimator$) ) THEN
        lost_to_col = track_state
      ENDIF
    ENDIF
  ENDIF

  IF(PRESENT(slix_lost)) THEN
    IF(track_state .ne. moving_forward$) THEN
      !particle is lost
      slix_lost = i+1
    ELSE
      !particle is not lost
      slix_lost = -1
    ENDIF
  ENDIF

  IF(PRESENT(plane_lost_at)) THEN
    IF(track_state .ne. moving_forward$) THEN
      IF( (slorbit(last_slice)%state .eq. lost_neg_x_aperture$).or. &
          (slorbit(last_slice)%state .eq. lost_pos_x_aperture$) ) THEN
        plane_lost_at = x_plane$
      ELSEIF( (slorbit(last_slice)%state .eq. lost_neg_y_aperture$).or. &
              (slorbit(last_slice)%state .eq. lost_pos_y_aperture$) ) THEN
        plane_lost_at = y_plane$
      ELSEIF( slorbit(last_slice)%state .eq. lost_pz_aperture$ ) THEN
        plane_lost_at = z_plane$
      ELSE
        plane_lost_at = 0  !unknown
      ENDIF
    ENDIF
  ENDIF

END SUBROUTINE

END MODULE slice_mod
