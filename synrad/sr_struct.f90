module sr_struct

  use bmad_struct
  use bmad_interface
  use precision_def


! The wall is specified by an array of points with straight lines (faces)
! between the points. the face between point i-1 and the point i is
! associated with the structure for point i. If a face represents an opening
! in the wall (eg the openings for the chess beam lines) then that face is
! call a phantom.

  type wall_pt_struct   ! struct for input points
    character(16) name  ! name of element (sliding_joint, etc.)
    real(rp) x, s       ! position of wall point
    integer n_seg       ! how many segments it will be broken up into
    integer ix_seg      ! index to seg(:) array. From ix_seg+1 to ix_seg+n_seg
    integer type        ! no_alley$, outer_wall$, etc.
    integer ix_pt       ! ordered (in s) index to pt() array
    integer closed_end_direct  ! For alleys: direction of the closed end
    logical phantom     ! is the face an opening?
  end type

! For the synch light power computation each face is broken up into segments.
! It is assumed that the power landing upon a segment is constant
! throughout the segment

  type source_struct
    integer ix_ele_source    ! element index at source
    type (coord_struct) start, now    ! 
  end type              

  type sr_power_struct       ! substruct for a segment 
    real(rp) power               ! total power on segment (Watts)
    real(rp) power_per_len       ! power density (Watts / m)
    real(rp) power_per_area      ! power density (Watts / m^2)
    real(rp) photons_per_sec     ! flux hitting segment in photons per sec
    integer ix_ele_source        ! element index for the largest source
    real(rp) s_source            ! s position of the largest source
    integer n_source             ! number of source points
    type (source_struct), pointer :: sources(:)
                                 ! list of source info for rays hitting this seg
  end type

  type wall_seg_struct       ! segment struct
    integer ix_pt            ! index to which point owns this segment
    real(rp) s, x            ! s, x position of the segment at the endpoint
    real(rp) s_mid, x_mid    ! s, x position of the segment at the midpoint
    real(rp) len             ! length of segment
    type (sr_power_struct) sr_power
  end type

  type alley_struct
    integer ix1, ix2, ix3    ! indexs of the three points
    logical where            ! closed_end$, open_end$, inbetween$
  endtype

! a wall is just a collection of points and segments

  type wall_struct           ! either positive_x or negative_x side wall
    integer side             ! positive_x$ or negative_x$ 
    type (wall_pt_struct), allocatable :: pt(:)
    type (wall_seg_struct), allocatable :: seg(:)
    type (alley_struct) :: alley(100)
    integer n_pt_tot
    integer n_seg_tot
    integer n_alley_tot
    integer ix_pt            ! last ix_pt
    real(rp) seg_len_max
  end type

! The computation tracks a set of synchrotron light rays from their source
! to the wall

  type ray_struct       ! struct for a light ray
    integer direction       ! direction of travel, +1 = forward direction
    integer ix_ele          ! index of element we are now tracking through
    integer alley_status
    real(rp) track_len      ! length of the track from the start
    type (coord_struct) start, old, now  ! coords
    logical crossed_end     ! ray crossed the lat end?
    integer ix_source       ! element index at source of ray
    integer ix_wall_pt      ! index of wall point where hit
    real(rp) p1_factor      ! factor for computing the power/length
    real(rp) p2_factor      ! factor for computing the power/area
    real(rp) r_wall         !
    real(rp) g_bend         ! g = |1/rho|  bending radius inverse at source point
    type (twiss_struct) x_twiss, y_twiss    ! twiss at source point
    type (wall_struct), pointer :: wall     ! points to which wall we hit
  end type

! misc stuff

  type ele_power_struct       ! power from a lat element
    real(rp) at_wall              ! power hitting the wall
    real(rp) radiated             ! power radiated
  end type

  type synrad_param_struct
    character(100) lat_file
    real(rp) epsilon_y     ! vertical emit
    real(rp) i_beam        ! beam current
    integer n_slice        ! # of slice per element or wiggler pole
  end type

  type outline_pt_struct
    character(16) name
    character(16) blueprint
    real(rp) s, x
    logical phantom
  end type

  type outline_struct
    character(16) name
    character(16) blueprint
    type (outline_pt_struct) in(200), out(200)
    integer n_out, n_in
    integer ix_out_slide, ix_in_slide
    logical has_alley
    logical zero_is_center
    real(rp) s_center
    logical overlay
  end type

  type wall_list_struct
    character(16) name
    integer ix_outline
    real(rp) s, len
  end type

  type concat_part_struct
    character(16) name
    integer direction
  end type

  type concat_struct
    character(16) name
    type (concat_part_struct) part(50)
    logical overlay
  end type

  character(20) wall_name(-1:1) / 'negative_x_side', '???',  'positive_x_side' /
  integer negative_x$ / -1 /, positive_x$ / 1 /

  integer possible_alley$ / -1 /
  integer no_alley$ / 0 /, inner_wall$ / 1 /, open_end$ / 2 /
  integer middle_wall$ / 3 /, closed_end$ / 4 /, outer_wall$ / 5 /

  integer no_local_alley$ / 0 /, out_of_alley$ / 1 /, in_alley$ / 2 /

end module
