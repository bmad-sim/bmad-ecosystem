module cesr_crossings_mod

  use bmad_struct
  use bmad_interface

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine CESR_CROSSINGS(i_train, j_car, species, n_trains_tot, n_cars, 
!                             cross_positions, n_car_spacing, train_spacing) 
!
! Subroutine to calculate all parasitic crossing points for a given bunch with 
!  every bunch circulating in the other direction. Sets up arrays of bunch 
!  positions and then calls LRBBI_crossings to do the actual computation. The 
!  crossing points are returned in cross_positions, with values measured cw
!  from the interaction point and scaled to be between 0 and 1. 
! 
! Modules needed:
!   use cesr_mod
!
! Input:
!   i_train         --  Integer: Position (index) of train containing the bunch
!   j_car           --  Integer: Position (index) of bunch
!   species         --  Integer: (+1=positrons, -1=electrons)
!   n_trains_tot    --  Integer: Total number of trains
!   n_cars          --  Integer: Number of bunches per train 
!   n_car_spacing(:)   --  Integer, Optional; if not specified subroutine
!                           uses 14ns=7buckets as spacing between bunches.
!                           Otherwise, specify up to ten times between
!                           consecutive bunches, but be sure to stop before
!                           the pattern repeats. I.e., if trains are 
!                           separated by 6ns, 6ns, 8ns, 6ns, 6ns, 8ns, etc.
!                           only put in 3 terms!
!   train_spacing(:)   --  Integer, Optional: if not specified subroutine uses
!                           140, 140, 147 buckets, repeating, between trains.
!                           This works the same way as bunch spacing - only
!                           enter the repeating pattern of time separations
!                           once, and leave the other elements blank or =0.
!                           It is assumed times are measured from the first
!                           bunch of one train to the first bunch of the 
!                           following train.
!
! Output:
!   cross_positions(:) --  Real(rp): array of positions of parasitic 
!                            crossings, where values lie between 0 and 1.
!      
!-

#include "CESR_platform.inc"

subroutine cesr_crossings(i_train, j_car, species, n_trains_tot, n_cars, &
                                 cross_positions, n_car_spacing, train_spacing)

  implicit none

  integer, intent(in) :: i_train, j_car, species, n_trains_tot, n_cars
  integer, optional, intent(in) :: train_spacing(1:10)
  integer, optional, intent(in) ::  n_car_spacing(1:10)
  integer :: i, j, k, tnumber, bnumber, trlength, ierr
  integer :: bunch_tot
  integer :: n_bucket

  real(rp), dimension(:), allocatable :: oppos_buckets
  real(rp), dimension(:), allocatable :: trtimes, btimes
  real(rp), dimension(:), intent(out) :: cross_positions

! Allocate arrays.

  allocate(trtimes(1:max(9,n_trains_tot)), STAT=ierr)
  if (ierr .ne. 0) THEN
    print*, "TRTIMES: ALLOCATION REQUEST DENIED."
    call err_exit
  endif

  bunch_tot = n_trains_tot*n_cars

  allocate(btimes(1:bunch_tot), stat=ierr)
  if (ierr .ne.  0) then
    print*, "TRTIMES: ALLOCATION REQUEST DENIED."
    call err_exit
  endif

  allocate(oppos_buckets(1:bunch_tot), stat=ierr)
  if (ierr .ne. 0) then
    print*, "OPPOS_BUCKETS: ALLOCATION REQUEST DENIED."
    call err_exit
  endif

! If not otherwise specified, use 7 buckets between bunches and 140, 140, 147,
! repeating, buckets between trains.

  btimes = 0
  trtimes = 0


  if (present (train_spacing)) then
    do i = 1, size(train_spacing)
      if (train_spacing(i) == 0) exit
    enddo
    tnumber = i - 1
    trtimes(1:tnumber) = train_spacing(1:tnumber)
  else
    tnumber = 3
    trtimes(1) = 140  !times are in buckets, 1 bucket=2ns
    trtimes(2) = 140
    trtimes(3) = 147
  endif

  if (present (n_car_spacing)) then
    do i = 1, size(n_car_spacing)
      if (n_car_spacing(i) == 0) exit
    enddo
    bnumber = i - 1
      btimes(1:bnumber) = n_car_spacing(1:bnumber)
  else
    bnumber = 1
    btimes(1) = 7
  endif

!Make an array, oppos_buckets, with the bucket number of each consecutive
! bunch.

  i = 1
  do j = 1, bnumber
    do i = 1, (n_cars - 1)/bnumber
      btimes(bnumber*i + j) = btimes(j)
    enddo
  enddo

  btimes(n_cars) = 0

  j = 1
  do i = 1, bnumber
    do j = 2, n_trains_tot
      do k = 1, n_cars-1
        btimes(n_cars*(j-1)+k) = btimes(k)
      enddo
    enddo
  enddo

  trlength = 0
  do i = 1, n_cars-1
    trlength = trlength + btimes(i)
  enddo

  trtimes = trtimes - trlength

  j = 1
  do i = 1, tnumber
    do j = 1, n_trains_tot/tnumber - 1
      trtimes(tnumber * j + i) = trtimes(i)
    enddo
  enddo

  k = 1
  do i = 1, bunch_tot
    if (btimes(i) == 0) then
      btimes(i) = trtimes(k)
      k = k + 1
    endif
  enddo

  oppos_buckets(1) = 0
  do i = 2, bunch_tot
    oppos_buckets(i) = oppos_buckets(i-1) + btimes(i-1)
  enddo

  n_bucket=(i_train-1)*n_cars+j_car

  if (species .eq. +1) then        !These signs work for positrons
    oppos_buckets=-oppos_buckets        !circulating ccw and electrons cw.
    n_bucket=-oppos_buckets(n_bucket)
  else
    n_bucket=-oppos_buckets(n_bucket)
  endif

!Calculate the crossing points of the bunch in n_bucket with the bunches in
!  oppos_bucket.

  call LRBBI_crossings(n_bucket, oppos_buckets, cross_positions)

!Clean up -- deallocate allocatable arrays.

  deallocate(btimes, stat=ierr)
  if (ierr .ne. 0) then
    print*, "BTIMES: DEALLOCATION REQUEST DENIED."
    call err_exit
  endif

  deallocate(trtimes, stat=ierr)
  if (ierr .ne. 0) then
    print*, "TRTIMES: DEALLOCATION REQUEST DENIED."
    call err_exit
  endif

  deallocate(oppos_buckets, stat=ierr)
  if (ierr .ne. 0) then
    print*, "OPPOS_BUCKETS: DEALLOCATION REQUEST DENIED."
    call err_exit
  endif

end subroutine cesr_crossings

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine LRBBI_crossings(n_bucket, oppos_buckets, cross_positions)
!
! Subroutine to calculate the location of the parasitic crossing points given
!   a bunch and an array of positions of the bunches it will cross.
!
! Modules needed:
!   use cesr_mod
!
! Input:
!   n_bucket          -- Real(rp): The bucket (position) of the bunch to find
!                           crossing points for.
!   oppos_buckets(:) -- Integer: Array of buckets (positions) of oppositely
!                           circulating bunches.
!
! Output:
!   cross_positions -- Real(rp): Array of positions of crossing points,
!     measured clockwise from IP:=0, and normalized to be between 0 and 1.
!-

#include "CESR_platform.inc"

subroutine LRBBI_crossings(n_bucket, oppos_buckets, cross_positions)

  implicit none

  real(rp), dimension(:), intent(in) :: oppos_buckets
  real(rp), dimension(:), intent(inout) :: cross_positions
  real(rp) :: bucket_tot

  integer :: bunch_tot
  integer, intent(in) :: n_bucket
  integer :: i, j

!

  bucket_tot = 1281

  bunch_tot = size(oppos_buckets)

  do i = 1, bunch_tot
    cross_positions(i) = mod((oppos_buckets(i) + n_bucket) / 2, bucket_tot)
    do while (cross_positions(i) .lt. 0)
      cross_positions(i) = bucket_tot + cross_positions(i)
    enddo
      cross_positions(i) = mod(cross_positions(i) + bucket_tot / 2, bucket_tot)
  enddo

  do j = bunch_tot + 1, 2 * bunch_tot
    cross_positions(j) = mod(cross_positions(j - bunch_tot) + bucket_tot / 2,&
              bucket_tot)
  enddo

! Normalize values to be between 0 and 1.

  do j = 1, 2 * bunch_tot
    cross_positions(j) = cross_positions(j)/bucket_tot
  enddo

end subroutine LRBBI_crossings

end module
