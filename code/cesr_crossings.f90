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
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   i_train         --  Integer: Position (index) of train containing the bunch
!   j_car           --  Integer: Position (index) of bunch
!   species         --  Integer: (+1=positrons, -1=electrons)
!   n_trains_tot    --  Integer: Total number of trains
!   n_cars          --  Integer: Number of bunches per train 
!   n_car_spacing(:)   --  Integer, Optional; if not specified subroutine
!       															uses 14ns=7buckets as spacing between bunches.
!       															Otherwise, specify up to ten times between
!       															consecutive bunches, but be sure to stop before
!       															the pattern repeats. I.e., if trains are 
!       															separated by 6ns, 6ns, 8ns, 6ns, 6ns, 8ns, etc.
!       															only put in 3 terms!
!   train_spacing(:)   --  Integer, Optional: if not specified subroutine uses
!       															140, 140, 147 buckets, repeating, between trains.
!       															This works the same way as bunch spacing - only
!       															enter the repeating pattern of time separations
!       															once, and leave the other elements blank or =0.
!       															It is assumed times are measured from the first
!       															bunch of one train to the first bunch of the 
!       															following train.
!
! Output:
!   cross_positions(:) --  Real: array of positions of parasitic 
!																			crossings, where values lie between 0 and 1.
!      
!-

!$Id$
!$Log$
!Revision 1.3  2002/01/08 21:44:37  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:49  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine cesr_crossings(i_train, j_car, species, n_trains_tot, n_cars, &
                                 cross_positions, n_car_spacing, train_spacing)

  use bmad_struct
  use bmad_interface

  implicit none

  integer, intent(in) :: i_train, j_car, species, n_trains_tot, n_cars         
  integer, optional, intent(in) :: train_spacing(1:10)
	integer, optional, intent(in) ::	n_car_spacing(1:10)
  integer :: i, j, k, p, pp, tnumber, bnumber, trlength, ierr            
  integer :: bunch_tot                                  
  real :: n_bucket
                                                                               
  real :: length                                                               
  real, dimension(:), allocatable :: oppos_buckets                
  real, dimension(:), allocatable :: trtimes, btimes                           
  real, dimension(:), intent(out) :: cross_positions            
                                                                               
! Allocate arrays.

  allocate(trtimes(1:n_trains_tot), STAT=ierr)            
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

