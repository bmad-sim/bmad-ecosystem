module cesr_crossings_mod

  use bmad_struct
  use bmad_interface

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine CESR_CROSSINGS(i_train, j_car, species, n_trains_tot, n_cars, 
!                           cross_positions, ptrain, pcar, &
!                           n_car_spacing, train_spacing )
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
!   n_car_spacing(1:10)   --  Integer, Optional; if not specified subroutine
!                          uses 14ns=7buckets as spacing between bunches.
!                          Otherwise, specify up to ten times between
!                          consecutive bunches, but be sure to stop before
!                          the pattern repeats. I.e., if trains are 
!                          separated by 6ns, 6ns, 8ns, 6ns, 6ns, 8ns, etc.
!                          only put in 3 terms!
!   train_spacing(1:10)   --  Integer, Optional: if not specified subroutine uses
!                          140, 140, 147 buckets, repeating, between trains.
!                          This works the same way as bunch spacing - only
!                          enter the repeating pattern of time separations
!                          once, and leave the other elements blank or =0.
!                          It is assumed times are measured from the first
!                          bunch of one train to the first bunch of the 
!                          following train.
!
! Output:
!   cross_positions(:) --  Real(rp): array of positions of parasitic 
!                          crossings, where values are >=0 and < 1,
!                          with both 0 and 1 being the south IP.
!   ptrain(1:900)       --  Integer: Train associated with this crossing.
!   pcar(1:900)         --  Integer: Car associated with this crossing.
!      
!-

#include "CESR_platform.inc"

subroutine cesr_crossings(i_train, j_car, species, n_trains_tot, n_cars, &
                          cross_positions, ptrain, pcar, &
                          n_car_spacing, train_spacing )

  implicit none

  integer, intent(in) :: i_train, j_car, species, n_trains_tot, n_cars         
  integer, optional, intent(in) :: train_spacing(1:10)
  integer, optional, intent(in) :: n_car_spacing(1:10)
  integer, optional, intent(out) :: ptrain(1:900)
  integer, optional, intent(out) :: pcar(1:900)
  integer :: i, j, k, p, pp, tnumber, bnumber, trlength, ierr            
  integer :: bunch_tot                                  
  real(rp) :: n_bucket
  
  integer :: ib,it                                                                             
  real(rp) :: length                                                               
  real(rp), dimension(:), allocatable :: oppos_buckets                
  real(rp), dimension(:), allocatable :: trtimes, btimes                           
  real(rp), dimension(:), intent(out) :: cross_positions            

!  integer ptrain(1:90)            

!  integer pcar(1:90)
                                                                               
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

     write(6,1100)n_car_spacing,train_spacing
1100 format(' CESR_CROSSINGS called with n_car_spacing=',10i5,'   train_spacing=',10i5)

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

  if ( n_bucket .gt. bunch_tot ) then
    print*, "OPPOS_BUCKETS array size exceeded. N_BUCKET,BUNCH_TOT= ", n_bucket, bunch_tot
    call err_exit
  endif

  if (species .eq. +1) then        !These signs work for positrons 
    oppos_buckets=-oppos_buckets        !circulating ccw and electrons cw.
    n_bucket=-oppos_buckets(n_bucket)
  else 
    n_bucket=-oppos_buckets(n_bucket)
  endif

!Calculate the crossing points of the bunch in n_bucket with the bunches in
!  oppos_bucket. 

  call LRBBI_crossings(n_bucket, oppos_buckets, cross_positions)

! Calculate train and car associated with each p.c.  20 aug 04 jac
! 2562*cross_positions looks like this for 9x5 e+ and t1.b1 e-:
!
! 1281.0 1274.0 1267.0 1260.0 1253.0
! 1141.0 1134.0 1127.0 1120.0 1113.0
! 1001.0  994.0  987.0  980.0  973.0
!  854.0  847.0  840.0  833.0  826.0
!  714.0  707.0  700.0  693.0  686.0
!  574.0  567.0  560.0  553.0  546.0
!  427.0  420.0  413.0  406.0  399.0
!  287.0  280.0  273.0  266.0  259.0
!  147.0  140.0  133.0  126.0  119.0
!    0.0 2555.0 2548.0 2541.0 2534.0
! 2422.0 2415.0 2408.0 2401.0 2394.0
! 2282.0 2275.0 2268.0 2261.0 2254.0
! 2135.0 2128.0 2121.0 2114.0 2107.0
! 1995.0 1988.0 1981.0 1974.0 1967.0
! 1855.0 1848.0 1841.0 1834.0 1827.0
! 1708.0 1701.0 1694.0 1687.0 1680.0
! 1568.0 1561.0 1554.0 1547.0 1540.0
! 1428.0 1421.0 1414.0 1407.0 1400.0
!
! So each e+ bunch has two parasitic crossings.
! So loop over each half of the array separately.

  if ( present ( pcar ) .and. present (ptrain) ) then
   do j = 1, bunch_tot + 1, bunch_tot
    do i = j, j+bunch_tot-1
     ib = mod ( i - j, n_cars ) + 1
     it = mod ( ( i - j ) / n_cars, n_trains_tot ) + 1
     if ( i .gt. 900 ) then
       print *,' CESR_CROSSINGS: WARNING. PCAR array size limit of 900 exceeded. i=',i
       cycle
     endif
     pcar ( i ) = ib
     ptrain ( i ) = it
!     write(6,2000)i,it,ib,2562*cross_positions(i)
2000 format(' CESR_CROSSINGS: i, train, bunch, cross_position= ',3i6,900f8.1)
    enddo
   enddo
  endif
                                                                           
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
!   n_bucket      		-- Real(rp): The bucket (position) of the bunch to find
!                      		 crossing points for.
!   oppos_buckets(:) -- Integer: Array of buckets (positions) of oppositely
!                      		 circulating bunches.
!
! Output:
!   cross_positions -- Real(rp): Array of positions of crossing points,
!   	measured clockwise from IP:=0, and normalized to be between 0 and 1.
!-

#include "CESR_platform.inc"

subroutine LRBBI_crossings(n_bucket, oppos_buckets, cross_positions)


  implicit none
                                                                              
  integer :: bunch_tot
  real(rp), intent(in) :: n_bucket
  real(rp), dimension(:), intent(in) :: oppos_buckets                        
  real(rp), dimension(:), intent(inout) :: cross_positions                      
  integer :: i, j
  real(rp) :: bucket_tot
                                                                        
!

  bucket_tot = 1281

	bunch_tot = size(oppos_buckets)

!     print *,' LRBBI_CROSSINGS: n_bucket= ',n_bucket
                                                                           
  do i = 1, bunch_tot
    cross_positions(i) = mod((oppos_buckets(i) + n_bucket) / 2, bucket_tot)
    do while (cross_positions(i) .lt. 0)                                    
      cross_positions(i) = bucket_tot + cross_positions(i)
    enddo                                                                   
      cross_positions(i) = mod(cross_positions(i) + bucket_tot / 2, bucket_tot)
!     print *,' LRBBI_CROSSINGS: i, cross_position= ',i,cross_positions(i)
  enddo                                                                     
                                                                              
  do j = bunch_tot + 1, 2 * bunch_tot
    cross_positions(j) = mod(cross_positions(j - bunch_tot) + bucket_tot / 2,& 
              bucket_tot)   
!     print *,' LRBBI_CROSSINGS: j, cross_position= ',i,cross_positions(j)
  enddo                                                                   
                                                                              
! Normalize values to be between 0 and 1.
                                                                              
  do j = 1, 2 * bunch_tot                                                       
    cross_positions(j) = cross_positions(j)/bucket_tot                     
  enddo                                                                     

  print *
  print *,'Rtn LRBBI_crossings calculates cross_positions array:'
  write(6,'(5f7.1)')cross_positions*bucket_tot*2
!  print *
!  print *,'Rtn LRBBI_crossings uses array oppos_buckets:'
!  write(6,'(5f7.1)')oppos_buckets
                                                                              
end subroutine LRBBI_crossings                                            
                                                                              


end module
