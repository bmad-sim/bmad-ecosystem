!+
! Subroutine LRBBI_crossings(n_bucket, oppos_buckets, cross_positions)
!
! Subroutine to calculate the location of the parasitic crossing points given
!   a bunch and an array of positions of the bunches it will cross.
!
! Modules needed:
!   use bmad
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


  use precision_def

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
                                                                              
