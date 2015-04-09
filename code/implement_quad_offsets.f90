!+           
! Subroutine Implement_quad_offsets(ring,sigy,cutoff,seed) 
!
! Insert a vertical offset at each quad, with some random gaussian distribution 
!
! Input:
!   ring -- lat_struct: Offsets will be added
!   sig_y -- real: width of gaussian
!   cutoff -- real: cutoff of distribution
!   seed -- integer: random number seed
!
! Output:
!   ring -- lat_struct: ring with offsets inserted
!

  subroutine implement_quad_offsets(ring,sigy,cutoff,seed) 
  use bmad
  use random_mod
  
 implicit none

 type (lat_struct) ring

 integer i,j
 integer seed
 
 real(rp) sigy,cutoff
 real(rp), allocatable :: offsets(:)
 real(rp) sig_offset/0./
 
!Initialize random number generator that is used to create phase space distribution and multiple scattering
   call ran_seed_put(seed)
   call ran_seed_get(seed)
 print *,' seed from ran_seed_get =',seed

 
 j=0
  do i=1,ring%n_ele_track
   if(ring%ele(i)%key == quadrupole$)j=j+1
  end do

  allocate (offsets(1:j))
  call ran_gauss(offsets(:))
  j=0
  do i=1,ring%n_ele_track
     if(ring%ele(i)%key == quadrupole$)then
        j=j+1
        ring%ele(i)%value(y_offset$) = max(offsets(j)*sigy, cutoff*sigy)
        sig_offset = sig_offset + ring%ele(i)%value(y_offset$) **2
     endif
  end do
  print '(a,es12.4,1x,a,1x,i10)',' sig_offset = ', sqrt(sig_offset/j), ' Number of quads = ', j

  deallocate (offsets)

  return
  end subroutine
