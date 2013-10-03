!+           
! Subroutine Implement_Pathlength_Patch(path_length_patch,ring,delta_frf, frf) 
!
! Insert a patch element at each RF cavity to adjust path length 
!  due to frequency offset delta_frf
! Input:
!   ring -- lat_struct: New elements will be inserted
!   delta_frf: optional - change in rf frequency (Hz)
!   frf : optional - nominal rf frequency (Hz)
!   patch_length_patch: logical, true if patch elements are already inserted
!
! Output:
!   ring -- lat_struct: ring with new elements inserted
!

  subroutine implement_pathlength_patch(path_length_patch,ring, delta_frf, frf) 

  use bookkeeper_mod

 implicit none

 type (lat_struct) ring
 type (ele_struct), allocatable :: patch_ele(:)

 integer i,j
 integer j_max
 integer ix, ii
 integer, allocatable :: ix_ele(:)

 real(rp), allocatable, save :: s_rf(:) 
 real(rp), optional :: delta_frf, frf
 real(rp) delta_s

 character*16, allocatable, save :: rf_ele(:)
 logical path_length_patch
 character*2 number

! find rf cavities


  j=0
  do i=1,ring%n_ele_track
   if(ring%ele(i)%key == rfcavity$) then
    j=j+1
    if(j == 1) frf = ring%ele(i)%value(rf_frequency$)
   endif
  end do

  j_max = j
!  j_max = 1
  allocate (patch_ele(1:2*j), s_rf(1:2*j), rf_ele(1:2*j), ix_ele(1:2*j))

 if( .not. path_length_patch)then  
  j=0
  number = '00'
  do i=1,ring%n_ele_track
   if(ring%ele(i)%key == rfcavity$) then
    j=j+1
    if(j<10)write(number(2:2),'(i1)')j
    if(j>=10)write(number(1:2),'(i2)')j
    if(j>99) print *,' Too many cavities for implement path length patch'
    rf_ele(j) = ring%ele(i)%name
    call init_ele(patch_ele(j))
    patch_ele(j)%name = 'PATCH_'//number//ring%ele(i)%name(1:len_trim(ring%ele(i)%name))
!    print '(a,i10,1x,a)',' j, patch_ele(j)%name ', j, patch_ele(j)%name
    patch_ele(j)%key = custom$
    patch_ele(j)%s = ring%ele(i-1)%s
    patch_ele(j)%value(custom_attribute1$) = 0. 
    patch_ele(j)%tracking_method = custom$
    patch_ele(j)%mat6_calc_method = tracking$
   endif
  end do
 
  do j = 1,j_max
!    call element_locator(rf_ele(j), ring, ix)
    ix = element_at_s(ring,patch_ele(j)%s, .true.)
!   print *,j,rf_ele(j),ix, ring%ele(ix)%name
    call insert_element(ring,patch_ele(j),ix)
    if (.not. bmad_com%auto_bookkeeper) call lattice_bookkeeper(ring)
    call lat_make_mat6(ring,-1)
!    print *,' ring%n_ele_max, ix = ', ring%n_ele_max,ix
!    print '(i10,1x,2a16)',ix,ring%ele(ix)%name, ring%ele(ix+1)%name

  end do
 endif

  j=0
  do i = 1,ring%n_ele_track
   if(ring%ele(i)%key == custom$ .and. ring%ele(i+1)%key == rfcavity$)then
     j = j+1
     s_rf(j) = ring%ele(i)%s
!     print *,' j, s_rf(j) ', j, s_rf(j)
     ix_ele(j) = i
   endif
  end do

  do j = 1,j_max
    if( j > 1)then
     delta_s = (s_rf(j) - s_rf(j-1))
    else
     delta_s = s_rf(j) + ring%ele(ring%n_ele_track)%s - s_rf(j_max)
    endif
    if(present(frf))then
      ring%ele(ix_ele(j))%value(custom_attribute1$) = delta_s * delta_frf/frf
!     print '(a,3i10,1x,a16)',' j_max,j, ix_ele(j),ring%ele(ix_ele(j))%name ', j_max, j,ix_ele(j),ring%ele(ix_ele(j))%name   
     print '(a,es12.4)',' ring%ele(ix_ele(j))%value(custom_attribute1$) = ',ring%ele(ix_ele(j))%value(custom_attribute1$)

     else
      ring%ele(ix_ele(j))%value(custom_attribute1$) = 0.
    endif
  end do
print '(a,es12.4)',' Change in path length = ', ring%ele(ring%n_ele_track)%s * delta_frf/frf
print '(a,es12.4)',' RF frequency = ', frf
  call lat_make_mat6(ring, -1)
  
  path_length_patch = .true.

  deallocate (patch_ele, s_rf, rf_ele, ix_ele)

  return
  end

!-------------------------------------------------------
