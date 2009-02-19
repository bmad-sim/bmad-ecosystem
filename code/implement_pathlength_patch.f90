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

  use bmad_struct
  use bmad_interface

 implicit none

 type (lat_struct) ring
 type (ele_struct), allocatable :: patch(:)

 integer i,j
 integer j_max
 integer ix, ii
 integer, allocatable :: ix_ele(:)

 real(rp), allocatable, save :: s_rf(:) 
 real(rp), optional :: delta_frf, frf
 real(rp) delta_s

 character*16, allocatable, save :: rf_ele(:)
 logical path_length_patch

! find rf cavities


  j=0
  do i=1,ring%n_ele_track
   if(ring%ele(i)%key == rfcavity$) then
    j=j+1
   endif
  end do

  j_max = j
  allocate (patch(1:j), s_rf(1:j), rf_ele(1:j), ix_ele(1:j))

 if( .not. path_length_patch)then  
  j=0
  do i=1,ring%n_ele_track
   if(ring%ele(i)%key == rfcavity$) then
    j=j+1
    rf_ele(j) = ring%ele(i)%name
    call init_ele(patch(j))
    patch(j)%name = 'PATCH_'//ring%ele(i)%name(1:len_trim(ring%ele(i)%name))
    patch(j)%key = patch$
    patch(j)%s = ring%ele(i-1)%s
    patch(j)%value(z_offset$) = 0.
   endif
  end do
 
  do j = 1,j_max
    call element_locator(rf_ele(j), ring, ix)

    call insert_element(ring,patch(j),ix)
    call lat_make_mat6(ring,-1)
  end do
 endif

  j=0
  do i = 1,ring%n_ele_track
   if(ring%ele(i)%key == patch$ .and. ring%ele(i+1)%key == rfcavity$)then
     j = j+1
     s_rf(j) = ring%ele(i)%s
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
      ring%ele(ix_ele(j))%value(z_offset$) = delta_s * delta_frf/frf
     else
      ring%ele(ix_ele(j))%value(z_offset$) = 0.
    endif
  end do

!  do ii = 1, j_max
!    i = ix_ele(ii)
!    print *,' i, ring%ele(i)%name, ring%ele(i)%value(z_offset$), ring%ele(i)%s '
!    print '(i,1x,a,2e12.4)', i, ring%ele(i)%name, ring%ele(i)%value(z_offset$), ring%ele(i)%s
!  end do

  call lat_make_mat6(ring, -1)
  
  path_length_patch = .true.

  deallocate (patch, s_rf, rf_ele, ix_ele)

  return
  end

!-------------------------------------------------------
