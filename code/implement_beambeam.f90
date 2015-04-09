!+           
! Subroutine Implement_Beambeam(ring,sigx,sigy,sigz,density) 
!
! Insert a beambeam element at each bend as a electron cloud 
!
! Input:
!   ring -- lat_struct: New elements will be inserted
!   sig_x -- real: cloud width
!   sig_y -- real: cloud height
!   sig_z -- real: cloud length
!   density -- real: density of the cloud
!
! Output:
!   ring -- lat_struct: ring with new elements inserted
!

  subroutine implement_beambeam(ring,sigx,sigy,sigz,density) 

  use bookkeeper_mod

 implicit none

 type (lat_struct) ring
 type (ele_struct), allocatable :: beambeam_ele(:)

 integer i,j
 integer j_max
 integer ix, ii
 integer nslice/5/

 real(rp) sigx,sigy,sigz
  real(rp), save :: sigx_s,sigy_s,sigz_s
 real(rp) delta_s
 real(rp) density
 real(rp), allocatable, save :: length(:)
 
 character*2 number

 logical beambeam/.false./

! find bends
 sigx_s =sigx
 sigy_s = sigy
 sigz_s = sigz

 ring%param%n_part = sigx_s*sigy_s*sigz_s * density
 
 if( .not. beambeam)then
  j=0
  do i=1,ring%n_ele_track
   if(ring%ele(i)%key == sbend$ .or. ring%ele(i)%key == rbend$) then
    j=j+1
   endif
  end do

  j_max = j
  allocate (beambeam_ele(1:j), length(1:j))

  j=0
  number = '00'
  do i=1,ring%n_ele_track
   if(ring%ele(i)%key == sbend$ .or. ring%ele(i)%key == rbend$) then
    j=j+1
    length(j) = ring%ele(i)%value(l$)
    if(j<10)write(number(2:2),'(i1)')j
    if(j>=10)write(number(1:2),'(i2)')j
    if(j>99) print *,' Too many bends for implement beambeam'
    call init_ele(beambeam_ele(j))
    beambeam_ele(j)%name = 'BEAMBEAM_'//number//ring%ele(i)%name(1:len_trim(ring%ele(i)%name))
!    print '(a,i10,1x,a)',' j, beambeam_ele(j)%name ', j, beambeam_ele(j)%name
    beambeam_ele(j)%key = beambeam$
    beambeam_ele(j)%s = ring%ele(i-1)%s
    beambeam_ele(j)%value(sig_x$) = sigx
    beambeam_ele(j)%value(sig_y$) = sigy
    beambeam_ele(j)%value(sig_z$) = sigz
    beambeam_ele(j)%value(n_slice$) = nslice
    beambeam_ele(j)%value(charge$) = -length(j)/sigz

!    beambeam_ele(j)%tracking_method = custom$
    beambeam_ele(j)%mat6_calc_method = tracking$
   endif
  end do
 
  do j = 1,j_max
    ix = element_at_s(ring,beambeam_ele(j)%s, .true.)
!   print *,j,rf_ele(j),ix, ring%ele(ix)%name
    call insert_element(ring, beambeam_ele(j),ix)
    if (.not. bmad_com%auto_bookkeeper) call lattice_bookkeeper(ring)
    call lat_make_mat6(ring,-1)
!    print *,' ring%n_ele_max, ix = ', ring%n_ele_max,ix
!    print '(i10,1x,2a16)',ix,ring%ele(ix)%name, ring%ele(ix+1)%name

  end do

 beambeam =.true.
!deallocate (beambeam_ele)
endif ! if not beambeam



  return
  end subroutine

!-------------------------------------------------------
