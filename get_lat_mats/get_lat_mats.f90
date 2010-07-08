program get_lat_mats
!+
! [cesr.bmad.cesr]get_lat_mats
! extract 1 turn matrix and 
! tranfer matrices from one set of detectors to the next.
!
! mjf 2006.07.07
!-

  use bmad
  use cesr_basic_mod

  implicit none

  include 'det_list.inc'
  
  character*40 :: lattice, cur_lattice
  character*80 :: lat_file
  character*10 :: start
  character*14 ::  filename(2), outfilename
  real(8) :: xtune, ytune, dx, dy, startx, starty, beta(2), gamma
  real(8), allocatable :: dk1(:)
  real(8) :: Ga(2,2), Gb_inv(2,2), C(2,2), Cbar(2,2), temp(2,2)
  type (lat_struct) :: ring, rev_ring
  type (coord_struct), allocatable :: orb(:)
  logical :: ok, zero(2)
  integer, allocatable :: det_ix(:)
  integer :: i, ip_l3_ix, nfiles, file
  integer :: lun, ios
  real(rp) :: mat6(6,6), det
!  integer, parameter :: bpm_tot = 14
!  character(8), parameter :: use_bpm(bpm_tot) = (/ &
!       "DET_00W ", "DET_01W ", "DET_02W ", "DET_04W ", "DET_07W ", &
!       "DET_08W ", "DET_08AW", "DET_09W ", "DET_10W ", "DET_11W ", &
!       "DET_12W ", "DET_02E ", "DET_01E ", "DET_00E "/)
  integer :: bpm_loc(bpm_tot)
  logical :: err
  ok = .false.
  zero(1) = .false.
  zero(2) = .false.
  call mpm_init('BMAD')
  start = "IP_L0"
 ! call getlat( cur_lattice )

 ! call choose_cesr_lattice( lattice, lat_file, cur_lattice, ring )

  call bmad_parser( '/nfs/cesr/mnt/lattice/cesr/bmad/bmad_chess_20090225.lat', ring )

  allocate( orb(0:ring%n_ele_max) )
  allocate(dk1(ring%n_ele_max))
  orb(0)%vec(6) = 0
  call set_on_off(elseparator$, ring, off$)
!  call bmad_parser( "bmad.", ring )

  call twiss_at_start( ring )

  call twiss_propagate_all( ring )

  do i=1, bpm_tot
     call element_locator(use_bpm(i), ring, bpm_loc(i))
  enddo

!  do i=1, bpm_tot
!     if (bpm_loc(i)>0) then
!     Print *, "Beta:", ring%ele(bpm_loc(i))%a%beta, &
!          ring%ele(bpm_loc(i))%b%beta
!     Print *, "Alpha:", ring%ele(bpm_loc(i))%a%alpha, &
!          ring%ele(bpm_loc(i))%b%alpha
!     else
!        Print *, "Could not find", i
!     endif
!  enddo
!  Print *, "Starting point (ex DET_00W): "
!  Read (*, 12) start
!  start = trim(adjustl(start))
!12 format (a10)
!Print *, "Start: ", start

! Return indices for all ele's of type marker into
! the det_ix array
  call elements_locator( 'marker:*', ring, det_ix, err)
  if (err) then
     print *, "Decoding error in elements_locator using marker$"
  endif

! Return index of element in ring with name start
  call element_locator( start, ring, ip_l3_ix)
  
  do i=1, ring%n_ele_max
     if (ring%ele(i)%key==quadrupole$) then
        if (ring%ele(i)%value(k1$) < 0) then
           dk1(i) = -1
        else if (ring%ele(i)%value(k1$) > 0) then
           dk1(i) = 1
        else
           dk1(i) = 0
        endif
     else
        dk1(i) = 0
     endif
  enddo

  Print *, "How many files to generate?"
  read *, nfiles

  Print *, "Starting tune shift in X?"
  read *, startx
  Print *, "Starting tune shift in Y?"
  read *, starty

  Print *, "Change in X tune:"
  Read *, dx
  Print *, "Change in Y tune:"
  Read *, dy

  startx = startx + ring%ele(ring%n_ele_track)%a%phi
  starty = starty + ring%ele(ring%n_ele_track)%b%phi

35 format ('lat_mats0', i1, '.out')
65 format ('lat_mats', i2, '.out')

  do file=1, nfiles
     if (file < 10) then 
        write (outfilename, 35) file
     else
        write (outfilename, 65) file
     endif
     
     lun = lunget()
     open( unit = lun, file = outfilename, iostat = ios )
     if (ios/=0) then
        write(*) 'Error opening file: ', outfilename
        stop
     endif

     !  call twiss_at_start(ring)
     !  Print *, "twiss", ring%ele(0)%x
     Print *, "Lattice: ", ring%name
     Print *, "X tune in radians:", ring%ele(ring%n_ele_track)%a%phi
     Print *, "Y tune in radians:",  ring%ele(ring%n_ele_track)%b%phi

     xtune = startx + dx*(file-1)
     ytune = starty + dy*(file-1)
  
     call set_tune(xtune, ytune, dk1, ring, orb, ok)
     if (.not. ok) then
        Print *, "Tune was not successful."
     endif

     Print *, "New X tune:", ring%ele(ring%n_ele_track)%a%phi
     Print *, "New Y tune:", ring%ele(ring%n_ele_track)%b%phi

87   format ("lat", i2, ".beta")
86   format ("lat0", i1, ".beta")
89   format ("lat", i2, ".cb")
88   format ("lat0", i1, ".cb")
     if (file < 10) then
        write (filename(1), 86) file
        write (filename(2), 88) file
     else
        write (filename(1), 87) file
        write (filename(2), 89) file
     endif


     open (unit = 27, file=filename(1), action="write", position="rewind")
     open (unit = 57, file=filename(2), action="write", position="rewind")
34   format (5x, f11.7, ' ,', 2x, f11.7)
91   format (5x, f11.7, ' ,', 2x, f11.7, 2x, ' ,', f11.7)
     !  write (27, "(15x, a12)") "Beta"
     !Some blank lines (for ease of concatenating files)
     write (27, '(1x)')
     write (27, '(1x)')
     write (57, '(1x)')
     write (57, '(1x)')

     write (27, 34) ring%ele(ring%n_ele_track)%a%phi, &
          ring%ele(ring%n_ele_track)%b%phi
     write (57, 34) ring%ele(ring%n_ele_track)%a%phi, &
          ring%ele(ring%n_ele_track)%b%phi

     do i=1, bpm_tot
        C = ring%ele(bpm_loc(i))%c_mat
        
        beta(1) = ring%ele(bpm_loc(i))%a%beta
        beta(2) = ring%ele(bpm_loc(i))%b%beta

        Ga(1,1) = 1 / sqrt(beta(1))
        Ga(1,2) = 0
        Ga(2,1) = ring%ele(bpm_loc(i))%a%alpha / sqrt(beta(1))
        Ga(2,2) = sqrt(beta(1))
        
        Gb_inv(1,1) = sqrt(beta(2))
        Gb_inv(1,2) = 0
        Gb_inv(2,1) = -ring%ele(bpm_loc(i))%b%alpha / sqrt(beta(2))
        Gb_inv(2,2) = 1 / sqrt(beta(2))

        temp = matmul(Ga, C)
        Cbar = matmul(temp, Gb_inv)
        !gamma =  sqrt( 1 - det(C))
        gamma = sqrt(1 - (Cbar(1,1)*Cbar(2,2) - Cbar(1,2)*Cbar(2,1)))
        Cbar(1,2) = Cbar(1,2) / gamma
        Cbar(2,2) = Cbar(2,2) / gamma
        beta(1) = (gamma**2) * beta(1)
        beta(2) = (gamma**2) * beta(2)

        write (27, 34) beta(1), beta(2)
        write (57, 91) Cbar(1,1), Cbar(1,2), Cbar(2,2)
     enddo
     
     !  write (27, "(10x, a12)") "Alpha"
     !  do i=1, bpm_tot
     !write (27, 34) ring%ele(bpm_loc(i))%a%alpha, ring%ele(bpm_loc(i))%b%alpha
     !  enddo
     
     ! Calculate and write out the single turn 6x6 xfr matrix
     call transfer_matrix_calc( ring, .false., mat6 )
     call mat_type( mat6, -lun, "One turn transfer matrix" )
     write (lun,*) " "
     
     print *, "size of det_ix:", size(det_ix)

     do i=2,size(det_ix)
        !   if (det_ix(i) >= ip_l3_ix) exit
        call transfer_matrix_calc( ring, .false., mat6, &
             ix1 = det_ix(i-1), ix2 = det_ix(i) )
        !    print *, ring%ele(det_ix(i-1))%name, " to ", &
        !         ring%ele(det_ix(i))%name, "transfer matrix: "
        write (lun,*) ring%ele(det_ix(i-1))%name, " to ", &
             ring%ele(det_ix(i))%name, "transfer matrix: "
        call mat_type( mat6, -lun)
        write (lun,*) " "
        !    print *
     enddo

! No longer want to reverse the ring in the east,
! so remove the following:
!   call lat_reverse( ring, rev_ring )
!   call elements_locator( marker$, rev_ring, det_ix )
!   call element_locator( start, rev_ring, ip_l3_ix )
! 
!   do i=2,size(det_ix)
!     if (det_ix(i) >= ip_l3_ix) exit
!     call transfer_matrix_calc( rev_ring, .false., mat6, &
!          ix1 = det_ix(i-1), ix2 = det_ix(i) )
!     print *, rev_ring%ele(det_ix(i-1))%name, " to ", &
!          rev_ring%ele(det_ix(i))%name, "transfer matrix: "
!     write (lun,*) rev_ring%ele(det_ix(i-1))%name, " to ", &
!          rev_ring%ele(det_ix(i))%name, "transfer matrix: "
!     call mat_type( mat6, lun )
!     write (lun,*) " "
!     print *
!   enddo
! 
     close(lun)

     print *, "Wrote file: ", outfilename

  enddo

  call goodbye

end program get_lat_mats

  
