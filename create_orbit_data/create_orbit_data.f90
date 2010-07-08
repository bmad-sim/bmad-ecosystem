program make_orbit
  use precision_def
  use nr

  implicit none

  include 'det_list.inc'

  integer i,i1,j,k,ibpm, af
  integer pntr
  integer bpm_num, bpm_ptr

  character*2 ew
  character*4 bpm_char
  character*12 label
  character(8) :: labels(140)
  character*3 new_label
  character*5 file
  real(8) gerr, noise
  real(8) :: t(4,4) 
  real(8) :: x_det(4)
  real(8), allocatable :: x(:,:)
  real(8) :: xamp, yamp
  integer :: numturns
  character*10 :: formt
  character*20 turnschar
  type detector_transport
     real(8) :: det(4,4)
     character*12 string
  end type detector_transport

  type (detector_transport) :: ring(140)

  logical :: found
!  integer, parameter :: bpm_tot = 14
!  character(4),parameter :: use_bpm(bpm_tot) = (/ &
!       "00W ", "01W ", "02W ", &
!       "04W ", "07W ", "08W ", &
!       "08AW", "09W ", "10W ", &
!       "11W ", "12W ", "02E ", &
!       "01E ", "00E "/)
  integer :: bpm_to_ring(bpm_tot), lat
  character(14) :: outfilename

  Print *, "Number of turns (multiple of 2): "
  read *, numturns
  allocate (x(numturns, 4))
  write (turnschar,*) numturns
  turnschar = adjustl(turnschar)

35 format ('lat_mats0', i1, '.out')
65 format ('lat_mats', i2, '.out')
  call init_random_seed()

  Print *, "Lat_mats file number: "
  Read *, lat
  if (lat<10) then
     write (outfilename, 35) lat
  else
     write (outfilename, 65) lat
  endif
  open(unit=20, type="old", name=outfilename)

  read(20,1)
1 format(1x)

  do i=1,4
     read(20,2) (t(i,j),j=1,4)
2    format(3x,4f11.6)
  enddo
  read(20,1)
  read(20,1)

  do k = 1 , 140

     read(20,1)
     read(20,1)
     read(20,3) label
!3    format(5x,a6)
3    format(1x,a10)
     label = trim(adjustl(label))
     labels(k) = label
 !    print *, "label: ", labels(k)

!     if (label(1:1) == "_" .or. label(1:2) == "0_")then
!        bpm_num = 0
!        ew = " "
!     else
!        read(label,31) bpm_num,ew
!31      format(i2,a2)
!     endif
     if (label(1:4) == "DET_") then
        if (label(7:7) == "E" .or. label(7:7) == "W") then
           read(label,31) bpm_num,ew
31         format(4x,i2,a2)
!           Print *, bpm_num, ew
        else
           read(label,377) bpm_num!,ew
377        format(4x,i3)
           ew = 'W'
        end if
     else
        bpm_num = 0
        ew = " "
     endif

     if (bpm_num>99) then
        formt = "(i3)"
     else
        formt = "(i2)"
     end if
     write(bpm_char,formt) bpm_num

     new_label=trim(adjustl(bpm_char))//trim(ew)

!     write(ring(k)%string,*) new_label
     write(ring(k)%string,4) new_label
4    format("CESR BPM ",a)


     if(label(1:3)=="00E") then
        pntr = k
     endif

     do i=1,4
        read(20,2) (ring(k)%det(i,j),j=1,4)
     enddo

     read(20,1)
     read(20,1)

  enddo

  close(unit=20)

  !
  !   Calc transport from IP to each detector
  !

!  do k=2,pntr-1 ! detector number
  do k=2,140 ! detector number

     ring(k)%det = matmul(ring(k)%det,ring(k-1)%det)

  enddo

 !  do k=pntr+1,100 ! detector number
 ! 
 !     ring(k)%det = matmul(ring(k)%det,ring(k-1)%det)
 ! 
 !  enddo


  ! do i=1,4
  !  type 2, (t(i,j),j=1,4)
  ! enddo

  ! type 5,ring(1)%string
  !5 format(1x,a12)

  ! do i=1,4
  !  type 2, (ring(1)%det(i,j),j=1,4)
  ! enddo

  Print*, "Enter level of noise to add to data in mm"
  read *, noise

 ! Print *, "X shaking amplitude (mm) (0.0005): "
 ! read *, xamp
 ! Print *, "Y shaking amplitude (mm) (0.00001): "
 ! read *, yamp

  do i1=1,2	! horz then vert

     write(file,6) i1
!6    format("det",i1,".")
6    format(i1,".")
     open(unit=21,name=file,type="new")
     write(21,8), bpm_tot
8    format(1x/28x,i3)
     x=0

     !Amplitude
     if(i1==1) then
        x(1,1)= 0.0005 !xamp
     else
        x(1,3)= 0.000075 !yamp
     endif

     call track(t,x, numturns)

     !+
     do ibpm = 1,bpm_tot
        found = .false.
        do k=1,140 ! detector number
           if (labels(k) == use_bpm(ibpm)) then
              bpm_to_ring(ibpm) = k
              found = .true.
              exit
           end if
        end do
        if (.not. found) then
           print *, "ERROR, couldn't find ", use_bpm(ibpm)
        end if
     end do

     call init_random_seed()
     !   if( k<=15 .or. (k>=pntr .and. k<=pntr+2) ) then
Print *, "Number of turns: ", numturns
     do ibpm = 1, bpm_tot
        bpm_ptr = bpm_to_ring(ibpm)
        write (21, *) "x"
!        Print *, ring(bpm_ptr)%string
        write (21,12) "Location:", &
             ring(bpm_ptr)%string(10:len(ring(bpm_ptr)%string))
        do af=1, 32
           write (21, *) "x"
        end do
12      format (2x,a9, 11x, a)
        do j = 1, numturns  ! turn number

           call gasdev_s(gerr)

           x_det = matmul(ring(bpm_ptr)%det,x(j,:))

           do k = 1, 4
              x_det(k) = x_det(k) + (gerr * noise).   ! mm
           enddo
           
           write(21,7) (x_det(i),i=1,3,2) ! x, y position
7          format(35x,f9.6,2x,f9.6)

        enddo
        !   endif

     enddo

     close(unit=21) 

  enddo

  stop
end program make_orbit

SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
          
  CALL SYSTEM_CLOCK(COUNT=clock)
          
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
          
  DEALLOCATE(seed)
END SUBROUTINE init_random_seed


subroutine track(t,x, numturns)
  implicit none
  integer i,j,k, numturns
  character*4 label
  real(8) :: t(4,4) 
  real(8) :: x(numturns,4)

  do i=2, numturns
     x(i,:)=matmul(t,x(i-1,:))
  enddo

end subroutine track

