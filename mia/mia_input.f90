module mia_input

  use mia_types
 
  type (data_file), allocatable :: file(:) !File being used
contains

  subroutine read_bpm_data(data,iset,nset)

    !Subroutine opens a file of bpm data and reads it into the module.
    !Routine reads in how many processors are active and assigns that
    !number to the variable bpmproc

    !tempx and tempy are temporary holding places for the first x and y    
    !Values from each processor; these values get assigned to the first 
    !location in the x and and y array respectively for its processor

    implicit none
    type(data_set)data                     !Data set
    integer :: OpenStatus, InputStatus, &  !Input statuses
         proc,&                            !Counts through processors
         turns,&                           !Number of turns
         C, I, S, &                        !Counters
         places, power
    integer :: iset, nset                  !Number of current set
    real(rp) :: tempx, tempy, &                !Temporary variables in x and y
         sum, ave                          !Sum and average
    logical :: continue, &                 !Whether to ask for another filename
         trunc
    real(rp), allocatable :: newt(:), &    !Temporary array
         pos_temp(:,:)

    character(30) :: Plane(2)          !See next line:
    data Plane/' data file 1.',' data file 2.'/

    continue = .true.
  

    DO WHILE (continue == .true.)
177    format (a120)
       if (fileread(iset)==.false.) then
          Write (*, '(1x, 3a)', advance = "no") "Enter path or number",&
               " beginning with # (ex. #03081) of", Plane(iset)
          read (*,177), data%filename
          data%filename = trim(data%filename)
!          call strget('Enter filename: ', data%shortName)
       endif

!      If data files are in the same location with the same naming pattern,
!      this section can be modified and uncommented to make entering
!      filenames easier. 
!       if (scan(data%shortName, "#")>0) then
          !Remove first character from input (#) and concatenate with
          !the path for data files.
!          data%shortName = data%shortName(2:len_trim(data%shortName))
!          data%shortName = "cbpm_" // trim(data%shortName) // ".raw"
!          data%filename = "log3$disk:[cesr.cesrbpm.raw.03]" // &
!               data%shortName
!       else
!          data%filename = trim(data%shortName)
!       endif

       if (scan(data%filename, "#")>0) then

          data%shortName = data%filename(scan(data%filename, &
               "#")+1:len_trim(data%filename))
          data%filename = "./data/" // data%shortName
       else
          if (scan(data%filename, "/")>0) then
             data%shortname = data%filename(scan(data%filename, &
                  "/", .true.)+1:len_trim(data%filename))
          endif
       endif

       open (unit = 1, file = data%filename, status = "old", &
            iostat = openstatus)

       !Check for bad filenames
       if (openstatus > 0) then
          print *, "*** Cannot open file ", data%filename, " ***"
          print *, "*** Try again. ***"
          continue = .true.
          fileread(iset) = .false.
       else
          continue = .false.   
       endif
    enddo

!Truncation was for testing purpose. If needed, uncomment this section and 152.
!151    call logic_get( 'T', 'C', ' (T)runcate values or (C)ontinue', trunc)
!    if (trunc) then
!       Print *, "Input number of decimal places (no greater than 7)"
!       read *, power
!       places = 10**power
!    endif

10  format (/ 1x, t28, i3)   !reads in number of processors

40  format (10x, a14,  t29, f10.7, t41, f10.7, t125, i6) 
    !reads in number of 
    !turns and the first 
    !values for x and y

30  format (1x, t29, f10.7, t41, f10.7) !reads in rest of x and y

    Read (1, 10, iostat = inputstatus) data%bpmproc
    if (inputstatus > 0) stop "*** INPUT ERROR ***"
    if (inputstatus < 0) stop "*** Not enough data ***"
    allocate (data%cdata(data%bpmproc))      !allocates the cdata array to be
    !as large as the number of active
    !bpm processors
    if (.not. ALLOCATED(file)) then
       allocate(file(nset))
    endif
    allocate(file(iset)%proc(data%bpmproc))

    Processor: do proc = 1, data%bpmproc

       Read (1, 40, iostat = inputstatus) file(iset)%proc(proc)%label, &
            tempx, tempy, data%numturns
       if (inputstatus > 0) stop "*** INPUT ERROR ***"
       if (inputstatus < 0) stop "*** Not enough data ***"	

       allocate (data%cdata(proc)%x(0:data%numturns-1))
       !allocates the x and y arrays to be as large as the number of turns 
       allocate (data%cdata(proc)%y(0:data%numturns-1)) 
       data%cdata(proc)%x(0) = tempx !assigns the temp values to the 
                                     !first posisitions
       data%cdata(proc)%y(0) = tempy !of the x and y arrays

       turn: do turns = 1, data%numturns-1
          !starts at 1 and goes
          !to turns-1 because
          !temps are read and 
          !stored earlier
          read (1, 30, iostat = inputstatus) &
               data%cdata(proc)%x(turns), data%cdata(proc)%y(turns)
          if (inputstatus > 0) stop "*** INPUT ERROR ***"
          if (inputstatus < 0) stop "*** Not enough data ***"	
!Truncation was for testing. If needed, uncomment this if statement and 151.
!          if (trunc) then
!             data%cdata(proc)%x(turns) = &
!                  floor( data%cdata(proc)%x(turns) * places ) / (places*1.0)
!             data%cdata(proc)%y(turns) = &
!                  floor(data%cdata(proc)%y(turns)*places) / (places*1.0)
!          endif
       End do turn
    end do Processor

    close (1)

    allocate(pos_temp(data%numturns, 2*data%bpmproc)) 
    !allocates room for the 
    !posistion history matrix to
    !make it have rows = number
    !of turns and columns = 
    !twice the number of 
    !processors that are active
    do proc = 1, data%bpmproc

       pos_temp(:,2*proc-1) = data%cdata(proc)%x(:)  &
            /sqrt(float(data%numturns))  !feeds x arrays into odd columns

       pos_temp(:,2*proc) = data%cdata(proc)%y(:)    &
            /sqrt(float(data%numturns))  !feeds y arrays into even columns
    end do


    sum = 0.0
    allocate(data%poshis(data%numturns, 2*data%bpmproc)) 
    !readies room for new 
    !posistion history matrix
    allocate(newt(data%numturns))  
    !creates a temp. place for transitional data
    do C = 1, 2*data%bpmproc

       do I = 1, data%numturns
          sum = sum + pos_temp(I,C)    !sums the columns of x or y
       end do

       ave = sum/data%numturns     !takes average of the sum of a column

       do s = 1, data%numturns
          newt(s) = pos_temp(s,C) - ave  
          !subtracts average value from each 
          !value in column and assigns the new
          !value to a place
          !in the temp. array "newt"
       end do

       data%poshis(:,C) = newt    
       !assigns the new data in the temp. array to the
       !second position history matrix         
    end do

    deallocate(newt)        !Deallocate temporary arrays.
    deallocate(pos_temp)

  end subroutine read_bpm_data

  subroutine locate_bpm (iset, data)

    !
    !Creates a list of the bpms that are active.  Keeps track of bpms in West 
    !or East 
    !Assigns BPMs a number in addition to their label. This is positive
    !in the west and negative in the east
    !

    implicit none
    type(data_set)data       !Data set
    integer :: nchar, &      !Place were a certain string was found.
         n_bpm,&             !Number of BPMs
         iset,i, &              !Counters
         length              !Length of a string 
    character(5) :: temp     !Holds temporary file number

    do n_bpm = 1, NUM_BPMS

       nchar = scan(file(iset)%proc(n_bpm)%label, 'BPM') + &	
            scan(file(iset)%proc(n_bpm)%label, 'bpm')	
       if (nchar > 0 .and. nchar < 7 ) then
          length = len_trim(file(iset)%proc(n_bpm)%label)
          temp = trim(file(iset)%proc(n_bpm)%label(nchar+4:length))
       end if
       nchar = scan(temp, 'W') + &   !Scan temp so W can be truncated
            scan(temp, 'w')          !Scans label for W(est)
       if (nchar > 0 ) then 
          file(iset)%proc(n_bpm)%is_west = .true.
          temp = temp(1:nchar-1) !Remove W from temp
       else                               !Look for E so it can be removed
          nchar = scan(temp, 'E') + scan(temp, 'e')
          temp = "-" //  temp(1:nchar-1)    !Add negative sign
       endif

       nchar = scan(temp, 'A') + &
            scan(temp, 'a')  !Scan for A (ex BPM 8AW)
       if (nchar>0) then
          temp = temp(1:nchar-1) // ".5"  !Remove "A" from temp and add 0.5
       endif
       read (temp, *), file(iset)%proc(n_bpm)%number

       if (file(iset)%proc(n_bpm)%number == 0 .and. &
            .not. file(iset)%proc(n_bpm)%is_west) then
          file(iset)%proc(n_bpm)%number = -0.000001
       end if

    end do

  end subroutine locate_bpm


  subroutine find_L(nset,data)

    !
    !Uses list of bpms from prior subroutine to find l's that are included in 
    !the data file being studied
    !

    implicit none
    integer :: inputstatus, openstatus,&  !Status indicators
         num, i, j, i_file, nset, k, li, &!Counters and nset--number of sets
         nfile, rnum                    !Counter and actual number of BPM pairs
    logical :: found_one(2)               !If a BPM pair has been found
    type(data_set) data(*)                !All data
    type(known_spacings), allocatable ::  bpm_old(:)   !BPM data read in

    !Placed in different location because
    !a BPM may not be in use.

    !knownl.inp contains a list of BPMs with known spacing.
    !Change to accept different locations for knownl.inp...
    open (unit = 2, file = "./data/knownl.inp", &
         status = "old", iostat = openstatus)
    if (openstatus > 0) Stop "*** Cannot open file knownl.inp ***"

60  format (1x, i5,/)
61  format (1x, a13, 3x, a13, 2x, f7.3)

    read (2,60, iostat = inputstatus) num
    if (inputstatus > 0) stop "*** INPUT ERROR ***"
    if (inputstatus < 0) stop "*** Not enough data ***"	
    allocate (bpm_old(num))
    call allocate_bpm_pairs(num)

    bpm_old(1)%number = num
    !Change bpm_pairs%number not to be in every bpm of the array?
    bpm_old(:)%in_use = .false.
    rnum = 0                          !Actual number of BPMs in use

    do i = 1, num
       read (2,61, iostat = inputstatus) (bpm_old(i)%bpm_name(j),j=1,2), &
            bpm_old(i)%length
       if (inputstatus > 0) stop "*** INPUT ERROR ***"
       if (inputstatus < 0) stop "*** Not enough data ***"	
       bpm_old(i)%in_use = .true.
    end do

    !This could be made more efficient by finding a BPM number for 
    !the pairs and scanning for that.
    !Also, this checks all BPMs vs all pairs--change to stop when a 
    !match is found. 
    bpm_old(1)%has_one = .true.
    do i_file = 1, nset                       !File
       do i = 1, num                          !BPM pairs
          bpm_old(i)%bpm_pntr(1)=0
          bpm_old(i)%bpm_pntr(2)=0
          do j = 1, NUM_BPMS           !All BPMs
             do k = 1, 2                      !First and second of the BPM pair
                if (file(i_file)%proc(j)%label == bpm_old(i)%bpm_name(k)) then
                   bpm_old(i)%bpm_pntr(k) = j
                end if
             end do
          end do
          if (bpm_old(i)%bpm_pntr(1) == 0 .or. &
               bpm_old(i)%bpm_pntr(2) == 0) then
             bpm_old(i)%in_use = .false.
          else 
             bpm_old(i)%in_use = .true.
             found_one(i_file) = .true.
          end if
       end do
       !Unnecessary?
       bpm_old(1)%has_one = bpm_old(1)%has_one .and. found_one(i_file)
    end do

    !Copies BPMs pairs that were found into bpm_pairs.
    !Not all BPMs in knownl.inp may be in use.
    !    do nfile = 1, nset                       !File
    do li = 1, num                        !BPM pairs from knownl.inp
       if (bpm_old(li)%in_use == .true.) then
          rnum = rnum+1
          !Files 1 and 2 should be the same
          !Remove file from bpm_pairs?
          !Can check consistancy in previous do loops.
          bpm_pairs(rnum) = bpm_old(li)
          bpm_pairs(rnum) = bpm_old(li)
       endif
    enddo
    if (rnum > 0) then
       bpm_pairs(1)%has_one = .true.
    else
       bpm_pairs(1)%has_one = .false.
    endif

    !    enddo
    bpm_pairs(1)%number = rnum
    deallocate (bpm_old)    

  end subroutine find_L

  subroutine get_ele_sPos(iset)
    !
    !Reads in a file containing bpms with their s positions and
    !element numbers and matches them with bpms in the input file.
    !
    character(7) :: detName
    integer :: eleNum, openstatus, i, inStat, iset, bpmInt
    real(rp) :: sPos, bpmNum
    logical :: continue, theEnd
    open (unit = 273, file = "./data/one_ring.det", &
         status = "old", iostat = openstatus)
    if (openstatus > 0) Stop "*** Cannot open file one_ring.det ***"

188 format (1x,a7,5x,f10.6,5x,i3)
56  format (2x,i2)

    theEnd=.false. !It's just the beginning
    continue = .true.

    do while (continue)

       read (273,188,iostat=inStat) detName, sPos, eleNum
       if (inStat>0) stop "*** Error reading one_ring.det ***"
       if (inStat < 0) then
          continue = .false.
          cycle
       endif
       detName = trim(adjustL(detName))
       
       select case(detName(1:2))
       case ('EN')
          !Location of the end of the ring: the length of the machine
          theEnd=.true.
          endLoc=sPos
       case ('DT')
          read (detName,56) bpmInt
          bpmNum = bpmInt
          if (detName(5:5) == "A" .or. detName(6:6) == "2") then
             bpmNum = bpmNum+0.5
          endif
          if ( detName(4:4) == "E" .or. detName(5:5) == "E") then
             bpmNum  = -bpmNum
             if (bpmNum == 0) then
                bpmNum = -0.000001 !Since there is no -0...something close.
             endif
          endif
       case default
          Print *, "Error in reading BPM S Positions"
          exit
       end select

       if (.not. theEnd) then
          do i=1,num_bpms
             if (file(iset)%proc(i)%number == bpmNum) then
                file(iset)%proc(i)%sPos = sPos
                file(iset)%proc(i)%eleNum = eleNum
             end if
          end do
       end if
    end do

  end subroutine get_ele_sPos


  subroutine match_processors(data)

    !
    !This routine makes sure the data being compared from
    !each file has matching processors (if there are multiple files)
    !

    integer :: nbpm,&            !BPM number
         n_file_1, n_file_2, &   !File counters?
         i                       !counter
    type(data_set) :: data(*)    !Data set

    n_file_1 = 1
    n_file_2 = 1
    nbpm = 1

    do while (n_file_1 <= data(1)%bpmproc)
       if (file(1)%proc(n_file_1)%label == &
            file(2)%proc(n_file_2)%label) then
          data_struc%proc(nbpm) = file(1)%proc(n_file_1)
          n_file_1 = n_file_1 + 1 	
          n_file_2 = n_file_2 + 1
          nbpm = nbpm + 1 
       else 
          do i = n_file_2 +1, data(2)%bpmproc 
             if (file(1)%proc(n_file_1)%label == &
                  file(2)%proc(i)%label) then
                data_struc%proc(nbpm) = file(1)%proc(n_file_1)
                n_file_1 = n_file_1 + 1 	
                n_file_2 = i + 1
                nbpm = nbpm + 1
                go to 100
             end if
          end do
          n_file_1 = n_file_1 + 1
       end if
100    continue 
    end do


    if (data(1)%bpmproc /= nbpm - 1) then
       print *, "Number of Processors don't match; Mode not supported"
       stop
    end if

  end subroutine match_processors

  subroutine deall_file
    deallocate(file)
  end subroutine deall_file

end module mia_input
