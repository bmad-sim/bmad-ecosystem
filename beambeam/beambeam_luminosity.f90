!........................................................................
!+
! Program    : beambeam_luminosity
!
! Description: Program to compute luminosity
!
! Arguments  :
!
! Mod/Commons:
!
! Calls      :
!
! Author     : 
!
! Modified   :
!-
!........................................................................
!
! $Id$
!
! $Log$
! Revision 1.3  2007/01/30 16:14:30  dcs
! merged with branch_bmad_1.
!
! Revision 1.2.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.2  2005/07/11 14:56:37  sni2
! added damping control to infile, fixed elfile bug
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.h"

program beambeam_luminosity

!  use bmad_struct
!  use bmad_interface
  use bmad
  use bmadz_interface
  use bsim_interface
  use scan_parameters
 
  implicit none

!SIBREN
#if defined(CESR_LINUX) || defined(CESR_WINCVF)
#include "mpif.h"
#endif

  type (lat_struct) ring, ring_ok, ring_in
  type (coord_struct), allocatable :: orb(:)
  type (scan_params_struct) scan_params
  

  integer i, j, ie
  integer i_ele
  integer int_Q_x, int_Q_y, ix, ixx,ixxx
  integer kx
  integer ios
  integer version
  integer n_part, n_turn, particle, i_train, j_car, n_trains_tot, n_cars, slices 

  real(rp), allocatable :: dk1(:) 
  real(rp) phi_x, phi_y
  real(rp) a(6,6), wr(6), wi(6), tune(6)
  real(rp) test
  real(rp) Q_z, current
  real(rp) Q_x, Q_y
  real(rp) min_sig
  real(rp) Q_x_init, Q_y_init
  real(rp) coupling_sb, coupling_wb
  real(rp) sig_in(1:3)

  character lat_file*100, file_name*100, out_file*100, new_lat_file*100
  character date*10, time*10
  character*2 wordx
  character*120 line, last_line

  logical ok, found
  logical lrbbi, beambeam_ip, close_pretz, close_vert, rec_taylor
  logical radiation, parallel
  
  character*120 go

!SIBREN
#if defined(CESR_LINUX) || defined(CESR_WINCVF)
  integer status(MPI_STATUS_SIZE)
#endif
  integer ierr,rank,size
  integer sib_i, lun, damping

!ANDREW
  character*5 fit
  integer readstatus
  integer :: arg_num,iargc
  real(rp) :: final_pos_in(1:4)   ! array in beambeam.in file giving closed orbit coord

  namelist / parameters / parallel, fit, Q_x, Q_y, &
                Q_z, n_part, lat_file, n_turn, particle, &
                i_train, j_car, n_trains_tot, n_cars, current, &
                lrbbi, beambeam_ip, close_pretz,  &
                close_vert, slices, rec_taylor, n_part, min_sig, Q_x_init, Q_y_init, &
                coupling_sb, coupling_wb, radiation, sig_in, go, final_pos_in, damping


global_rank = 0
rank = 0
size=2


!Since the file name is now input form the command line, all the processes can do this at once
!Thus we can save a transmit as well as simplify the conversion to a parallel/non-parallel
!version
!  if(rank.eq.0) then
!     do 
!        print '(a, $)', ' Input command file <CR=beambeam.in>: '
!        read(*, '(a)') file_name
!        call string_trim (file_name, file_name, ix)
!        if (ix .eq. 0) file_name = 'beambeam.in'
!        file_name = 'beambeam.in'

! USE COMMAND LINE ARG TO INPUT FILENAME
        arg_num=iargc()
        if(arg_num==0) then
           file_name='beambeam.in'
        else
           call getarg(1, file_name)
        end if
        call string_trim (file_name, file_name, ix)
        open (unit= 1, file = file_name, status = 'old', iostat = ios)
        
!        if (ios == 0) then
!           !We have a good file, lets send the name to the other processes
!           do sib_i=1,size
!              if(scan_params%parallel) then
!                 call MPI_SEND(file_name,100,MPI_CHARACTER,sib_i,sib_i,MPI_COMM_WORLD,ierr)
!              end if
!           enddo
!           exit
!        else
        if(ios.ne.0)then
           type *
           type *, 'ERROR: CANNOT OPEN FILE: ', trim(file_name)
        endif
        
!     enddo
!  else
!     !we are not rank 0 so we can simply receive the file name
!     if(scan_params%parallel)then
!        call MPI_RECV(file_name,100,MPI_CHARACTER,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
!        call string_trim (file_name,file_name,ix)
!        open (unit=1, file = file_name, status = 'old', iostat = ios)
!     end if
!  endif

! read in the parameters
  parallel = .false.
  particle = positron$
  beambeam_ip = .false.
  close_pretz = .false.
  close_vert = .false.
  rec_taylor = .false.
  slices = 1
  min_sig = 0.
  coupling_sb = 0.01
  coupling_wb = 0.
  damping = 20000
  radiation = .false.
  scan_params%sig_in(1:3) = 0.0
  scan_params%final_pos_in%vec(1:4) = 0.0
  fit = "nope"  ! this default corresponds to no fitting
  go = ".true." ! this default corresponds to no element changes
  
  rewind(1)
  read(1, nml = parameters,iostat=readstatus)
  if(readstatus > 0) then
     if(rank.eq.0) print *,"CAN NOT READ IN SCAN_PARAMS FROM ",file_name
     stop
  end if

  scan_params%Q_z    =        Q_z 
  scan_params%n_part =        n_part 
  scan_params%lat_file =      lat_file 
  scan_params%n_turn =        n_turn 
  scan_params%particle =      particle
  scan_params%i_train =       i_train 
  scan_params%j_car =         j_car 
  scan_params%n_trains_tot =  n_trains_tot 
  scan_params%n_cars =        n_cars 
  scan_params%current =       current
  scan_params%lrbbi =         lrbbi
  scan_params%beambeam_ip =   beambeam_ip 
  scan_params%close_pretz =   close_pretz
  scan_params%close_vert =    close_vert 
  scan_params%slices =        slices 
  scan_params%rec_taylor =    rec_taylor
  scan_params%min_sig    =    min_sig
  scan_params%coupling_sb   = coupling_sb
  scan_params%coupling_wb =   coupling_wb
  scan_params%radiation  =    radiation
  scan_params%sig_in     =    sig_in
  scan_params%fit        =    fit
  scan_params%final_pos_in%vec(1:4) = final_pos_in(1:4)
  scan_params%parallel   =    parallel
  scan_params%damping    =    damping

!SIBREN
#if defined(CESR_LINUX) || defined(CESR_WINCVF)
if(scan_params%parallel)then
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
   global_rank = rank
   call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
   size = size-1
   if(rank.eq.0) print *,"Using parallel code"
else
print *,"Using non-parallel code"
end if
#else

if(scan_params%parallel)then
   print *,"Parallel operation requested on a platform that does not support it."
   print *,"Using non-parallel code."
end if
#endif

! init lattice
  if(lat_file(1:8) == 'digested')then
     call read_digested_bmad_file (lat_file, ring, version)
  else
     call fullfilename(lat_file,lat_file)
     call bmad_parser (lat_file, ring)
  endif

  if(.not.((trim(go)=='.true.').or.(trim(go)=='.false.'))) then
     lun=lunget()
     open(UNIT=lun,FILE=trim(go),STATUS="OLD")
  end if

  do while (.not. (trim(go)=='.true.'))

!SIBREN
     if(rank.eq.0) then
        !We only want to have to enter things once, so only rank 0 does this
        type '(a, $)', ' BEAMBEAM: element change or GO> '
        if(trim(go)=='.false.') then
           read(*, '(a)') line
        else
           read(lun, '(a)') line
        end if
        
        ix = index(line, '!')
        if (ix /= 0) line = line(:ix-1)        ! strip off comments
        
        call str_upcase(line, line)
        call string_trim(line, line, ix)
        
        if (ix == 0) then       ! nothing typed. do the same thing
           line = last_line
        endif
        

     end if
#if defined (CESR_LINUX) || defined(CESR_WINCVF)
     if(scan_params%parallel)then
        call MPI_BCAST(line,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
     end if
#endif
     last_line = line
     
     call str_upcase(line,line)
     if(line(1:1) .eq. 'G')exit
     if(line(1:2) == 'EX' .or. line(1:2) == 'QU')then
#if defined (CESR_LINUX) || defined(CESR_WINCVF)
        if(scan_params%parallel)then
           call MPI_FINALIZE(ierr)
        endif
#endif
        stop
     endif
     
     !  lat_file = line
     
     call find_change( line, ring)
     
  end do

  if(.not.((trim(go)=='.true.').or.(trim(go)=='.false.'))) then
     close(unit=lun)
  end if

  ix = index(file_name,'.')
  out_file = file_name(1:ix)//'hed'

  ixx = max(index(lat_file,'/'),1)
  ixxx = len_trim(lat_file)
  kx = ixxx-ixx+1
  new_lat_file = lat_file(ixx:ixxx)
  new_lat_file(kx+1:kx+1) = "'"

  call date_and_time(date, time)

!SIBREN
!this file only need be written once, let 0 do it
  if(rank.eq.0) then
     open(unit=21, file = out_file)
     write (21, '(a13,a100)') 'Lattice   = `', new_lat_file
     write (21, *) 'File_name = `', trim(file_name)// "'"
     write (21, '(a13,1x,a4,a1,a2,a1,a2,1x,a2,a1,a2,a1,a2,a1)') 'data_date = `', &
          date(1:4),'-',date(5:6),'-',date(7:8),time(1:2),':',time(3:4),':',time(5:6), "'"
     write (21,'(a11,i5)' )'  n_part = ',n_part
     write (21,'(a11, i7)')'  n_turn = ',n_turn
     write (21, *) 'Return'
     close( unit=21)
     
     out_file = file_name(1:ix)//'out'
     open(unit=21, file = out_file)
     write (21, '(a13,a100)') 'Lattice   = `', new_lat_file
     write (21, *) 'File_name = `', trim(file_name)// "'"
     write (21, '(a13,1x,a4,a1,a2,a1,a2,1x,a2,a1,a2,a1,a2,a1)') 'data_date = `', &
          date(1:4),'-',date(5:6),'-',date(7:8),time(1:2),':',time(3:4),':',time(5:6), "'"
     write (21, *) 'Return'
     write (21, *)
     write(21,1) Q_x, Q_y, Q_z
1    format(1x,'!','  Q_x = ',f5.3, '  Q_y = ',f5.3,  &
          ' Q_z = ',f5.3)
2    format(1x,'! Traj',i3,' start ',2x,6e10.2)
     
     
     write(21,*)
     
     write(21,'(a1,2a10,2x,2a10,1x,3a12,1x,3a12,2x,a5,2x,a4,2x,a9)')'!',' Q_x(sb)  ',' Q_y(sb)  ',&
          ' Q_x(bb)  ',' Q_y(bb)  ', &
          '  sigx_in  ','  sigy_in   ','  sigz_in   ', &
          '  sigx_out ','  sigy_out  ','  sigz_out  ', &
          'n_out',' lum', '  current'
  endif

  do i = 1, ring.n_ele_max
    if(ring.ele(i).value(x1_limit$) == 0.)ring.ele(i).value(x1_limit$) = 0.05
    if(ring.ele(i).value(x2_limit$) == 0.)ring.ele(i).value(x2_limit$) = 0.05
    if(ring.ele(i).value(y1_limit$) == 0.)ring.ele(i).value(y1_limit$) = 0.05
    if(ring.ele(i).value(y2_limit$) == 0.)ring.ele(i).value(y2_limit$) = 0.05
  enddo

  ring.param.aperture_limit_on = .true.

  allocate(dk1(ring%n_ele_max))
  allocate(orb(0:ring%n_ele_max))

  if(Q_x_init == 0)Q_x_init = Q_x
  if(Q_y_init == 0)Q_y_init = Q_y

! find which quads to change

  orb(0)%vec = 0.
  call twiss_at_start(ring)
  call closed_orbit_calc(ring, orb, 4)
  call lat_make_mat6(ring, -1, orb)
  call twiss_at_start(ring)
  call twiss_propagate_all(ring)
  call choose_quads(ring, dk1)

  int_Q_x = int(ring%ele(ring%n_ele_track)%a%phi / twopi)
  int_Q_y = int(ring%ele(ring%n_ele_track)%b%phi / twopi)
  phi_x = (int_Q_x + Q_x_init) * twopi
  phi_y = (int_Q_y + Q_y_init) * twopi

  call custom_set_tune(phi_x, phi_y, dk1, ring, orb, ok)    

    phi_x = (int_Q_x + Q_x) * twopi
    phi_y = (int_Q_y + Q_y) * twopi


    scan_params%file_name = file_name(1:ix-1)

    call out_io(s_info$,"BEAMBEAM",'scan_params%file_name  ' // scan_params%file_name)
    ring_in = ring
    if(size.eq.0) then
       call out_io(s_abort$,"BEAMBEAM", "You must have at least 2 processes running to continue",&
            "(one to track particles, at least one to calculate luminocity)",&
            "beambeam aborting after having digested lattice file")
    else
       call beambeam_scan(ring_in, scan_params, phi_x, phi_y)
       if(rank.eq.0) then
          write(21,'(2f10.5,2x,2f10.5,3e12.4,2x,3e12.4,2x,i5,2x,e12.4,2x,e10.2)')Q_x,Q_y, &
               ring_in%a%tune/twopi, ring_in%b%tune/twopi, &
               scan_params%sig_in(1:3), &
               scan_params%sig_out(1:3), scan_params%n_part_out, &
               scan_params%lum, scan_params%current
          
          print *,current,Q_x, Q_y, scan_params%lum
          
          
          close(unit = 21)
          
       endif
    endif
    deallocate(orb)
#if defined (CESR_LINUX) || defined(CESR_WINCVF)
    if(scan_params%parallel)then
       call MPI_FINALIZE(ierr)
    end if
#endif
    
  end program beambeam_luminosity
