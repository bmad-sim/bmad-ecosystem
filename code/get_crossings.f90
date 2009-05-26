
subroutine get_crossings(lat,string,ncross, cross)

 use bmad
 use bunchcross_mod
 use constraints_mod

  implicit none

 type(lat_struct) lat
 type (pc_struct) pc
 type(constraint_struct) con


  integer trains, bunches, spacing
  integer ncross
  integer test
  integer ix
  integer i,l
  character*120 string
  character*64 file
  real(rp)  cross(1000)
  
   file = 'NULL'
   call string_trim(string, string, ix)
     print '(a)',string
     if(string(1:4) /= 'PRET')then 
      print *,' SOMETHING WRONG IN GET_CROSSINGS '
      stop
     else
      call string_trim(string(ix+1:),string,ix)
!      print *,string(1:ix)
      read(string(1:ix),'(i)',err=99)test
      print *,test
      if(test > 0 .and. test < 20)then 
         trains = test
         call string_trim(string(ix+1:),string,ix)
         read(string(1:ix),'(i)')bunches
         call string_trim(string(ix+1:),string,ix)
         read(string(1:ix),'(i)')spacing
 print *,' trains = ', trains
 print *,' bunches = ', bunches
 print *,' spacing = ', spacing  
       else
99        file = string(1:ix)
      endif
     endif
 print *,' file = ', file

   con%n_trains = trains
   con%n_cars = bunches

   con%n_14ns_space = spacing

   con%BunchPattern = file

   call bunchcross(lat,con,pc)

   ncross = pc%total_pc
   cross(1:ncross)  = pc%cross(1:ncross)%ele%s
  print *,' ncross = ',ncross
  l = ncross/10
  do i=1,l
   print '(10f8.2)',cross(10*(i-1)+1:min(i*10,ncross))
  enddo
  if(l*10 < ncross)print '(10f8.2)',cross(10*l+1:ncross)
 return
 end



  






