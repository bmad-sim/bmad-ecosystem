    
    subroutine adjust_optics(lat)
    use bmad
    implicit none
    type (lat_struct) lat

      character(120) line, last_line
      character*16 ele_name
      integer ix
      integer ix_ele

      do while(line(1:1) /= 'G')
10    write (*, '(a)', advance = 'NO') '  element change or GO> '
      read(5, '(a)',err=10) line
     
  ix = index(line, '!')
  if (ix /= 0) line = line(:ix-1)        ! strip off comments

  call str_upcase(line, line)
  call string_trim(line, line, ix)

  if (ix == 0) then       ! nothing typed. do the same thing
      line = last_line
  endif

   last_line = line

   call str_upcase(line,line)
   ele_name = line(1:ix)
   if(line(1:1) /= 'G')then
    call find_change( line, lat)

    write (*, '(a)', advance = 'NO') ' Call type_ele ?'
    read(5,'(a)') line
    call str_upcase(line, line)
    call string_trim(line, line, ix)
     if(line(1:ix) == 'Y')then
      call element_locator(ele_name,lat,ix_ele)

      call type_ele(lat%ele(ix_ele),.true.,6)
     endif
    endif
   end do

   return
  end
