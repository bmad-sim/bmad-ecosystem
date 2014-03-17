!	call sparse_delim(ics) !ics is integer vec of chars to use as delim
!	should end with ics(i).lt.0  (ie 0 may be a delim)
!       eg To use space, tab, comma, null as delim, ics=[9,32,44,0,-1]
	subroutine sparse_delim(ics)
        use csr_sparse_mod
	implicit none
	integer:: ics(*),i
	do i=1,128
	 if(ics(i).lt.0) return 	!null=0 can be a delim
	 whitesp(i:i)=char(ics(i))  	!convert to char
	enddo
	return
	end			
