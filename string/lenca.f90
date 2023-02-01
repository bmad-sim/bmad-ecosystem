!       n_notblank=lenca(strvec,len) !return # of last non-blank string
!       from vector of strings. Note if first str is blank or null, no others
!       are checked, but otherwise strvec is checked from last to first,
!       so results with imbedded null/blank can be inconsistent.
	integer*4 function lenca(ca,lll)
	implicit none
	integer*4:: lll,i
	character*(*):: ca(lll)
	lenca=0       !if 1st ele is blank, or null, return 0
	if(ca(1).eq.' '.or.ichar(ca(1)(1:1)).eq.0) return 
	do i=lll,1,-1
         if((ca(i).ne.' ').and.(ichar(ca(i)(1:1)).ne.0)) then
          lenca=i ; return
         endif 
	enddo
	return
	end
