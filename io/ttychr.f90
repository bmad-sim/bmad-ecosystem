         integer function ttychr()
         use sim_utils    !for get_a_char
         use input_mod
         !wrapper to get 1 char, non-blocking  
         implicit none
         character*1 ic,ignore_this(2)
#if defined (CESR_VMS) 
          ttychr=vms_ttychr()
#else
         ignore_this(1)=' '
         call get_a_char(ic,.false.,ignore_this)  !dont wait, ignore blank
         ttychr=ichar(ic)
#endif
         return
         end
