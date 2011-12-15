! stub to allow ttychr to link on VMS
        integer function vms_ttychr
        character*1 ic
        print *,' non-blocking version later: now enter 1 char and <cr> '
        read(*,1) ic
1       format(a)
        vms_ttychr=ichar(ic)
        return
        end
