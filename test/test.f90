program test
  implicit none

  integer lun, lunget, i
  logical err
  character*7 fname

  do i=1,10
    lun = lunget()
    print *, "lun: ", lun
    write (fname,'(a,i2.2)') "test", i
    call openwr(lun, fname, err,0)
  end do

end program test
