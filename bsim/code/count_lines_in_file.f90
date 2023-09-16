MODULE count_lines_in_file_mod

CONTAINS

SUBROUTINE count_lines_in_file(file_name,lines)

IMPLICIT NONE

!CHARACTER*200, INTENT(IN) :: file_name
CHARACTER(*), INTENT(IN) :: file_name
INTEGER, INTENT(OUT) :: lines

CHARACTER*266 bitbucket
INTEGER i, read_status

OPEN(100, FILE=file_name,STATUS='OLD',ACTION='READ')
i = 0
read_status = 0
DO WHILE (read_status .eq. 0)
  READ(100,FMT='(A)',IOSTAT=read_status) bitbucket
    IF (read_status .eq. 0) THEN
      i = i+1
    ENDIF
ENDDO
CLOSE(100)

lines = i

END SUBROUTINE count_lines_in_file

END MODULE count_lines_in_file_mod






