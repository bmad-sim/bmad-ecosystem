program tune_scan

use ts_mod
!$ use omp_lib

implicit none

type (ts_params_struct) ts
type (ts_com_struct) ts_com
type (ts_data_struct), allocatable, target :: ts_dat(:,:,:)

integer ja, jb, jz, omp_threads

!---------------------------------------------
! Init

omp_threads = -1
!$ omp_threads = omp_get_max_threads()
if (omp_threads == -1) then
  print *, 'Note: tune_scan program compiled without OpenMP threading.'
else
  print *, 'OpenMP number of threads:', omp_threads
endif

call ts_init_params (ts, ts_com)
allocate (ts_dat(0:ts_com%n_a, 0:ts_com%n_b, 0:ts_com%n_z))

!---------------------------------------------
! Main loop

!$OMP PARALLEL DO COLLAPSE(3)
do ja = 0, ts_com%n_a
do jb = 0, ts_com%n_b
do jz = 0, ts_com%n_z
  call ts_track_particle (ts, ts_com, ja, jb, jz, ts_dat(ja,jb,jz))
enddo
enddo
enddo
!$OMP END PARALLEL DO

!---------------------------------------------
! Write results

call ts_write_results (ts, ts_com, ts_dat)

end program
