program srdt_lsq_soln

use bmad
use srdt_mod

implicit none

character(200) in_file, lat_file

type(lat_struct) lat
type(coord_struct), allocatable :: co(:)
type(summation_rdt_struct) srdt
type(summation_rdt_struct), allocatable :: per_ele_rdt(:)
type(ele_pointer_struct), allocatable :: eles(:)

real(rp), allocatable :: ls_soln(:)
real(rp) chrom_x, chrom_y, weight(10)
real(rp) wgt_chrom_x, wgt_chrom_y, wgt_h20001, wgt_h00201, wgt_h10002
real(rp) wgt_h21000, wgt_h30000, wgt_h10110, wgt_h10020, wgt_h10200

integer i, j, k, nVar
integer status
integer gen_slices, sxt_slices
integer, allocatable :: var_indexes(:)

character(200) var_names

namelist /srdt_lsq/ lat_file, var_names, gen_slices, sxt_slices, chrom_x, chrom_y, &
                    wgt_chrom_x, wgt_chrom_y, wgt_h20001, wgt_h00201, wgt_h10002, &
                    wgt_h21000, wgt_h30000, wgt_h10110, wgt_h10020, wgt_h10200

gen_slices = 60
sxt_slices = 120
chrom_x = 1.0
chrom_y = 1.0
var_names = ''

in_file = 'srdt_lsq_soln.in'
if (iargc() > 0) call getarg(1, in_file)

wgt_chrom_x = 0;  wgt_chrom_y = 0;  wgt_h20001 = 0;  wgt_h00201 = 0;  wgt_h10002 = 0
wgt_h21000 = 0;  wgt_h30000 = 0;  wgt_h10110 = 0;  wgt_h10020 = 0;  wgt_h10200 = 0

open (unit = 10, file = in_file, action='read')
read (10, nml = srdt_lsq)
close (10)

write(*,*) "Preparing lattice..."

call bmad_parser(lat_file, lat)
call set_on_off(rfcavity$,lat,off$)
call twiss_and_track(lat, co, status)

call lat_ele_locator(var_names, lat, eles, nVar)
allocate(var_indexes(nVar))
do i = 1, nVar
  var_indexes(i) = eles(i)%ele%ix_ele
enddo

weight = [wgt_chrom_x, wgt_chrom_y, wgt_h20001, wgt_h00201, wgt_h10002, &
                        wgt_h21000, wgt_h30000, wgt_h10110, wgt_h10020, wgt_h10200]

call srdt_lsq_solution(lat, var_indexes, ls_soln, gen_slices, sxt_slices, chrom_x, chrom_y, weight)

open(45,file='srdt_lsq_soln.var')
do i=1,size(ls_soln)
  !write(45,'(a,a,es16.8)') trim(lat%ele(var_indexes(i))%name), '[k2]=', ls_soln(i)
  j = var_indexes(i)
  write(45,'(a,a,es16.8)') trim(lat%ele(j)%name), '[k2]=', ls_soln(i)
enddo
close(45)

end program


