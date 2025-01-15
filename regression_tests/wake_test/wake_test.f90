program wake_test

use beam_mod
use bmad

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, ele_z1, ele_p1, ele_p2
type (wake_sr_mode_struct), pointer :: w
type (coord_struct) :: orb0
type (coord_struct) :: p1, p2
type (beam_init_struct) :: beam_init
type (bunch_struct) :: bunch, bunch_init, bunch0, bunch2
type (bunch_params_struct) bparams

real(rp) dz, z0

integer i, nargs, n, ie, i0, i1

logical print_extra, err_flag

character(100) :: lat_file

!------------------------------------------

lat_file = 'wake_test.bmad'

print_extra = .false.
nargs = command_argument_count()
if (nargs == 1)then
   call get_command_argument(1, lat_file)
   print *, 'Using ', trim(lat_file)
   print_extra = .true.
elseif (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

call bmad_parser (lat_file, lat)
!call write_bmad_lattice_file('lat2.bmad', lat)
!call bmad_parser('lat2.bmad', lat)

call ran_seed_put(123456)

open (1, file = 'output.now')

!---------------------------------
! Transfer test

ele_z1 => lat%branch(1)%ele(1)   ! z1 element
ele_p1 => lat%ele(1)             ! p1 element
ele_p2 => lat%ele(2)             ! p2 element

write (1, '(a, 2l1, a)') '"Wake-Transfer" STR "', ele_z1%wake%sr == ele_p1%wake%sr, ele_z1%wake%lr == ele_p2%wake%lr, '"'

! Short range wake test.

ele => lat%ele(1)
write (1, '(a, 2i4, f10.2)') '"SR-Size" ABS 0' , &
                size(ele%wake%sr%long), size(ele%wake%sr%trans), ele%wake%sr%z_max

do i = 1, size(ele%wake%sr%trans)
  w => ele%wake%sr%trans(i)
  write (1, '(a, i0, a, 4es14.6, 2i4)') '"SR-T', i, '" ABS 0', &
                                w%amp, w%damp, w%k, twopi*w%phi, w%polarization, w%position_dependence
enddo

do i = 1, size(ele%wake%sr%long)
  w => ele%wake%sr%long(i)
  write (1, '(a, i0, a, 4es14.6, i4)') '"SR-L', i, '" ABS 0', &
                               w%amp, w%damp, w%k, twopi*w%phi, w%position_dependence
enddo

!

z0 = 7e-3_rp

p1%vec = 0
p1%vec(1) = 1e-3
p1%vec(3) = -2e-3
p1%vec(5) = z0

call init_coord (p1, p1%vec, lat%ele(0), element_end = upstream_end$)
call init_coord (p2, p1%vec, lat%ele(0), element_end = upstream_end$)

p1%charge = 1
p2%charge = 0

ele => lat%ele(1)

! Make wake from first particle

ele%wake%sr%long%b_sin = 0
ele%wake%sr%long%b_cos = 0
ele%wake%sr%long%a_sin = 0
ele%wake%sr%long%a_cos = 0
ele%wake%sr%z_ref_long = p1%vec(5)

ele%wake%sr%trans%b_sin = 0
ele%wake%sr%trans%b_cos = 0
ele%wake%sr%trans%a_sin = 0
ele%wake%sr%trans%a_cos = 0
ele%wake%sr%z_ref_trans = p1%vec(5)

call sr_longitudinal_wake_particle (ele, p1)
call sr_transverse_wake_particle (ele, p1)

! Add in a second wake

p1%vec(5) = z0 + lat%particle_start%vec(5)
call sr_longitudinal_wake_particle (ele, p1)
call sr_transverse_wake_particle (ele, p1)

!

dz = -2e-3

do i = 1, 2
  p2%vec = 0
  p2%vec(5) = i*dz + z0
  call sr_longitudinal_wake_particle (ele, p2)
  call sr_transverse_wake_particle (ele, p2)
  write (1, '(a, i0, a, 3es20.9)') '"SR', i, '" REL 1E-8' , p2%vec(2:6:2)
enddo

! Beam tracking

beam_init%n_particle = 100
beam_init%random_engine = 'quasi'
beam_init%a_norm_emit = 3e-9
beam_init%b_norm_emit = 1e-9
beam_init%dPz_dz = 0.0
beam_init%bunch_charge = 1 !100.0e-12
beam_init%sig_pz = 1e-12
beam_init%sig_z = 5.99585e-3  ! 200 ps * cLight

call init_bunch_distribution (lat%ele(0), lat%param, beam_init, 0, bunch_init)
bunch = bunch_init

call calc_bunch_params(bunch_init, bparams, err_flag)
write (1, '(a, 6es16.8)') '"BP-Charge" ABS 1E-8', bparams%charge_live, bparams%charge_tot
write (1, '(a, 6es16.8)') '"BP-ST" ABS 1E-8', bparams%s, bparams%t
write (1, '(a, 6es16.8)') '"BP-Centroid" ABS 1E-8', bparams%centroid%vec
write (1, '(a, 6es16.8)') '"BP-Sig1" ABS 1E-8', bparams%sigma(1,:)
write (1, '(a, 6es16.8)') '"BP-Sig2" ABS 1E-8', bparams%sigma(2,:)
write (1, '(a, 6es16.8)') '"BP-Sig3" ABS 1E-8', bparams%sigma(3,:)
write (1, '(a, 6es16.8)') '"BP-Sig4" ABS 1E-8', bparams%sigma(4,:)
write (1, '(a, 6es16.8)') '"BP-Sig5" ABS 1E-8', bparams%sigma(5,:)
write (1, '(a, 6es16.8)') '"BP-Sig6" ABS 1E-8', bparams%sigma(6,:)
write (1, '(a, 6es16.8)') '"BP-Amode" ABS 1E-8', bparams%a%beta, bparams%a%alpha, bparams%a%emit, bparams%a%norm_emit
write (1, '(a, 6es16.8)') '"BP-Bmode" ABS 1E-8', bparams%b%beta, bparams%b%alpha, bparams%b%emit, bparams%b%norm_emit
write (1, '(a, 6es16.8)') '"BP-Xmode" ABS 1E-8', bparams%x%beta, bparams%x%alpha, bparams%x%emit, bparams%x%norm_emit
write (1, '(a, 6es16.8)') '"BP-Ymode" ABS 1E-8', bparams%y%beta, bparams%y%alpha, bparams%y%emit, bparams%y%norm_emit
write (1, '(a, 6es16.8)') '"BP-Zmode" ABS 1E-8', bparams%z%beta, bparams%z%alpha, bparams%z%emit, bparams%z%norm_emit

!

call track1_bunch (bunch, ele, err_flag)

bmad_com%sr_wakes_on = .false.
bunch0 = bunch_init
call track1_bunch (bunch0, ele, err_flag)

write (1, '(a, 6es18.9)') '"SR-P20" REL 1E-8', bunch%particle(20)%vec - bunch0%particle(20)%vec
write (1, '(a, 6es18.9)') '"SR-P40" REL 1E-8', bunch%particle(40)%vec - bunch0%particle(40)%vec

! Long range wake test

ele => lat%ele(2)

bunch = bunch_init

do n = 1, 3
  call track1_bunch (bunch, ele, err_flag)
  do i = 1, size(bunch%particle)
    bunch%particle(i)%t = bunch%particle(i)%t + 1e-7
  enddo
enddo

bmad_com%lr_wakes_on = .false.

bunch0 = bunch_init

do n = 1, 3
  call track1_bunch (bunch0, ele, err_flag)
  do i = 1, size(bunch0%particle)
    bunch0%particle(i)%t = bunch0%particle(i)%t + 1e-7
  enddo
enddo

write (1, '(a, 6es18.9)') '"LR-P20" ABS 1E-19' , bunch%particle(20)%vec - bunch0%particle(20)%vec
write (1, '(a, 6es18.9)') '"LR-P40" ABS 1E-19' , bunch%particle(40)%vec - bunch0%particle(40)%vec

! Long range wake with superimposed element.

bmad_com%lr_wakes_on = .true.
ele => lat%ele(3)

bunch2 = bunch_init

do n = 1, 3
  do ie = 3, 5
    call track1_bunch (bunch2, lat%ele(ie), err_flag)
  enddo
  do i = 1, size(bunch2%particle)
    bunch2%particle(i)%t = bunch2%particle(i)%t + 1e-7
  enddo
enddo

write (1, '(a, 6es18.9)') '"dB-LR-Pipe-P20" ABS 1E-20' , dvec(bunch2%particle(20)%vec - bunch%particle(20)%vec)
write (1, '(a, 6es18.9)') '"dB-LR-Pipe-P40" ABS 1E-20' , dvec(bunch2%particle(40)%vec - bunch%particle(40)%vec)

! Long range wake in rf cavity

bunch0 = bunch_init
do n = 1, 3
  call track1_bunch (bunch0, lat%ele(6), err_flag)
enddo

bunch2 = bunch_init
do n = 1, 3
  call track1_bunch (bunch2, lat%ele(7), err_flag)
enddo

bunch2%particle%vec(5) = bunch2%particle%vec(5) - 6d-17       ! Correction due to rf2 having finite length
bunch2%particle%vec(6) = bunch2%particle%vec(6) + 6.9414d-12  ! Correction due to rf2 having finite length

write (1, '(a, 6es18.9)') '"dB-LR-RF-P20" ABS 1E-15' , dvec(bunch2%particle(20)%vec - bunch0%particle(20)%vec)
write (1, '(a, 6es18.9)') '"dB-LR-RF-P40" ABS 1E-15' , dvec(bunch2%particle(40)%vec - bunch0%particle(40)%vec)

!---------------------------------
! Sort test

bunch%particle(1)%state = lost$
bunch%particle(3)%state = lost$
bunch%particle(10)%state = lost$
call order_particles_in_z (bunch)

do i = 1, size(bunch%particle) - 1
  i0 = bunch%ix_z(i)
  i1 = bunch%ix_z(i+1)

  if (bunch%particle(i1)%state /= alive$) then
    if (i /= bunch%n_live) then
      print *, 'Sort problem1'
      call err_exit
    endif
    exit
  endif

  if (bunch%particle(i0)%vec(5) < bunch%particle(i1)%vec(5)) then
    print *, 'Sort problem2'
    call err_exit
  endif
enddo

!------------------------------------------------------------------
! SR z_long wake test

ele => lat%branch(2)%ele(1)
bunch = bunch_init
bmad_com%sr_wakes_on = .true.
call track1_bunch (bunch, ele, err_flag)

write (1, '(a, 6es18.9)') '"SRZ-A" REL 1E-8', bunch%particle(1:5)%vec(6) - bunch_init%particle(1:5)%vec(6)
write (1, '(a, 6es18.9)') '"SRZ-B" REL 1E-8', bunch%particle(21:25)%vec(6) - bunch_init%particle(21:25)%vec(6)

!------------------------------------------------------------------
contains

! vec(6) is much larger than the other components so rescale to be of the same magnitude

function dvec(vec_in) result (vec_out)
real(rp) vec_in(6), vec_out(6)
vec_out = vec_in
vec_out(6) = 1d-4 * vec_in(6)
end function dvec

end program
