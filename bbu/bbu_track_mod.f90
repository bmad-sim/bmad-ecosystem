module bbu_track_mod

use bmad
use beam_mod

type bbu_stage_struct
  integer :: ix_ele_lr_wake 
  integer :: ix_head_bunch
end type

type bbu_beam_struct
  type (bunch_struct), allocatable :: bunch(:)  ! Bunches in the lattice
  integer, allocatable :: ix_ele_bunch(:)       ! element where bunch is 
  integer ix_bunch_head       ! Index to head bunch(:)
  integer ix_bunch_end        ! Index of the end bunch(:). -1 -> no bunches.
  integer n_bunch_in_lat      ! Number of bunches transversing the lattice.
  type (bbu_stage_struct), allocatable :: stage(:)
end type

type bbu_param_struct
  character(80) lat_file_name
  logical hyberdize
  real(rp) limit_factor
  real(rp) low_power_lim, high_power_lim
  real(rp) simulation_time, bunch_freq, init_hom_amp
  real(rp) current
  real(rp) rel_tol
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine bbu_setup (lat, ds_bunch, init_hom_amp, bbu_beam)

implicit none

type (lat_struct), target :: lat
type (bbu_beam_struct) bbu_beam
type (beam_init_struct) beam_init
type (ele_struct), pointer :: ele

integer i, j, ih

real(rp) ds_bunch, init_hom_amp, rr(4)

character(16) :: r_name = 'bbu_setup'

! Size bbu_beam%bunch

if (ds_bunch == 0) then
  call out_io (s_fatal$, r_name, 'DS_BUNCH IS ZERO!')
  call err_exit
endif

bbu_beam%n_bunch_in_lat = (lat%param%total_length / ds_bunch) + 1

if (allocated(bbu_beam%bunch)) deallocate (bbu_beam%bunch)
allocate(bbu_beam%bunch(bbu_beam%n_bunch_in_lat+1))
bbu_beam%ix_bunch_head = 1
bbu_beam%ix_bunch_end = -1  ! Indicates No bunches

call re_allocate (bbu_beam%ix_ele_bunch, bbu_beam%n_bunch_in_lat+1)

! Find all elements that have a lr wake.

j = 0
do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (.not. associated(ele%wake)) cycle
  if (size(ele%wake%lr) == 0) cycle
  j = j + 1
  do ih = 1, size(ele%wake%lr)
    call ran_gauss (rr)
    ele%wake%lr%b_sin = init_hom_amp * rr(1) / 2
    ele%wake%lr%b_cos = init_hom_amp * rr(2) / 2
    ele%wake%lr%a_sin = init_hom_amp * rr(3) / 2
    ele%wake%lr%a_cos = init_hom_amp * rr(4) / 2
    ele%wake%lr%b_sin = 0
    ele%wake%lr%b_cos = init_hom_amp 
    ele%wake%lr%a_sin = 0
    ele%wake%lr%a_cos = 0
  enddo
enddo

if (j == 0) then
  call out_io (s_fatal$, r_name, 'NO LR WAKES FOUND IN LATTICE!')
  call err_exit
endif

if (allocated(bbu_beam%stage)) deallocate (bbu_beam%stage)
allocate (bbu_beam%stage(j))

! Bunches have to go through a given physical lcavity in the correct order. 
! To do this correctly, the lattice is divided up into stages.
! A given stage has one and only one lcavity.
! The first stage starts at the beginning of the lattice.
! The last stage ends at the end of the lattice.
! All stages except the last end at an lcavity.

! bbu_beam%stage(i)%ix_ele_lr_wake holds the index in lat%ele(:) of the lcavity of the i^th stage.

! bbu_beam%stage(i)%ix_head_bunch holds the index in bbu_beam%bunch(:) of the head bunch for the i^th stage.
! The head bunch is the next bunch that will be propagated through the stage.

bbu_beam%stage%ix_head_bunch = -1    ! Indicates there are no bunches ready for a stage.
bbu_beam%stage(1)%ix_head_bunch = 1

j = 0
do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (.not. associated(ele%wake)) cycle
  if (size(ele%wake%lr) == 0) cycle
  j = j + 1
  bbu_beam%stage(j)%ix_ele_lr_wake = i
enddo

end subroutine bbu_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! This routine tracks one bunch through one stage.

subroutine bbu_track_a_stage (lat, bbu_beam, lost)

implicit none

type (lat_struct) lat
type (bbu_beam_struct) bbu_beam

real(rp) best_finish_time, finish_time

integer i_stage, i_stage_track, ix_bunch, ix_ele
integer ix_ele_start, ix_ele_end, j, ib, ib2

character(20) :: r_name = 'bbu_track_a_stage'

logical err, lost

! Look at each stage track the bunch with the earliest time to finish.

best_finish_time = 1e30 ! Something large
i_stage_track = -1

do i_stage = 1, size(bbu_beam%stage)

  ! ix_bunch is the bunch index of the next bunch waiting to go through.
  ! If no bunch is waiting to do this stage then cycle to the next stage

  ix_bunch = bbu_beam%stage(i_stage)%ix_head_bunch
  if (ix_bunch < 0) cycle

  ! ix_ele is the index of the element at the end of the stage.

  ix_ele = bbu_beam%stage(i_stage)%ix_ele_lr_wake
  finish_time = bbu_beam%bunch(ix_bunch)%t_center + lat%ele(ix_ele)%ref_time

  if (finish_time < best_finish_time) then
    best_finish_time = finish_time
    i_stage_track = i_stage
  endif

enddo

if (i_stage_track == -1) then
  call out_io (s_fatal$, r_name, 'NO BUNCHES TO PROPAGATE!')
  call err_exit
endif

! The bunch is now cleared for tracking through this stage

ix_ele_end = bbu_beam%stage(i_stage_track)%ix_ele_lr_wake
if (i_stage_track == size(bbu_beam%stage)) ix_ele_end = lat%n_ele_track

ib = bbu_beam%stage(i_stage_track)%ix_head_bunch
ix_ele_start = bbu_beam%bunch(ib)%ix_ele

do j = ix_ele_start+1, ix_ele_end
  call track1_bunch (bbu_beam%bunch(ib), lat, j, bbu_beam%bunch(ib), err)
  if (.not. all(bbu_beam%bunch(ib)%particle%ix_lost == not_lost$)) then
    lost = .true.
    return
  endif
enddo

! If the next stage does not have any bunches waiting to go through then the
! tracked bunch becomes the head bunch for that stage.

if (i_stage_track /= size(bbu_beam%stage)) then  ! If not last stage
  if (bbu_beam%stage(i_stage_track+1)%ix_head_bunch == -1) &
                      bbu_beam%stage(i_stage_track+1)%ix_head_bunch = ib
endif

! If the bunch upstream from the tracked bunch is at the same stage as was tracked through,
! then this bunch becomes the new head bunch for the stage. Otherwise there is no head bunch
! for the stage.

if (ib /= bbu_beam%ix_bunch_end) then
  ib2 = modulo (ib, size(bbu_beam%bunch)) + 1 ! Next bunch upstream
  if (bbu_beam%bunch(ib2)%ix_ele == ix_ele_start) then
    bbu_beam%stage(i_stage_track)%ix_head_bunch = ib2
  else
    bbu_beam%stage(i_stage_track)%ix_head_bunch = -1  ! No one waiting to go through this stage
  endif
endif

lost = .false.

end subroutine bbu_track_a_stage

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine bbu_add_a_bunch (lat, bbu_beam, beam_init)

implicit none

type (lat_struct) lat
type (bbu_beam_struct) bbu_beam
type (beam_init_struct) beam_init

integer ixb, ix0

character(20) :: r_name = 'bbu_add_a_bunch'

! Init bunch

if (bbu_beam%ix_bunch_end == -1) then ! if no bunches
  ixb = bbu_beam%ix_bunch_head
else
  ixb = modulo (bbu_beam%ix_bunch_end, size(bbu_beam%bunch)) + 1 ! Next empty slot
  if (ixb == bbu_beam%ix_bunch_head) then
    call out_io (s_fatal$, r_name, 'BBU_BEAM%BUNCH ARRAY OVERFLOW')
    call err_exit
  endif
endif

call init_bunch_distribution (lat%ele(0), beam_init, bbu_beam%bunch(ixb))

! If this is not the first bunch need to correct some of the bunch information

if (ixb /= bbu_beam%ix_bunch_head) then
  ix0 = bbu_beam%ix_bunch_end
  bbu_beam%bunch(ixb)%ix_bunch = bbu_beam%bunch(ix0)%ix_bunch + 1
  bbu_beam%bunch(ixb)%z_center = bbu_beam%bunch(ix0)%z_center - beam_init%ds_bunch
  bbu_beam%bunch(ixb)%t_center = -bbu_beam%bunch(ixb)%z_center * &
                         lat%ele(0)%value(p0c$) / (c_light * lat%ele(0)%value(e_tot$))
endif

bbu_beam%ix_bunch_end = ixb

end subroutine bbu_add_a_bunch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine bbu_remove_head_bunch (bbu_beam)

implicit none

type (bbu_beam_struct) bbu_beam
character(20) :: r_name = 'bbu_remove_head_bunch'

! If there are no bunches then there is an error

if (bbu_beam%ix_bunch_end == -1) then
  call out_io (s_fatal$, r_name, 'TRYING TO REMOVE NON-EXISTANT BUNCH!')
  call err_exit
endif

! mark ix_bunch_end if we are poping the last bunch.

if (bbu_beam%ix_bunch_end == bbu_beam%ix_bunch_head) bbu_beam%ix_bunch_end = -1

! Update ix_bunch_head 

bbu_beam%ix_bunch_head = modulo(bbu_beam%ix_bunch_head, size(bbu_beam%bunch)) + 1

end subroutine bbu_remove_head_bunch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine bbu_hom_power_calc (lat, bbu_beam, hom_power)

implicit none

type (lat_struct), target :: lat
type (bbu_beam_struct) bbu_beam
type (lr_wake_struct), pointer :: lr

real(rp) hom_power

integer i, j, ix, ix_pass

!

hom_power = 0
n_hom = 0

do i = 1, size(bbu_beam%stage)
  ix = bbu_beam%stage(i)%ix_ele_lr_wake
  call multipass_chain (i, lat, ix_pass)
  if (ix_pass > 1) cycle
  do j = 1, size(lat%ele(ix)%wake%lr)
    lr => lat%ele(ix)%wake%lr(j)
    hom_power = hom_power + lr%b_sin**2 + lr%b_cos**2  + lr%a_sin**2 + lr%a_cos**2
  enddo
enddo

hom_power = sqrt(hom_power)

end subroutine bbu_hom_power_calc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_power, lost) 

implicit none

type (lat_struct) lat
type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (beam_init_struct) beam_init

real(rp) hom_power

integer i

logical lost

! Setup.

call bbu_setup (lat, beam_init%ds_bunch, bbu_param%init_hom_amp, bbu_beam)

do i = 1, size(bbu_beam%bunch)
  call bbu_add_a_bunch (lat, bbu_beam, beam_init)
enddo

! Track

do

  call bbu_track_a_stage (lat, bbu_beam, lost)
  if (lost) return

  ! If the head bunch is finished then remove it and seed a new one.

  if (bbu_beam%bunch(bbu_beam%ix_bunch_head)%ix_ele == lat%n_ele_track) then
    call bbu_remove_head_bunch (bbu_beam)
    call bbu_add_a_bunch (lat, bbu_beam, beam_init)
    if (bbu_beam%bunch(bbu_beam%ix_bunch_end)%t_center > bbu_param%simulation_time) return
  endif

  ! Compute average power

  call bbu_hom_power_calc (lat, bbu_beam, hom_power)

  if (hom_power < bbu_param%low_power_lim) return
  if (hom_power > bbu_param%high_power_lim) return

enddo

end subroutine

end module
