module bbu_track_mod

use bmad
use beam_mod

type bbu_stage_struct
  integer :: ix_ele_lr_wake   ! Element index of element with the wake 
  integer :: ix_ele_start     ! Exit end of this element is beginning of stage.
  integer :: ix_pass          ! Pass index when in multipass section
  integer :: ix_stage_pass1   ! Index of corresponding stage on first pass
  integer :: ix_head_bunch
  integer :: ix_hom_max
  integer :: ix_stage_to_track ! Ordered (in time) list of which stages to track. 
  real(rp) :: hom_power_max
end type

type bbu_beam_struct
  type (bunch_struct), allocatable :: bunch(:)  ! Bunches in the lattice
  integer ix_bunch_head       ! Index to head bunch(:)
  integer ix_bunch_end        ! Index of the end bunch(:)
  integer n_bunch_in_lat      ! Number of bunches transversing the lattice.
  integer ix_stage_power_max
  integer n_stage_to_track
  integer ix_last_stage_tracked
  type (bbu_stage_struct), allocatable :: stage(:)
  real(rp) hom_power_max
  real(rp) time_now
  real(rp) one_turn_time
  real(rp) dt_bunch
  real(rp) bunch_charge
end type

type bbu_param_struct
  character(80) lat_file_name
  logical hybridize
  logical write_hom_info
  logical keep_overlays_and_groups  ! Keep when hybridizing?
  logical keep_all_lcavities        ! Keep when hybridizing?
  real(rp) limit_factor
  real(rp) simulation_turns_max, bunch_freq
  real(rp) init_hom_amp
  real(rp) current
  real(rp) rel_tol
  logical drscan
  logical use_interpolated_threshold
  integer elindex
  character*40 elname
  integer nstep
  real(rp) begdr
  real(rp) enddr
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine bbu_setup (lat, bbu_param, bbu_beam)

use nr

implicit none

type (lat_struct), target :: lat
type (bbu_beam_struct), target :: bbu_beam
type (bbu_param_struct) bbu_param
type (ele_struct), pointer :: ele
type (lr_wake_struct), pointer :: lr
type (bunch_struct), pointer :: bunch

integer i, j, k, ib, ix_pass, n_links, ix_ele, n_bunch
integer, allocatable, save :: ix_chain(:)

real(rp), allocatable :: time_at_wake_ele(:)
real(rp) rr(4)

character(16) :: r_name = 'bbu_setup'

! Size bbu_beam%bunch

if (bbu_beam%dt_bunch == 0) then
  call out_io (s_fatal$, r_name, 'DT_BUNCH IS ZERO!')
  call err_exit
endif

bbu_beam%n_bunch_in_lat = (lat%ele(lat%n_ele_track)%ref_time / bbu_beam%dt_bunch) + 1
bbu_beam%one_turn_time = lat%ele(lat%n_ele_track)%ref_time

! Find all elements that have a lr wake.

j = 0
do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (.not. associated(ele%wake)) cycle
  if (size(ele%wake%lr) == 0) cycle
  j = j + 1
enddo

if (j == 0) then
  call out_io (s_fatal$, r_name, 'NO LR WAKES FOUND IN LATTICE!')
  call err_exit
endif

if (allocated(bbu_beam%stage)) deallocate (bbu_beam%stage)
allocate (bbu_beam%stage(j), time_at_wake_ele(j))

! Bunches have to go through a given physical lcavity in the correct order. 
! To do this correctly, the lattice is divided up into stages.
! A given stage has one and only one lcavity.
! The first stage starts at the beginning of the lattice.
! The last stage ends at the end of the lattice.
! All stages except the last end at an lcavity.

! bbu_beam%stage(i)%ix_ele_lr_wake holds the index in lat%ele(:) of the lcavity of the i^th stage.

! bbu_beam%stage(i)%ix_head_bunch holds the index in bbu_beam%bunch(:) of the head 
! bunch for the i^th stage.
! The head bunch is the next bunch that will be propagated through the stage.

j = 0
do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (.not. associated(ele%wake)) cycle
  if (size(ele%wake%lr) == 0) cycle
  j = j + 1
  bbu_beam%stage(j)%ix_ele_lr_wake = i
  time_at_wake_ele(j) = modulo(lat%ele(i)%ref_time, bbu_beam%dt_bunch)
  do k = 1, size(lat%ele(i)%wake%lr)
    lr => lat%ele(i)%wake%lr(k)
    lr%a_sin = bbu_param%init_hom_amp
    lr%a_cos = bbu_param%init_hom_amp
    lr%b_sin = bbu_param%init_hom_amp
    lr%b_cos = bbu_param%init_hom_amp
  enddo
  call multipass_chain (i, lat, ix_pass, n_links, ix_chain)
  bbu_beam%stage(j)%ix_pass = ix_pass
  if (ix_pass > 0) then
    do k = 1, j
      if (bbu_beam%stage(k)%ix_ele_lr_wake == ix_chain(1)) then
        bbu_beam%stage(j)%ix_stage_pass1 = k
        exit
      endif
    enddo
    if (k > j) then
      call out_io (s_fatal$, r_name, 'BOOKKEEPING ERROR.')
      call err_exit
    endif
  else
    bbu_beam%stage(j)%ix_stage_pass1 = j
  endif

  if (j == 1) then
    bbu_beam%stage(j)%ix_ele_start = 0
  else
    bbu_beam%stage(j)%ix_ele_start = bbu_beam%stage(j-1)%ix_ele_lr_wake
  endif

enddo

! which stage to track through next is determined by the time

call indexx (time_at_wake_ele, bbu_beam%stage%ix_stage_to_track)
bbu_beam%n_stage_to_track = 1
bbu_beam%time_now = 0

! bunch setup: Start at t = 0

if (allocated(bbu_beam%bunch)) deallocate (bbu_beam%bunch)
allocate(bbu_beam%bunch(bbu_beam%n_bunch_in_lat+10))

bbu_beam%stage%ix_head_bunch = -1    ! Indicates there are no bunches ready for a stage.

j = size(bbu_beam%stage)
ix_ele = bbu_beam%stage(j)%ix_ele_lr_wake
n_bunch = int(lat%ele(ix_ele)%ref_time / bbu_beam%dt_bunch) + 1

j = 1
do i = -1, -n_bunch, -1
  ib = i + 1 + n_bunch
  bunch => bbu_beam%bunch(ib)
  allocate (bunch%particle(1))
  bunch%particle(1)%charge = bbu_beam%bunch_charge
  bunch%ix_bunch = i
  bunch%t_center = (i+1) * bbu_beam%dt_bunch 
  do
    ix_ele = bbu_beam%stage(j)%ix_ele_lr_wake
    if (bunch%t_center + lat%ele(ix_ele)%ref_time > 0) exit
    j = j + 1
  enddo
  bunch%ix_ele = bbu_beam%stage(j)%ix_ele_start
  bbu_beam%stage(j)%ix_head_bunch = ib
enddo

bbu_beam%ix_bunch_head = 1
bbu_beam%ix_bunch_end = n_bunch 

end subroutine bbu_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine bbu_track_all (lat, bbu_beam, bbu_param, hom_power_gain, growth_rate, lost) 

implicit none

type (lat_struct) lat
type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param

real(rp) hom_power0, hom_power_sum, r_period, r_period0, power
real(rp) hom_power_gain, growth_rate

integer i, n_period, n_count, n_period_old, ix_ele

logical lost

character(16) :: r_name = 'bbu_track_all'

! Setup.

call bbu_setup (lat, bbu_param, bbu_beam)

call lattice_bookkeeper (lat)
bmad_com%auto_bookkeeper = .false. ! To speed things up.

! init time at hom and hom_power

bbu_beam%stage%hom_power_max = 0
bbu_beam%hom_power_max = 0
bbu_beam%ix_stage_power_max = 1
growth_rate = real_garbage$

! Track...
! To decide whether things are stable or unstable we look at the HOM power integrated
! over one turn period. The power is integrated over one period since there can 
! be a large variation in HOM power within a period.
! And since it is in the first period that all the stages are getting populated with bunches,
! We don't start integrating until the second period.

n_period_old = 0
hom_power_sum = 0
n_count = 0

do

  call bbu_track_a_stage (lat, bbu_beam, lost)
  if (lost) exit

  r_period = bbu_beam%time_now / bbu_beam%one_turn_time
  n_period = int(r_period)

  if (r_period > bbu_param%simulation_turns_max) then
    call out_io (s_warn$, r_name, 'Simulation_turns_max exceeded. ending tracking. \f10.2\ ', &
                                                                        r_array = (/ hom_power_gain /) )
    exit
  endif

  ! Compute integrated power. 
  ! The baseline hom power is computed over the 2nd period after the offset bunches 
  ! have passed through the lattice.

  if (n_period /= n_period_old) then
    if (n_period == 2) then
      hom_power0 = hom_power_sum / n_count
      r_period0 = r_period

    elseif (n_period > 2) then
      hom_power_gain = (hom_power_sum / n_count) / hom_power0
      growth_rate = log(hom_power_gain) / (r_period - r_period0)
      if (hom_power_gain < 1/bbu_param%limit_factor) exit
      if (hom_power_gain > bbu_param%limit_factor) exit      
    endif

    hom_power_sum = 0
    n_count = 0
    n_period_old = n_period
  endif

  call bbu_hom_power_calc (lat, bbu_beam)
  hom_power_sum = hom_power_sum + bbu_beam%hom_power_max
  n_count = n_count + 1
  
enddo

! Finalize

bmad_com%auto_bookkeeper = .true.


end subroutine bbu_track_all

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! This routine tracks one bunch through one stage.

subroutine bbu_track_a_stage (lat, bbu_beam, lost)

implicit none

type (lat_struct), target :: lat
type (bbu_beam_struct) bbu_beam
type (lr_wake_struct), pointer :: lr

real(rp) time_now

integer i, ix_stage, ix_bunch, ix_ele
integer ix_ele_start, ix_ele_end, j, ib, ib2

character(20) :: r_name = 'bbu_track_a_stage'

logical err, lost

! Look at each stage track the bunch with the earliest time to finish and track this stage.

ix_stage = bbu_beam%stage(bbu_beam%n_stage_to_track)%ix_stage_to_track
bbu_beam%n_stage_to_track = modulo (bbu_beam%n_stage_to_track, size(bbu_beam%stage)) + 1 

ib = bbu_beam%stage(ix_stage)%ix_head_bunch
if (ib == -1) then
  call out_io (s_fatal$, r_name, 'NO BUNCHES FOR THE TRACKING A STAGE. GET HELP!')
  call err_exit
endif

ix_ele_end = bbu_beam%stage(ix_stage)%ix_ele_lr_wake

time_now = bbu_beam%bunch(ib)%t_center + lat%ele(ix_ele_end)%ref_time
if (time_now < bbu_beam%time_now) then
  call out_io (s_fatal$, r_name, 'TIME ORDERING PROBLEM. GET HELP!')
  call err_exit
endif
bbu_beam%time_now = time_now

! Track

ix_ele_start = bbu_beam%bunch(ib)%ix_ele

do j = ix_ele_start+1, ix_ele_end
  call track1_bunch (bbu_beam%bunch(ib), lat, j, bbu_beam%bunch(ib), err)
  if (.not. all(bbu_beam%bunch(ib)%particle%ix_lost == not_lost$)) then
    lost = .true.
    return
  endif
enddo

! If the bunch tracked is the last bunch in the train, we need to add a new bunch for
! The next time this stage (which must have been stage #1) is tracked.

if (ib == bbu_beam%ix_bunch_end) then
  if (ix_stage /= 1) then
    call out_io (s_fatal$, r_name, 'TRACKING BOOKKEEPING PROBLEM. GET HELP!')
    call err_exit
  endif
  call bbu_add_a_bunch (lat, bbu_beam)
endif

! If this is the last stage then remove the bunch.
! Otherwiese, if the next stage does not have any bunches waiting to go 
! through then the tracked bunch becomes the head bunch for that stage.

if (ix_stage == size(bbu_beam%stage)) then  ! If last stage
  call bbu_remove_head_bunch (bbu_beam)
elseif (bbu_beam%stage(ix_stage+1)%ix_head_bunch == -1) then
  bbu_beam%stage(ix_stage+1)%ix_head_bunch = ib
endif

! If the bunch upstream from the tracked bunch is at the same stage as was tracked through,
! then this bunch becomes the new head bunch for the stage. Otherwise there is no head bunch
! for the stage.

ib2 = modulo (ib, size(bbu_beam%bunch)) + 1 ! Next bunch upstream
if (bbu_beam%bunch(ib2)%ix_ele == ix_ele_start) then
  bbu_beam%stage(ix_stage)%ix_head_bunch = ib2
else
  bbu_beam%stage(ix_stage)%ix_head_bunch = -1  ! No one waiting to go through this stage
endif

! Misc.

bbu_beam%ix_last_stage_tracked = ix_stage
lost = .false.

end subroutine bbu_track_a_stage

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine bbu_add_a_bunch (lat, bbu_beam)

implicit none

type (lat_struct) lat
type (bbu_beam_struct), target :: bbu_beam
type (bunch_struct), pointer :: bunch

integer i, ixb, ix0
real(rp) r(2)

character(20) :: r_name = 'bbu_add_a_bunch'

! Init bunch

ixb = modulo (bbu_beam%ix_bunch_end, size(bbu_beam%bunch)) + 1 ! Next empty slot
if (ixb == bbu_beam%ix_bunch_head) then
  call out_io (s_fatal$, r_name, 'BBU_BEAM%BUNCH ARRAY OVERFLOW')
  call err_exit
endif

bunch => bbu_beam%bunch(ixb)
if (.not. allocated(bunch%particle)) allocate (bunch%particle(1))
bunch%particle(1)%charge = bbu_beam%bunch_charge
bunch%ix_ele = 0

! If this is not the first bunch need to correct some of the bunch information

ix0 = bbu_beam%ix_bunch_end  ! Old bunch end
bunch%ix_bunch = bbu_beam%bunch(ix0)%ix_bunch + 1
bunch%t_center = bbu_beam%bunch(ix0)%t_center + bbu_beam%dt_bunch

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

if (bbu_beam%ix_bunch_end == bbu_beam%ix_bunch_head) then
  call out_io (s_fatal$, r_name, 'BOOKKEEPING PROBLEM. PLEASE GET HELP!')
endif

! Update ix_bunch_head 

bbu_beam%ix_bunch_head = modulo(bbu_beam%ix_bunch_head, size(bbu_beam%bunch)) + 1

end subroutine bbu_remove_head_bunch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Calculates power in mode with maximal amplitude

subroutine bbu_hom_power_calc (lat, bbu_beam)

implicit none

type (lat_struct), target :: lat
type (bbu_beam_struct) bbu_beam
type (lr_wake_struct), pointer :: lr

real(rp) hom_power_max, hom_power2

integer i, j, i1, ixm, ix

! Only need to update the last stage tracked

hom_power_max = 0
i = bbu_beam%ix_last_stage_tracked
i1 = bbu_beam%stage(i)%ix_stage_pass1
ix = bbu_beam%stage(i1)%ix_ele_lr_wake
do j = 1, size(lat%ele(ix)%wake%lr)
  lr => lat%ele(ix)%wake%lr(j)
  hom_power2 = max(lr%b_sin**2 + lr%b_cos**2, lr%a_sin**2 + lr%a_cos**2)
  if (hom_power_max < hom_power2) then
    hom_power_max = hom_power2
    bbu_beam%stage(i1)%ix_hom_max = j
  endif
enddo

! Update the new hom power.

hom_power_max = sqrt(hom_power_max)
bbu_beam%stage(i1)%hom_power_max = hom_power_max

! Find the maximum hom power in any element.

if (hom_power_max > bbu_beam%hom_power_max) then
  bbu_beam%ix_stage_power_max = i1
  bbu_beam%hom_power_max = hom_power_max

elseif (i1 == bbu_beam%ix_stage_power_max) then
  ixm = maxloc(bbu_beam%stage%hom_power_max, 1)
  bbu_beam%ix_stage_power_max = ixm
  bbu_beam%hom_power_max = bbu_beam%stage(ixm)%hom_power_max
endif

end subroutine bbu_hom_power_calc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine write_homs (lat, bunch_freq, trtb, currth)

! Write out information on lattice and HOMs
! Adapted from Changsheng's Get_Info.f90
!
! 19 March 2009 J.A.Crittenden

implicit none

type (bbu_param_struct) bbu_param

type (lat_struct) lat


!Arrays used to store lattice information
real(rp), Dimension(6,6) :: mat, oldmat,imat, testmat
real(rp), Dimension(:,:,:), allocatable :: erlmat
real(rp), Dimension(:),     allocatable :: erltime


integer i, j, k, kk, l, u
real(rp)  time, Vz, P0i, P0f, gamma
logical anavalid
logical judge
integer :: matrixsize = 4
real(rp) cnumerator,trtb,currth,currthc,rovq,matc,poltheta
real(rp) epsilon,kappa
real(rp) stest
integer nr
real(rp) bunch_freq
allocate(erlmat(800, matrixsize, matrixsize))
allocate(erltime(800))

!

! Initialize the identity matrix
call mat_make_unit(imat)

! Find the time between each cavity in the low energy lattice
! The time between two cavities is stored in erltime(k)

time=0
k=1     

do i = 1, lat%n_ele_track

  gamma=lat%ele(i)%value(E_TOT$)/m_electron       ! Calculate the z component of the velocity
  Vz=c_light*sqrt(1.-1./gamma**2)
   
  if (lat%ele(i)%key == LCAVITY$ ) then
    time=time+(lat%ele(i)%value(l$))/c_light*(1+0.5/(gamma*lat%ele(i-1)%value(E_TOT$)/m_electron))
    judge=.false.  ! True if the rf cavity has wake fields
    if (associated(lat%ele(i)%wake)) then
      do j=1, size(lat%ele(i)%wake%lr)
        if (lat%ele(i)%wake%lr(j)%R_over_Q > 1E-10) judge=.true.
      enddo
    endif
    if (judge) then
      print *, ' Storing time for HOM cavity',k,time
      erltime(k)=time
      time=0
      k=k+1       
    endif
  elseif (lat%ele(i)%key == hybrid$) then
    time = time + lat%ele(i)%value(delta_ref_time$)
  else
    time=time+(lat%ele(i)%value(l$)) / Vz
  endif
   
enddo



!Calculate Transport Matrices for the low energy lattice
!Matrix elements are stored in erlmat(i,j,k)
oldmat=imat
testmat = imat
P0i=lat%ele(0)%value(p0c$)     ! Get the first longitudinal reference momentum
k=1
kk=1
judge =.false.
anavalid = .true.

!print '(a)', ' Cavity    HOM      Ith(A)  Ith_coup(A)      tr       homfreq      RoverQ        Q       Pol Angle       T12         T14         T32        T34   sin omega*tr     tr/tb'
print '(a)', ' Cavity    HOM      Ith(A)  Ith_coup(A)      tr       homfreq  R/Q(ohms/m^2)     Q       Pol Angle       T12         T14         T32        T34   sin omega*tr     tr/tb'

do i=0, lat%n_ele_track
   
   mat=matmul(lat%ele(i)%mat6, oldmat) 
  
   if (lat%ele(i)%key == LCAVITY$ ) then
     if(associated(lat%ele(i)%wake)) then
       do j=1, size(lat%ele(i)%wake%lr)
          if(lat%ele(i)%wake%lr(j)%R_over_Q >1E-10) then            
            judge =.true.
          endif
       enddo
     endif

     if (judge) then

      P0f=lat%ele(i)%value(p0c$)
      mat(1,2)=mat(1,2)/P0i
      mat(1,4)=mat(1,4)/P0i
      mat(2,1)=mat(2,1)*P0f
      mat(2,2)=mat(2,2)*P0f/P0i
      mat(2,3)=mat(2,3)*P0f
      mat(2,4)=mat(2,4)*P0f/P0i
      mat(3,2)=mat(3,2)/P0i
      mat(3,4)=mat(3,4)/P0i
      mat(4,1)=mat(4,1)*P0f
      mat(4,2)=mat(4,2)*P0f/P0i
      mat(4,3)=mat(4,3)*P0f
      mat(4,4)=mat(4,4)*P0f/P0i

      do l=1, matrixsize
        do u=1, matrixsize
           erlmat(k,l,u)=mat(l,u)
        enddo
      enddo


      ! Don't print k=1, which is just the interval from the beginning
      ! of the lattice to the first cavity with a HOM

      if(k.gt.1.and.erltime(k).gt.0.)then

        do j=1, size(lat%ele(i)%wake%lr)

! Analytic approximation is not valid if any cavity has more than one HOM
           if (j.gt.1)anavalid=.false.

           ! This code uses the "linac definition" of R/Q, which is
           ! a factor of two larger than the "circuit definition." 
           ! The HOM files are in the "circuit definition" and
           ! the R/Q values are Ohms/m^2, whereas lat%ele(i)%wake%lr(j)%R_over_Q is in Ohms.
           rovq = 2*lat%ele(i)%wake%lr(j)%R_over_Q * (c_light/(2*pi*lat%ele(i)%wake%lr(j)%freq))**2
!           print *,' RovQ in Ohms',rovq

          ! Follow PRSTAB 7, 054401 (2004) (some bug here)
 !          kappa   = 2 * c_light * bunch_freq / (  rovq * (2*pi*lat%ele(i)%wake%lr(j)%freq)**2 )

           cnumerator = 2  * c_light / ( rovq * lat%ele(i)%wake%lr(j)%Q * 2*pi*lat%ele(i)%wake%lr(j)%freq )

           stest = mat(1,2) * sin ( 2*pi*lat%ele(i)%wake%lr(j)%freq*erltime(k) ) 

           epsilon =  2*pi*lat%ele(i)%wake%lr(j)%freq / ( bunch_freq * 2*lat%ele(i)%wake%lr(j)%Q )
           nr = int(erltime(k)*bunch_freq) + 1
           trtb = erltime(k)*bunch_freq


           if (nr*epsilon.lt.0.5)then
            if (stest.le.0.)then
! Threshold current for case epsilon * nr <<1 and T12*sin omega_lambda*tr < 0
              currth = - cnumerator / stest
            else
! Threshold current for case epsilon * nr <<1 and T12*sin omega_lambda*tr > 0
              currth = cnumerator / ( epsilon * abs(mat(1,2)) )
              if ( mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k), pi ) .le. pi/2 )then
               currth = currth * sqrt ( epsilon**2 + ( mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k), pi ) / nr )**2 )
!               print *,' First half. Currth, Mod =',currth,mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k), pi )
              else
               currth = currth * sqrt ( epsilon**2 + ( (pi - mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k), pi )) / nr )**2 )
!               print *,' Second half. Currth, Mod =',currth,mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k), pi )
              endif
            endif
           elseif (nr*epsilon.gt.2.)then
! Threshold current for case epsilon * nr >> 1
              currth = cnumerator / ( epsilon * abs(mat(1,2)) )
!
!              if ( mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) - sign(1.,stest)*pi/2, 2*pi )  .le. pi )then
!               currth = currth * sqrt ( epsilon**2 + ( mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) - sign(1.,stest)*pi/2, 2*pi ) / nr )**2 )
!               print *,' First half. Tr/tb, T12*sin, Currth, Mod =',trtb,stest,currth,mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) - sign(1.,stest)*pi/2, 2*pi )
!              else
!               currth = currth * sqrt ( epsilon**2 + ( (2*pi - mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) - sign(1.,stest)*pi/2, 2*pi )) / nr )**2 )
!               print *,' Second half. Tr/tb, T12*sin, Currth, Mod =',trtb,stest,currth,mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) - sign(1.,stest)*pi/2, 2*pi )
!              endif
!
              if ( mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) + pi/2, 2*pi )  .le. pi )then
               currth = currth * sqrt ( epsilon**2 + ( mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) + pi/2, 2*pi ) / nr )**2 )
               print *,' First half. Tr/tb, T12*sin, Currth, Mod =',trtb,stest,currth,mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) + pi/2, 2*pi )
              else
               currth = currth * sqrt ( epsilon**2 + ( (2*pi - mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) + pi/2, 2*pi )) / nr )**2 )
               print *,' Second half. Tr/tb, T12*sin, Currth, Mod =',trtb,stest,currth,mod( 2*pi*lat%ele(i)%wake%lr(j)%freq * erltime(k) + pi/2, 2*pi )
              endif
!!!!!!!!             currth = abs ( cnumerator / mat(1,2) )
           else
            currth=0.
          endif


          ! Threshold current for coupling
          poltheta = 2*pi*lat%ele(i)%wake%lr(j)%angle
          matc = mat(1,2)*cos(poltheta)**2 + ( mat(1,4) + mat(3,2) )*sin(poltheta)*cos(poltheta) + mat(3,4)*sin(poltheta)**2
          currthc = currth * abs ( mat(1,2) / matc )

          print '(i4, i9, 3x, 2es11.2, es12.5, 9es12.3, es14.5)', kk, j, currth, currthc, erltime(k), &
                           lat%ele(i)%wake%lr(j)%freq, lat%ele(i)%wake%lr(j)%R_over_Q,lat%ele(i)%wake%lr(j)%Q,lat%ele(i)%wake%lr(j)%angle, &
                           mat(1,2),mat(1,4),mat(3,2),mat(3,4), &
                           sin (2*pi*lat%ele(i)%wake%lr(j)%freq*erltime(k)),trtb

          print '(13x,a,es12.5,3x,a,es12.5,3x,a,es12.5)','R/Q= ', rovq, ' Ohms  epsilon= ', epsilon, 'nr*epsilon= ', nr*epsilon

          if(kk.ne.1)anavalid=.false.

        enddo
        kk=kk+1 ! Count nr of cavities for which analytic approximations are calculated

      endif   

      
      mat=imat     ! Re-initialize the transfer matrix
      k=k+1        ! Count cavities with HOMs
      P0i=P0f 

      endif ! End of judge selection
   endif ! End of cavity selection
   
   judge =.false.
   oldmat = mat
  
enddo

! Analytic approxmation is valid only for a single HOM in a single cavity
if (anavalid.eq..false.)then
 currth = 0.
endif

deallocate (erlmat, erltime)

end subroutine write_homs


end module
