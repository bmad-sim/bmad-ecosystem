!+
 ! Subroutine track1_gun_space_charge (bunch, ele, err, to_s_coords)
 !
 ! Subroutine to track a bunch of particles through an e_gun.
 !
 ! Input:
 !   bunch       -- bunch_struct: Starting bunch position.
 !   ele         -- ele_struct: E_gun element to track through. Must be part of a lattice.
 !   to_s_coords -- logical, optional: Default is True. If False leave bunch in time coords at end of tracking.
 !
 ! Output:
 !   bunch     -- bunch_struct: Ending bunch position.
 !   err       -- logical: Set true if there is an error. EG: Too many particles lost for a CSR calc.
 !-
 
 subroutine track1_bunch_space_charge (bunch, ele, err, to_s_coords)
 
 use space_charge_mod, dummy => track1_bunch_space_charge
 
 implicit none
 
 type (bunch_struct), target :: bunch
 type (ele_struct), target :: ele
 type (branch_struct), pointer :: branch
 type (coord_struct), pointer :: p
 
 integer i
 logical err, finished, include_image
 logical, optional :: to_s_coords
 real(rp) :: dt_step, dt_next, t_now, t_end
 
 integer, parameter :: fixed_time_step$ = 1, adaptive_step$ = 2 ! Need this in bmad_struct
 
 character(*), parameter :: r_name = 'track1_bunch_space_charge'
 
 ! Initialize variables
 
 branch => pointer_to_branch(ele)
 dt_step = space_charge_com%dt_track_step  ! Init time step.
 dt_next = dt_step
 include_image = (ele%space_charge_method == cathode_fft_3d$) ! Include cathode image charge?
 
 ! Don't track zero-length elements
 if (ele%value(l$) == 0) then
   if (logic_option(.true., to_s_coords)) call drift_to_s(bunch, ele%s, bunch)
   err = .false.
   return
 endif
 
 ! Drift bunch to the same time
 if (bunch%t0 == real_garbage$) then
   bunch%t0 = minval(bunch%particle%t, bunch%particle%state==alive$ .or. bunch%particle%state==pre_born$) 
 endif
 t_now = bunch%t0
 call drift_to_t(bunch, bunch%t0, bunch)
 
 ! Convert to t-based coordinate
 do i = 1, size(bunch%particle) 
   p => bunch%particle(i)
   call convert_particle_coordinates_s_to_t(p, s_body_calc(p, branch%lat%ele(p%ix_ele)), branch%lat%ele(p%ix_ele)%orientation)
 enddo
 
 ! Track
 do
   if (ele%tracking_method==fixed_step_time_runge_kutta$) then
     call sc_step(bunch, ele, include_image, t_now+dt_step)
   else
     call sc_adaptive_step(bunch, ele, include_image, t_now, dt_step, dt_next)
   end if
 
   t_now = t_now + dt_step
   dt_step = dt_next
 
   ! Check if all particles are finished
   finished = .true.
   do i = 1, size(bunch%particle)
     p => bunch%particle(i)
     if (p%state == pre_born$ .or. (p%s < ele%s - 0.1_rp * bmad_com%significant_length .and. p%state == alive$)) then
       finished = .false.
       exit
     endif
   enddo
   if (finished) exit
 enddo
 
 bunch%t0 = minval(bunch%particle%t, bunch%particle%state==alive$)
 do i= 1, size(bunch%particle) 
   p => bunch%particle(i)
   call convert_particle_coordinates_t_to_s(p, branch%lat%ele(p%ix_ele))
 enddo
 
 
 if (logic_option(.true., to_s_coords)) call drift_to_s(bunch, ele%s, bunch)
 err = .false.
 
 end subroutine track1_bunch_space_charge