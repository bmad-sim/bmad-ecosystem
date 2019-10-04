!+
! subroutine set_tune3(lat, target_tunes, everything_ok, verbose_in, write_out_in, regex_in)
!
! wrapper for set_tune and set_z_tune together
!
! Modules Needed:
!    use bsim
!
! Input:
!   lat               -- lat_struct:
!   target_tunes(1:3) -- Real(rp): Integer + fractional tunes for a, b, z modes
!   verbose_in        -- logical, optional: Verbose output toggle; default == .false.
!   write_out_in      -- logical: Whether to write out quads after qtuneing. 0 == do not write,
!                                 > 0 means write out. Default == 0.
!   regex_in          -- character(40), optional: Regular expression mask for matching quads to
!                                 use in qtuneing. If regex string begins with 'match::', will
!                                 look for a match element to use rather than varying quads.
!   
!
! Output:
!   lat                          -- with adjusted quads and RF to match desired tunes
!     %ele(i_rf)%value(voltage$) -- Voltage on the cavity.
!   everything_ok                -- logical.  If present, returns true or false if set was
!                                             successful.  

!-
subroutine set_tune3(lat, target_tunes, everything_ok, verbose_in, write_out_in, regex_in)

  use bmad
  use z_tune_mod

  implicit none

  type(lat_struct), target :: lat
  type(ele_struct), pointer :: ele, match_ele
  type(ele_pointer_struct), allocatable :: match_eles(:)
  real(rp) target_tunes(3), new_tunes(3)
  type(coord_struct), allocatable :: co(:)
  real(rp), allocatable, save :: dk1(:), on_off_save(:)
  real(rp) :: dtilt = 1.e-2, eps = 5.e-4 ! convergence radius epsilon
  real(rp) dphi_a, dphi_b
  logical everything_ok, err
  logical verbose
  logical :: use_match_ele = .false.
  logical, optional :: verbose_in
  integer, optional :: write_out_in
  integer write_out, lun, param, ix, n_loc
  integer :: max_iter = 10
  type (all_pointer_struct) a_ptr
  
  character(40) :: param_name = '', write_name = '', regex = ''
  character(40), optional :: regex_in
  character(47) regex_match_check

  everything_ok = .true.

  ! set verbose logical; default to "true"
  if (present(verbose_in)) then
     verbose = verbose_in
  else
     verbose = .true.
  endif
  if (present(write_out_in)) then
     write_out = write_out_in
  else
     write_out = 0
  endif

  if (present(regex_in)) then
     regex = regex_in
  else
     regex = ''
  endif

  if (verbose .eqv. .true.) then
     write(*,'(a)') "============================================="
     write(*,'(a)') "Q_TUNING LATTICE..."
  endif

  
  if (all(target_tunes .lt. 1)) then
     write(*,'(a)') "Only fractional tunes given for target_tunes! Must supply integer + fractional tunes."
     write(*,'(a)') "Stopping here..."
     stop
  endif

  call set_on_off(rfcavity$, lat, save_state$, saved_values = on_off_save)
  call set_on_off(rfcavity$, lat, on$)
  call calc_z_tune(lat)
  call set_on_off(rfcavity$, lat, off$)

  if (verbose .eqv. .true.) write(*,'(a,3f10.5)') "INITIAL tunes: ", lat%ele(lat%n_ele_track)%a%phi/(2.*pi), &
       lat%ele(lat%n_ele_track)%b%phi/(2.*pi), lat%z%tune/(2.*pi)

  ! check if user provided regex corresponding to a match element:
  if (regex .ne. '') then
     write(regex_match_check, '(a7,a40)') 'match::',trim(regex)
     call lat_ele_locator(regex_match_check, lat, match_eles, n_loc, err)
     if (n_loc .gt. 1) then
        write(*,*) 'set_tune3 - too many match elements with user-supplied name. Using quads instead.'
     elseif (n_loc .eq. 1) then
        use_match_ele = .true.
        match_ele => match_eles(1)%ele
     endif
  endif

  

  ! if user has not specified one or more tunes, set target
  ! equal to the design. otherwise, translate fractional
  ! tune to radians.
  if (target_tunes(1) .lt. 1.e-12) target_tunes(1) = lat%ele(lat%n_ele_track)%a%phi / (2.*pi)
  if (target_tunes(2) .lt. 1.e-12) target_tunes(2) = lat%ele(lat%n_ele_track)%b%phi / (2.*pi)
  if (abs(target_tunes(3)) .lt. 1.e-12) target_tunes(3) = lat%z%tune / (2.*pi)
  
  new_tunes(1:3) = target_tunes(1:3) * 2.*pi


  
  if (use_match_ele) then ! set dphi.a and dphi.b for match element:
     do ix = 1, max_iter ! allow max_iter retries
          call pointer_to_attribute(match_ele, 'dphi_a', .true., a_ptr, err)
          a_ptr%r = a_ptr%r + (new_tunes(1) - lat%ele(lat%n_ele_track)%a%phi)
          call set_flags_for_changed_attribute(match_ele, a_ptr)
     
          call pointer_to_attribute(match_ele, 'dphi_b', .true., a_ptr, err)
          a_ptr%r = a_ptr%r + (new_tunes(2) - lat%ele(lat%n_ele_track)%b%phi)
          call set_flags_for_changed_attribute(match_ele, a_ptr)
     
          call lattice_bookkeeper(lat)
     
          call closed_orbit_calc(lat, co)
          call lat_make_mat6(lat, -1, co)
          call twiss_at_start(lat)
          call twiss_propagate_all(lat)
          
          dphi_a = abs(new_tunes(1) - lat%ele(lat%n_ele_track)%a%phi)
          dphi_b = abs(new_tunes(2) - lat%ele(lat%n_ele_track)%b%phi)
          if ((dphi_a .lt. eps) .and. (dphi_b .lt. eps)) then
               exit
          elseif (ix .eq. max_iter) then
               everything_ok = .false.
               exit
          else
               continue           
          endif

     enddo

  else ! using array of quads, and not a single match element
     
     if (.not. allocated(dk1)) then ! first call to set_tune3; allocate and find quads
        allocate(dk1(lat%n_ele_max))
        call choose_quads_for_set_tune(lat, dk1, regex)
     endif

     call set_tune(new_tunes(1), new_tunes(2), dk1, lat, co, everything_ok)

  endif
  
  call set_on_off(rfcavity$, lat, on$)
  call set_z_tune(lat, new_tunes(3))
  call set_on_off(rfcavity$, lat, off$)

  if (verbose .eqv. .true.) then
     write(*,'(a,3f10.5)') "DESIRED tunes: ", target_tunes(:)
     write(*,'(a,3f10.5)') "    SET tunes: ", lat%ele(lat%n_ele_track)%a%phi/(2.*pi), &
          lat%ele(lat%n_ele_track)%b%phi/(2.*pi), lat%z%tune/(2.*pi)
     write(*,'(a)')        "============================================="
  endif
  
  if (everything_ok .eqv. .false.) then
     if (global_com%exit_on_error .eqv. .true.) then 
        write(*,*) "Could not set tunes! Stopping here..."
        stop
     else
        write(*,*) "Could not set tunes! Continuing..."
     endif
  endif

  if (write_out .gt. 0) then ! write out elements used to qtune:
     lun = lunget()
     open(unit=lun, file='qtune3.bmad', status='replace')
     do ix = 1, lat%n_ele_max
        ele => lat%ele(ix)
        if (.not. ((ele%key .eq. quadrupole$) .or. (ele%key .eq. rfcavity$))) cycle 

        if (ele%key .eq. quadrupole$) then
           if (abs(abs(value_of_attribute(ele,'TILT',err)) - pi/4.) .lt. dtilt) cycle ! don't write out skew quads; 
                                        ! can't just look for pi/4 because of misaligned skews
           param = k1$
           param_name = 'K1'
        elseif (ele%key .eq. rfcavity$) then
           param = voltage$
           param_name = 'VOLT'
        endif

        if (ele%slave_status .eq. super_slave$) cycle 
   
        ! remove 'by_index' vs. 'by_name', in favor of using '<ele_name>##n' indexing.
        ! note: requires first calling assign_unique_ele_ids! 
        if (ele%ix_pointer .eq. 0) then
           write_name = ele%name
        else
           write(write_name,'(2a,i0)') trim(ele%name), '##', ele%ix_pointer
        endif
        write(lun,'(4a,es18.9)') trim(write_name), '[', trim(param_name), &
             '] := ', value_of_attribute(ele,param_name,err)
        
     enddo
     close(lun)
  endif ! write_out .gt. 0

  call set_on_off(rfcavity$, lat, restore_state$, saved_values = on_off_save)


end subroutine set_tune3
