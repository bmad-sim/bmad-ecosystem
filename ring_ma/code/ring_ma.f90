program ring_ma
  use bmad
  use correct_ring_mod
  use dr_misalign_mod
!  use ran_state, only: ran_seed
  implicit none


  integer, parameter :: n_corrections_max = 4
  integer, parameter :: n_ma_max = 30

!==========================================================================
! Parameters for the input file

  character*200 lattice_file, output_file, comment
  real(rp) :: det_abs_res, det_diff_res, det_rot_res
  real(rp) :: alignment_multiplier=1.
  real(rp) :: sigma_cutoff=3.
  integer :: seed, n_iterations, n_lm_iterations, n_bad_seeds
  logical :: write_orbits=.false.
  real(rp) :: key_value1, key_value2
  type(ma_struct) :: ma(n_ma_max)
  type(correct_struct) :: correct(n_corrections_max)

  namelist /ring_ma_init/ lattice_file, seed, n_iterations, output_file, &
       comment, n_lm_iterations, correct, &
       write_orbits, key_value1, key_value2, ma, &
       alignment_multiplier, sigma_cutoff, &
       det_abs_res, det_diff_res, det_rot_res


  type data_struct
     real(rp) emit_x, emit_y, rms_x, rms_y, rms_eta_y, rms_cbar12, &
          rms_phi_x, rms_phi_y, rms_param
  end type data_struct
  type(data_struct), allocatable :: datablock(:,:)

  real(rp) rms_param

  type(lat_struct) :: design_ring        ! straight from the lattice file
  type(lat_struct) :: ma_ring            ! misaligned

  type(coord_struct), allocatable :: co(:)
  type(normal_modes_struct) modes
  integer rad_cache, i_ele, i_correction, ix_suffix, time(8), i_iter

  character*200 init_file, base_name, out_file, opt_file
  logical everything_ok

!==========================================================================
! Get init file from command line

  if (cesr_iargc() .ne. 1) then
     write(*,*) "usage: ring_ma <init_file>"
     stop
  end if
  call cesr_getarg(1, init_file)
  open(1, file=init_file)
  read(1, nml=ring_ma_init)
  close(1)

!==========================================================================
! Transfer the parameters from the input file to the appropriate structures

  dr_misalign_params%alignment_multiplier = alignment_multiplier
  dr_misalign_params%sigma_cutoff         = sigma_cutoff
  correct_ring_params%sigma_cutoff        = sigma_cutoff
  correct_ring_params%n_lm_iterations     = n_lm_iterations
  correct_ring_params%det_abs_res         = det_abs_res 
  correct_ring_params%det_diff_res        = det_diff_res
  correct_ring_params%det_rot_res         = det_rot_res

  dr_misalign_params%accumulate_errors    = .false.
  dr_misalign_params%tie_dup_ele          = .true.
  correct_ring_params%eta_delta_e_e       = 1.e-3
  correct_ring_params%write_elements      = .true.
  correct_ring_params%skip_dup_ele        = .true.

!==========================================================================
! Setup file names

  ix_suffix = index(init_file, ".", .true.)
  if (ix_suffix > 0) then
     base_name = init_file(:ix_suffix-1)
  else
     base_name = init_file
  end if
  out_file = trim(base_name) // ".out"

!==========================================================================
! Read in design ring from lattice and do some intialization


  if (match_wild(lattice_file, "*.xsif")) then
     call fullfilename(lattice_file, lattice_file)
     call xsif_parser(lattice_file, design_ring)
  else
     call bmad_parser(lattice_file, design_ring)
  end if

  do i_ele = 1, design_ring%n_ele_max
     if (design_ring%ele(i_ele)%key == wiggler$) &
          design_ring%ele(i_ele)%map_with_offsets = .false.
  end do
  call reallocate_coord(co, design_ring%n_ele_track)
  call lat_make_mat6(design_ring)
  call twiss_and_track(design_ring, co)
  if (write_orbits) then
     opt_file = trim(base_name)//".design"
     write(*,*) "Writing design orbit file: ", trim(opt_file)
     open(2, file=opt_file, recl=250)
     call write_opt(design_ring, co)
     close(2)
  end if

  allocate(datablock(n_iterations, 0:n_corrections_max))
  call reallocate_coord(cr_model_co, design_ring%n_ele_track)

!==========================================================================
! Initialize random-number generator from time or with specific value

  if (seed < 1) then
     call date_and_time(values=time)
     seed = time(8)
  end if
  call ran_seed(seed)

!==========================================================================
! Start iterations

  bmad_status%exit_on_error = .false.
  open(1,file=out_file, recl=250)
  write(1,'("#",A)') trim(comment)
  write(1,'("# Random seed:",i)') seed
  rad_cache = 0
  do i_iter = 1, n_iterations
     write(*,*) "============================================"
     write(*,*) "Starting iteration:        ", i_iter

     ! Catch when things are too far off
     n_bad_seeds = 0
     do
        ma_ring = design_ring
        call dr_misalign(ma_ring, ma)
        co(0)%vec = 0.
        call twiss_and_track (ma_ring, co, everything_ok)
        if (everything_ok) exit
        n_bad_seeds = n_bad_seeds + 1
        if (n_bad_seeds > 5) then
           write(*,'(A)') "Couldn't find any good seeds."
           write(*,'(A)') "Perhaps the misalignment parameters are too severe."
           stop
        end if
     end do

     call radiation_integrals(ma_ring, co, modes, rad_cache)
     call release_rad_int_cache(rad_cache)

     i_correction = 0
     if (write_orbits) then
        write(opt_file, '(A,".opt.",i3.3)') trim(base_name), i_iter
        write(*,*) "Writing optimization file: ", trim(opt_file)
        open(2, file=opt_file, recl=250)
        call write_opt(ma_ring, co)
     end if

     write(*,*) "Emittance:", modes%a%emittance, modes%b%emittance
 
     ! Store initial values
     call ring_to_data(ma_ring, co, modes, datablock(i_iter, 0))
     call write_data(datablock, i_iter, 0)

     ! Apply correction(s)
     do i_correction = 1, n_corrections_max
        rms_param = 0.
        if (all(correct(i_correction)%cor(:)%param == 0)) cycle
        cr_model_ring = design_ring
        call correct_ring(ma_ring, correct(i_correction), rms_param)
        call twiss_and_track(ma_ring, co)

        if (write_orbits) call write_opt(ma_ring, co)

        call radiation_integrals(ma_ring, co, modes, rad_cache)
        call release_rad_int_cache(rad_cache)
        call ring_to_data(ma_ring, co, modes, datablock(i_iter, i_correction))
        call write_data(datablock, i_iter, i_correction)
        write(*,*) "Emittance:", modes%a%emittance, modes%b%emittance
     end do
     if (write_orbits) close(2)
     correct_ring_params%write_elements = .false.
  end do

  ! Write summary
  call write_data(datablock)
  close(1)

contains
!==========================================================================
! Routine for writing out the complete element-by-element
! parameters for every seed. Can generate a LOT of output.

  subroutine write_opt(ring, co)
    implicit none
    type(lat_struct) :: ring
    type(coord_struct) :: co(:)
    integer i_ele
    real(rp) cbar(2,2)

    write(2,'(2A5,9A13,A18)') '# cor', 'ele', 's', 'x', 'y', 'eta_x', 'eta_y', &
         'phi_a', 'phi_b', 'cbar12', 'length', 'name'
    do i_ele = 1, ring%n_ele_track
       call c_to_cbar(ring%ele(i_ele), cbar)
       write(2,'(2i5,9e13.4,a18)') i_correction, i_ele, ring%ele(i_ele)%s, &
            co(i_ele)%vec(1), co(i_ele)%vec(3), &
            ring%ele(i_ele)%x%eta, ring%ele(i_ele)%y%eta, &
            ring%ele(i_ele)%a%phi - design_ring%ele(i_ele)%a%phi, &
            ring%ele(i_ele)%b%phi - design_ring%ele(i_ele)%b%phi, &
            cbar(1,2), ring%ele(i_ele)%value(l$), trim(ring%ele(i_ele)%name)
    end do
    write(2,*)
    write(2,*)
  end subroutine write_opt

!==========================================================================
! Routine to transfer relevant parameters from a ring to the DATABLOCK

  subroutine ring_to_data(ring, co, modes, data)
    implicit none
    type(lat_struct), intent(in) :: ring
    type(coord_struct), intent(in) :: co(:)
    type(normal_modes_struct), intent(in) :: modes
    type(data_struct), intent(out) :: data
    real(rp) l, ltot, cbar(2,2)
    integer i_ele

    ! Emittances
    data%emit_x = modes%a%emittance
    data%emit_y = modes%b%emittance

    ! RMS values
    data%rms_x = 0.
    data%rms_y = 0.
    data%rms_eta_y = 0.
    data%rms_cbar12 = 0.
    data%rms_phi_x = 0.
    data%rms_phi_y = 0.
    ltot = 0.
    do i_ele = 1, ring%n_ele_track
       call c_to_cbar(ring%ele(i_ele), cbar)
       l = ring%ele(i_ele)%value(l$)
       ltot = ltot + l
       data%rms_x      = data%rms_x      + l * co(i_ele)%vec(1)**2
       data%rms_y      = data%rms_y      + l * co(i_ele)%vec(3)**2
       data%rms_eta_y  = data%rms_eta_y  + l * ring%ele(i_ele)%y%eta**2
       data%rms_cbar12 = data%rms_cbar12 + l * cbar(1,2)**2
       data%rms_phi_x  = data%rms_phi_x  + l * (mod(ring%ele(i_ele)%a%phi - design_ring%ele(i_ele)%a%phi, twopi))**2
       data%rms_phi_y  = data%rms_phi_y  + l * (mod(ring%ele(i_ele)%b%phi - design_ring%ele(i_ele)%b%phi, twopi))**2
    end do
    data%rms_x      = sqrt(data%rms_x / ltot)
    data%rms_y      = sqrt(data%rms_y / ltot)
    data%rms_eta_y  = sqrt(data%rms_eta_y / ltot)
    data%rms_cbar12 = sqrt(data%rms_cbar12 / ltot)
    data%rms_phi_x  = sqrt(data%rms_phi_x / ltot)
    data%rms_phi_y  = sqrt(data%rms_phi_y / ltot)
    data%rms_param  = rms_param
  end subroutine ring_to_data

!==========================================================================
! Routine to write out the results of an individual seed, OR to write
! a summary table

  subroutine write_data(data, iter, cor)
    implicit none
    type(data_struct), allocatable, intent(in) :: data(:,:)
    integer, intent(in), optional :: iter, cor
    integer i_cor

    if (present(iter)) then
       if (iter==1 .and. cor==0) then
          write(1,'(2a6,7a14)') "# iter", "cor", "emit_y", "rms_y", &
               "rms_eta_y", "rms_cbar12", "rms_phi_x", "rms_phi_y", "param_rms"
       end if
       write(1,'(2i6,7e14.5)') iter, cor, &
            data(iter, cor)%emit_y, &
            data(iter, cor)%rms_y, &
            data(iter, cor)%rms_eta_y, &
            data(iter, cor)%rms_cbar12, &
            data(iter, cor)%rms_phi_x, &
            data(iter, cor)%rms_phi_y, &
            data(iter, cor)%rms_param
    else
       write(1,*)
       write(1,'(a)') "# Summary"
       write(1, '("#",a6,a10,6a14)') "cor", "param", "key_val1", "key_val2", "mean", "sigma", "50pct", "95pct"
       do i_cor = 0, n_corrections_max
          if (i_cor > 0) then
             if (all(correct(i_cor)%cor(:)%param == 0)) cycle
          end if
          write(1, '(" ",i6,a10,6e14.5)') i_cor, "emit_y",  key_value1, key_value2, data_line(data(:,i_cor)%emit_y)
          write(1, '(" ",i6,a10,6e14.5)') i_cor, "orbit_y", key_value1, key_value2, data_line(data(:,i_cor)%rms_y)
          write(1, '(" ",i6,a10,6e14.5)') i_cor, "eta_y",   key_value1, key_value2, data_line(data(:,i_cor)%rms_eta_y)
          write(1, '(" ",i6,a10,6e14.5)') i_cor, "cbar12",  key_value1, key_value2, data_line(data(:,i_cor)%rms_cbar12)
          write(1, '(" ",i6,a10,6e14.5)') i_cor, "phi_x",   key_value1, key_value2, data_line(data(:,i_cor)%rms_phi_x)
          write(1, '(" ",i6,a10,6e14.5)') i_cor, "phi_y",   key_value1, key_value2, data_line(data(:,i_cor)%rms_phi_y)
          write(1, '(" ",i6,a10,6e14.5)') i_cor, "param",   key_value1, key_value2, data_line(data(:,i_cor)%rms_param)
       end do
    end if
  end subroutine write_data

!==========================================================================
! Routine to compute statistical summary information for an array

  function data_line(array)
    use nr
    implicit none

    real(rp) :: data_line(4)
    real(rp), intent(in) :: array(:)
    integer indx(size(array)), k50, k95

    ! Compute mean and standard deviation
    call avevar(array, data_line(1), data_line(2))
    data_line(2) = sqrt(data_line(2))

    ! Compute quantiles if we have enough seeds
    if (size(array) > 10) then
       call indexx(array, indx)
       k50 = ceiling(.50 * size(array))
       k95 = ceiling(.95 * size(array))
       data_line(3) = array(indx(k50))
       data_line(4) = array(indx(k95))
    else
       data_line(3:4) = -1
    end if
  end function data_line

end program ring_ma
