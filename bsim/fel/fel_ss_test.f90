!+
! Program fel_ss_test
!
! Single-slice steady-state FEL tracker, validated against Genesis 1.3 Version 4. This is
! deliverable 3 of the FEL port; see bsim/fel/README.md and the design brief.
!
! The program walks a Bmad lattice and applies the seam of the design (brief section 4.1):
!
!   - Elements whose name starts with UND are FEL segments. They are stepped internally in
!     delz with the transcribed Genesis physics (fel_track_mod): transverse push with
!     natural focusing, RK4 ponderomotive advance, source deposition and FFT field solve.
!     Bmad tracking is not used inside them.
!
!   - Every other element: the bunch is converted from the packed FEL arrays to
!     coord_structs and tracked by Bmad (track1_bunch). The packed arrays ARE Bmad
!     coordinates (see fel_beam_mod), so the conversion is a plain copy; the only phase
!     bookkeeping is one advance of the common reference phase phi0 per element. The
!     radiation field is drifted through free space by wavefront_drift.
!
! The starting state is a pair of Genesis dumps (&write of beam and field), so both codes
! track from bitwise-identical initial conditions. Diagnostics matching Genesis's
! definitions are recorded at the same z positions Genesis records them: once at the start
! and once after every integration step, one step per interlude element.
!
! Input is a namelist file:
!
!   &fel_ss_params
!     lat_file = "aramis.bmad"                 ! Bmad lattice.
!     beam_file = "Aramis-initial.par.h5"      ! Genesis particle dump to start from.
!     field_file = "Aramis-initial.fld.h5"     ! Genesis field dump to start from.
!     out_root = "fel_ss"                      ! Prefix for the three output files.
!     gamma0 = 11357.82                        ! Genesis's reference gamma.
!     delz = 0.045                             ! Target integration step inside undulators [m].
!     und_aw = 0.84853                         ! Undulator parameter (rms).
!     und_lambdau = 0.015                      ! Undulator period [m].
!     und_kx = 0.5, und_ky = 0.5               ! Natural focusing, deck convention (before ku^2).
!     und_helical = T
!     interlude_model = "bmad"                 ! "bmad" (the seam, default) or "genesis".
!     split_weights = F                        ! Weight-invariance test mode; see below.
!   &end
!
! interlude_model selects how the field-free elements are handled. "bmad" is the
! deliverable's architecture: track1_bunch for the particles, the exact theta mapping from
! Bmad's z, wavefront_drift for the field. "genesis" instead uses the transcribed Genesis
! interlude step (fel_track_interlude_genesis) everywhere, which prices what the seam
! changes: with it the whole run should agree with Genesis at transcription level, and the
! difference between the two modes is the transport model difference, measured rather than
! argued about.
!
! split_weights = T replaces each imported particle by two copies at identical
! coordinates carrying 1/3 and 2/3 of its weight. Every collective observable -- power,
! bunching, the field itself -- must be identical to the unsplit run, because the
! dynamics is per particle and the sources and reductions are linear in the weight. The
! benchmark harness runs this against the unsplit run to test the weighted paths, which
! nothing Genesis produces can test: the Genesis dump format carries no weights, so a
! Genesis comparison only ever sees the uniform case.
!
! Outputs: <out_root>.diag.txt (one row per record: z, field and beam diagnostics),
! <out_root>-final.fld.h5 and <out_root>-final.par.h5 (Genesis-format dumps of the end
! state, for field-by-field comparison).
!-

program fel_ss_test

use fel_track_mod
use wavefront_hdf5_mod
use beam_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (fel_beam_struct), target :: fbeam
type (fel_slice_struct), pointer :: sl
type (wavefront_struct) wf
type (bunch_struct) bunch
type (fel_und_struct) und
type (fel_slice_diag_struct) bdiag

real(rp) :: gamma0 = 0, delz = 0, und_aw = 0, und_lambdau = 0
real(rp) :: und_kx = 0.5_rp, und_ky = 0.5_rp
logical :: und_helical = .true.
logical :: split_weights = .false.
character(400) :: lat_file = '', beam_file = '', field_file = '', out_root = 'fel_ss'
character(16) :: interlude_model = 'bmad'

real(rp) unit_scale, z_now, ks, power, on_axis, qf
integer ie, istep, n_arg, iu_diag, iu_nml
logical err

character(400) param_file
character(*), parameter :: r_name = 'fel_ss_test'

namelist / fel_ss_params / lat_file, beam_file, field_file, out_root, gamma0, delz, &
                           und_aw, und_lambdau, und_kx, und_ky, und_helical, interlude_model, split_weights

! Read parameters.

n_arg = command_argument_count()
if (n_arg /= 1) then
  print '(a)', 'Usage: fel_ss_test <param_file>'
  stop 1
endif
call get_command_argument (1, param_file)

open (newunit = iu_nml, file = param_file, status = 'old', action = 'read')
read (iu_nml, nml = fel_ss_params)
close (iu_nml)

if (gamma0 <= 0 .or. delz <= 0 .or. und_aw <= 0 .or. und_lambdau <= 0) then
  print '(a)', 'fel_ss_test: gamma0, delz, und_aw and und_lambdau must all be set and positive.'
  stop 1
endif

if (interlude_model /= 'bmad' .and. interlude_model /= 'genesis') then
  print '(a)', 'fel_ss_test: interlude_model must be "bmad" or "genesis", got: ' // trim(interlude_model)
  stop 1
endif

! Read the lattice and the shared starting state.

call bmad_parser (lat_file, lat)
branch => lat%branch(0)

call fel_read_genesis4_beam (fbeam, beam_file, gamma0, err)
if (err) stop 1
if (split_weights) call do_split_weights (fbeam)
if (size(fbeam%slice) /= 1) then
  print '(a, i0)', 'fel_ss_test: steady state needs exactly one slice, got ', size(fbeam%slice)
  stop 1
endif
sl => fbeam%slice(1)

call wavefront_read_genesis4 (wf, field_file, err)
if (err) stop 1

! To Genesis internal units for the whole run; see fel_track_mod's header.

unit_scale = fel_field_unit_scale(wf)
wf%Ex = wf%Ex * unit_scale
ks = twopi / wf%wavelength

! Undulator segment parameters, constant for every segment in this benchmark. kx, ky get
! Genesis's unroll scaling by ku^2 (Lattice.cpp:412-413).

und%aw = und_aw
und%ku = twopi / und_lambdau
und%kx = und_kx * und%ku**2
und%ky = und_ky * und%ku**2
und%helical = und_helical

! Diagnostics file, one row per record at Genesis's record positions.

open (newunit = iu_diag, file = trim(out_root) // '.diag.txt', action = 'write')
write (iu_diag, '(a)') '#         z                  power         on_axis_intensity        bunching        ' // &
      'bunchingphase          energy          energyspread            xsize                 ysize'

z_now = 0
call write_diag_row()      ! Initial record, matching Genesis's diag before the first step.

! Walk the lattice.

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)

  ! Zero length elements (Bmad's end marker, for one) get no step and no diagnostic
  ! record: Genesis's unrolled lattice has no counterpart for them.

  if (ele%value(l$) == 0) cycle

  if (ele%name(1:3) == 'UND') then

    ! FEL segment: Genesis's unroll, nstep = round(l/delz), equal steps.

    und%nstep = nint(ele%value(l$) / delz)
    if (und%nstep == 0) und%nstep = 1
    und%dz = ele%value(l$) / und%nstep

    do istep = 1, und%nstep
      call fel_track_und_step (und, fbeam, wf, err)
      if (err) stop 1
      z_now = z_now + und%dz
      call write_diag_row()
    enddo

  elseif (interlude_model == 'bmad') then

    ! The seam: Bmad tracks the bunch (coordinate copies in and out), wavefront_drift
    ! moves the field, and the common phase phi0 advances by the reference rate with
    ! Genesis's drift surrogate ks/(2*gamma0^2) as the reference wavenumber.

    call fel_slice_to_bunch (fbeam, sl, ele, bunch, err)
    if (err) stop 1
    call track1_bunch (bunch, ele, err)
    if (err) then
      print '(2a)', 'fel_ss_test: tracking error in element ', trim(ele%name)
      stop 1
    endif
    call fel_bunch_to_slice (bunch, ele, sl, err)
    if (err) stop 1

    fbeam%phi0 = fbeam%phi0 + ele%value(l$) * &
                    fel_phi0_rate(ks, ks * 0.5_rp / gamma0 / gamma0, fbeam%p0_mc)

    call wavefront_drift (wf, ele%value(l$), err)
    if (err) stop 1

    z_now = z_now + ele%value(l$)
    call write_diag_row()

  else

    ! Genesis's own interlude model, transcribed, for pricing what the seam changes.

    qf = 0
    if (ele%key == quadrupole$) qf = ele%value(k1$)
    call fel_track_interlude_genesis (qf, ele%value(l$), fbeam, wf, err)
    if (err) stop 1

    z_now = z_now + ele%value(l$)
    call write_diag_row()
  endif
enddo

close (iu_diag)

! Final dumps in Genesis format. The field goes back to V/m for the wavefront writer;
! the writer's dfl conversion then lands within a couple of ulp of Genesis's own dump
! scale (the composition of the two constants is exact in real arithmetic).

wf%Ex = wf%Ex / unit_scale
call wavefront_write_genesis4 (wf, trim(out_root) // '-final.fld.h5', err, 'x')
if (err) stop 1

call fel_write_genesis4_beam (fbeam, trim(out_root) // '-final.par.h5', err)
if (err) stop 1

print '(a)', 'fel_ss_test done.'
print '(a)', '  ' // trim(out_root) // '.diag.txt'
print '(a)', '  ' // trim(out_root) // '-final.fld.h5'
print '(a)', '  ' // trim(out_root) // '-final.par.h5'

call wavefront_fft_free()

!------------------------------------------------------------------------------
contains

subroutine do_split_weights (beam)

! Replace each particle by two coincident copies with weights w/3 and 2w/3. The order --
! all first copies, then all second copies -- keeps the original particles' storage
! order, which keeps the RK4 arithmetic per copy identical to the unsplit run.

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sp
integer is, ip, n0

do is = 1, size(beam%slice)
  sp => beam%slice(is)
  n0 = sp%n
  call fel_slice_reallocate (sp, 2*n0)
  do ip = 1, n0
    sp%x(n0+ip) = sp%x(ip);  sp%px(n0+ip) = sp%px(ip)
    sp%y(n0+ip) = sp%y(ip);  sp%py(n0+ip) = sp%py(ip)
    sp%z(n0+ip) = sp%z(ip);  sp%pz(n0+ip) = sp%pz(ip)
    sp%weight(n0+ip) = 2 * sp%weight(ip) / 3
    sp%weight(ip) = sp%weight(ip) / 3
  enddo
  sp%n = 2*n0
enddo

end subroutine do_split_weights

!------------------------------------------------------------------------------

subroutine write_diag_row ()

call fel_field_diag (wf, power, on_axis)
call fel_slice_diag (fbeam, sl, ks, bdiag)

write (iu_diag, '(9es24.16)') z_now, power, on_axis, bdiag%bunching, bdiag%bunchingphase, &
                              bdiag%energy, bdiag%energyspread, bdiag%xsize, bdiag%ysize

end subroutine write_diag_row

end program fel_ss_test
