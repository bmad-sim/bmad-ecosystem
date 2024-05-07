# Changelog

Log started 2024-01-01.

Types of entries:
- `Added` for new features.
- `Changed` for changes in existing functionality.
- `Deprecated` for soon-to-be removed features.
- `Removed` for now removed features.
- `Fixed` for any bug fixes.
- `Security` in case of vulnerabilities.

- 2024-05-06 Fixed Tao `cut` command orbit setting.

- 2024-04-24 Removed: srdt_lsq_soln program since it is not supported and Tao has this functionality.

- 2024-04-24 Fixed: Foil tracking when there is an offset.

- 2024-04-23 Changed: HomeBrew and MacPorts now automatically detected.

- 2024-04-19 Removed: Gprof references in compile scripts.

- 2024-04-19 Changed: Tweak to speed up offset_particle when there are no offsets.

- 2024-04-16 Fixed: `set_ele_attribute` routine will now set err_flag True on error with unknown variable in expression.

- 2024-04-16 Added: New Tao global: `global%beam_dead_cutoff`

- 2024-04-16 Fixed: Bug traced to OMP and multipole caching. Problem is that if auto_bookkeeper is
turned off initially, the cache is not use. Track1_bunch_csr turns it on for speed but
then tracking with OMP gives a race condition when the cache gets initialized during tracking.

- 2024-04-16 Fixed: Ramper bookkeeping when rampers control overlays and groups. 
The new bookkeeping is more efficient in terms of computation time.

- 2024-04-15 Added: New Tao command: `show rampers`.

- 2024-04-06 Changed: Foil element attribute drel_thickness_dx changed to dthickness_dx.

- 2024-03-29 Fixed: Corrected sbend changed attribute bookkeeping for k1 and k2. 

- 2024-03-28 Fixed: m56 calc for standing_wave lcavity.

- 2024-03-28 Fixed: Fix setting of lat_sigma_calc_needed logical in Tao.

- 2024-03-28 Fixed: Tao now checks for call file infinite loop.

- 2024-03-28 Fixed: Fix expression eval when there is an evaluation range in Tao.

- 2024-03-27 Added: Added `ltt_com%ltt_tracking_happening_now` in long_term_tracking for use in custom code.

- 2024-03-27 Added: "Bunch0" combined bunch for averages output in long_term_tracking.

- 2024-03-27 Fixed:  Reduced unnecessary output for tune_scan program.

- 2024-03-25 Fixed: Now MPI messages will go through out_io for long_term_tracking.

- 2024-03-25 Fixed: Added negative thickness error checking for foil element.

- 2024-03-23 Added: `SWAVE` parameter to translation between Bmad and MAD8.

- 2024-03-22 Fixed: Test of particle outside of RF bucket ignoring closed orbit z.

- 2024-03-21 Added: `ltt%print_info_messages` parameter for long_term_tracking.

- 2024-03-20 Fixed: Particle z-calc with beam init and multiple bunches. 

- 2024-03-18 Added: output_only_last_turns and output_combined_bunches parameters to long_term_tracking. 

- 2024-03-12 Fixed: long_term_tracking extraction tracking. 

- 2024-03-06 Fixed: Tao confusion with multiple universes with the same lattice and Rf is to be turned off.

- 2024-03-06 Fixed: Fix radiation calc for taylor element with finite length.

- 2024-03-06 Fixed: Tao now checks for call file infinite loop.

- 2024-02-26 Fixed: Correction to `track_a_lcavity` ref time calc.

- 2024-02-25 Fixed: Spin tracking will not respect element is_on = False setting.

- 2024-02-25 Fixed: Now chrom.w_a, etc. datums can be used with open lattice. 

- 2024-02-18 Fixed: Now long_term_tracking with energy ramping will properly init particles.

- 2024-02-16 Fixed: Eigen anal from sigma matrix.

- 2024-02-13 Fixed: Phase trombone tune set.

- 2024-02-12 Fixed: Reinstated phase_trombone in closed geometry lattice.

- 2024-02-13 Fixed: Add mode flip warning in Tao.

- 2024-02-13 Fixed: Phase trombone tune set. 

- 2024-02-12 Fixed: Updated tune_scan program to properly insert phase trombone element. 

- 2024-02-11 Added: `t_center` to `beam_init_struct`. 

- 2024-02-10 Added: Basic control_lord (for feedback elements) parsing is done. (#796)

- 2024-02-09 Added: New dispersion derivative. 

- 2024-02-09 Fixed: Fixed aperture_type set for super_slaves. (#792)

- 2024-02-05 Added: Added energy kick to beambeam tracking. (#766)

- 2024-01-18 Added: Code for pointing at cartesian_map(N)%term components.

- 2024-01-17 Fixed: Expression parsing of ...+d+....

- 2024-01-15 Fixed: Sliced crab_cavity phase calc.

- 2024-01-14 Added: *INDIVIDUAL* mode to the long_term_tracking program for resonant extraction simulations.

- 2024-01-11 Fixed: Linear beambeam spin tracking.

- 2024-01-11 Fixed: Corrected `ltt_init_tracking` logic for when beam init is needed.

- 2024-01-11 Fixed: SAD quad to bmad sad_mult translation.

- 2024-01-11 Added: *foil* lattice element for simulating things like charge stripping, emittance smoothing, etc.
