# Changelog

Log started 2024-01-01.

Types of entries:
- `Added` for new features.
- `Changed` for changes in existing functionality.
- `Deprecated` for soon-to-be removed features.
- `Removed` for now removed features.
- `Fixed` for any bug fixes.
- `Security` in case of vulnerabilities.


- 2024/01 Added: extraction line tracking to the long_term_tracking program.

- 2024/01 Added: *INDIVIDUAL* mode to the long_term_tracking program for resonant extraction simulations.

- 2024/01 Added: *foil* lattice element for simulating things like charge stripping, emittance smoothing, etc.



- 2024-03-

- 2024-03-

- 2024-03-

- 2024-03-

- 2024-03-

- 2024-03-06 Fixed: Tao now checks for call file infinite loop.

- 2024-03-06 Fixed: Fix radiation calc for taylor element with finite length.

- 2024-03-06 Fixed: Tao confusion with multiple universes with the same lattice and Rf is to be turned off.

- 2024-03-12 Fixed: long_term_tracking extraction tracking. 

- 2024-03-18 Added: output_only_last_turns and output_combined_bunches parameters to long_term_tracking. 

- 2024-03-20 Fixed: Particle z-calc with beam init and multiple bunches. 

- 2024-03-21 Added: ltt%print_info_messages parameter for long_term_tracking.

- 2024-03-22 Fixed: Test of particle outside of RF bucket ignoring closed orbit z.

- 2024-03-23 Added: SWAVE parameter to translation between Bmad and MAD8.

- 2024-03-25 Fixed: Added negative thickness error checking for foil element.

- 2024-03-25 Fixed: Now MPI messages will go through out_io for long_term_tracking.

- 2024-03-27 Fixed:  Reduced unnecessary output for tune_scan program.

- 2024-03-27 Added: "Bunch0" combined bunch for averages output in long_term_tracking.

- 2024-03-27 Added: Added ltt_com%ltt_tracking_happening_now in long_term_tracking for use in custom code.

- 2024-03-28 Fixed: Fix expression eval when there is an evaluation range in Tao.

- 2024-03-28 Fixed: Tao now checks for call file infinite loop.

- 2024-03-28 Fixed: Fix setting of lat_sigma_calc_needed logical in Tao.

- 2024-03-28 Fixed: m56 calc for standing_wave lcavity.

- 2024-03-29 Fixed: Corrected sbend changed attribute bookkeeping for k1 and k2. 

- 2024-04-06 Changed: Foil element attribute drel_thickness_dx changed to dthickness_dx.

- 2024-04-15 Added: New Tao command: `show rampers`.

- 2024-04-16 Fixed: Ramper bookkeeping when rampers control overlays and groups. 
The new bookkeeping is more efficient in terms of computation time.

- 2024-04-16 Fixed: Bug traced to OMP and multipole caching. Problem is that if auto_bookkeeper is
turned off initially, the cache is not use. Track1_bunch_csr turns it on for speed but
then tracking with OMP gives a race condition when the cache gets initialized during tracking.

