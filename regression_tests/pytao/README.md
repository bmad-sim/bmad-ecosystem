# Pytao Tests
This directory contains tests of Bmad using the pytao library as an interface.
The python module `pytest` is used to run tests and test discovery follows its conventions.

## Organization
- `bookkeeping`: Tests related to lattice bookkeeping status flags
- `change_cmd`: Tests of the `change` command switches and `@`, `d`, `%` prefixes
- `completion`: Tests of the `pipe complete` tab completion command
- `kicks`: Small kick and multipole response of every element type compared with the first order formulas
- `lcavity`: Tests related to lcavity element
- `set_data_var`: Tests of `set data` / `set var` and the settable flags in `pipe data` / `pipe var`
- `smoke_tests`: Add lattices here to confirm they at least load with `tao` without errors
- `sr_wakes`: Tests related to wake

## Filenames Note
Filenames for the pytao tests here should not be repeated due to the problem mentioned in PR #1620.
