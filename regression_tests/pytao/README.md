# Pytao Tests
This directory contains tests of Bmad using the pytao library as an interface.
The python module `pytest` is used to run tests and test discovery follows its conventions.

## Organization
- `lcavity`: Tests related to lcavity element
- `smoke_tests`: Add lattices here to confirm they at least load with `tao` without errors
- `sr_wakes`: Tests related to wake

## Filenames Note
Filenames for the pytao tests here should not be repeated due to the problem mentioned in PR #1620.
