# Documentation for Regression Testing of Bmad

Two methods for running regression tests currently exist in bmad-ecosystem:

- A pytest-based test suite which
  - Automates running of separate Fortran binaries and comparing the results to
    `output.correct` files. Each of these tests are contained in a directory
    here `*_test` (for example, `pipe_test`).
  - Relies on the downstream [`PyTao`](https://github.com/bmad-sim/pytao/) library
    to run an additional suite of regression tests, verifying that example
    lattice simulation data has not changed between releases.
- The original `scripts/run_tests.py` which only supports automates running
  of separate Fortran binaries and comparing the results to `output.correct` files
  - This only automates running the separate Fortran binaries and does not support
    any PyTao-related tests.

# Overview of new Python and pytest-based tests

## Requirements

In order to install requirements to run these tests, use Python with `conda` or
`virtualenv` to create an isolated testing environment with the necessary
dependencies.

### With conda

```
conda create -n bmad-test python pip
conda activate bmad-test
python -m pip install -r regression_tests/requirements.txt
```

### With virtualenv

```
python -m venv test-env
source ./test-env/bin/activate
python -m pip install -r regression_tests/requirements.txt
```

## Running the tests

With the environment from the previous step activated (`conda activate` or `source venv/bin/activate`):

- Run all of the tests

  ```
  cd regression_tests
  pytest -v
  ```

- Run all of the tests with a specific binary directory:

  ```
  cd regression_tests
  pytest -v --bmad-bin=
  ```

- Run a test with `xyz` in the name:

  ```
  pytest -v -k xyz
  ```

- Run all tests and show standard output (i.e., `print`s) from each test:

  ```
  pytest -v -s
  ```

- Exit after the first failing test:

  ```
  pytest -v -x
  ```

- See more flags with `--help` to pytest.

  ```
  pytest --help
  ```

Consider using `python -m pytest` in place of `pytest` above if you are having
trouble using the correct `pytest` from your `PATH`.

## Implementation details

In `/regression_tests`, we have:

- `conftest.py`: this defines some directories and settings for all test suite tests.
  It allows you to configure additional `pytest` arguments, specifically `--bmad-bin`
  in order to check the results from different builds.
- `test_fortran.py`: this runs special test suite binaries defined in
  `./*_test` and compares the results based on `output.correct` files.
- `test_snapshots.py`: this uses PyTao to validate bmad-doc lattice examples
  against a snapshot version.

## `test_snapshots.py`

This uses PyTao to validate bmad-doc lattice examples against a snapshot
version.

To run these tests specifically,

```
pytest -v test_snapshots.py
```

### Updating example lattice snapshots

To update the example lattice snapshots, run:

```
python test_snapshots.py
```

Snapshots are stored in `/regression_tests/snapshots/{example_name}`.

# Overview of original testing method (`scripts/run_tests.py`)

The regression testing works by running programs and comparing the program
output to a file which has the "correct" output.

The main "regression_tests" directory contains a file called "TESTS.LIST". In
this file is a list of subdirectories of "regression_tests". These
subdirectories are called the "regression" subdirectories.

When the regression test script is run, each regression subdirectory is visited
in turn. In each regression subdirectory, the appropriate program is run. The
name of the program is assumed to be the same as the name of the regression
subdirectory. However, if the regression subdirectory has a file name "run.py",
this script is run instead of any program. If 'run.py' is executed, the
directory where the regression test script is looking for executable program
files is passed as the first (and only) argument to run.py.

The program will create an output file called "output.now". The regression
script will compare this file to an existing file in the subdirectory named
"output.correct". The regression test is successful if the numbers in two files
are the same within the tolerances given in the "output.now" file. See the
"Constructing a new test" section for details on the syntax of these files.

The results of the regression test is stored in a file named
"regression.results".

For the purposes of assigning an overall "passing" or "failing" grade to the
regression tests when the run_test.py script is run within a job script, for
each regression subdirectory there is a maximum allowable number of tests that
can be failed. This number is set in TESTS.LIST and the default is zero. If the
number of tests failed for any regression subdirectory is greater than the
maximum allowable for that subdirectory, the run_test.py script will exit at
the end with a non-zero exit code indicating overall regression failure.
Setting a non-zero maximum failure number should only be done as a temparary
measure to allow job scripts to run smoothly. Note: If run_test.py is run with
the -test option, it is assumed that the testing is not done within a script an
the maximum allowable failures is take to be zero since the pass/fail grade
does not matter in this case.

### Running the Tests

0. Do everything in the `regression_tests` directory.

1. If needed, compile and link the programs:

```
   make DO_EXTRA_MAKES=Y
   The executables are created in:
   ../bin/
```

2. Run the testing script `scripts/run_test.py`:

```
  scripts/run_test.py {-bin <exe_dir>} {-test <test_dir>} {-list <test_list_file>} {-debug}
```

Defaults:

```
<exe_dir> = "../bin" ! This is relative to current directory.
<test_dir> = "" ! For running a single test. Overrides using a test_list_file.
<test_list_file> = "test.list" ! For running multiple tests.
```

<exe_dir> is the directory where all the programs are. If <exe_dir> is a relative path name, it
must be relative to any subdirectory of regression_tests. <exe_dir> is optional and, if not
present, will default to "../../bin"

3. The results will be saved in a file "regression.results"

### Constructing a new test

1. Create a new subdirectory of regression_tests. For this documentation call this subdirectory
   "zzz".

2. Add a line to the "regression_tests/TESTS.LIST" file. Each line in this file is a subdirectory
   name indicating the subdirectory containing the test code. Note: Lines starting with an exclamation
   mark "!" are printed to the terminal but are otherwise ignored. After the subdirectory name, an
   optional number indicates the maximum number of tests that can be failed without triggering an error
   exit code at the end of all the tests (see above). For example, the number "4" would indicate that
   it is acceptable if the number of tests failed was four or less.

3. If needed: In this subdirectory put your testing program code and any input files needed for
   running the program. Additionally, a file "output.correct" needs to be present containing the
   expected output of the program. Generation of this file is detailed below.

4. When the program is run, it must generate a file called "output.now". This file will be compared
   to "output.correct".

5. The format of the "output.now" and "output.correct" files are as follows: Lines that
   are blank or that begin with a "!" exclamation character are ignored. All other
   "specification" lines have the syntax:
   "<tag_string>" <data_id> <tolerance> <value_1> <value_2> ...
   where:
   <tag_string> = Identification string associated with the line.
   <data_id> = How to interpret the data.
   <tolerance> = Allowable tolerance to pass test. Absent for STR data.
   <value_1>, etc. = List of values. Must use double quotes for STR data.
   Double quote marks must be used to surround the <tag_string> field.
   The <data_id> field can be one of:
   STR String or logical data.
   REL Real data with relative tolerance.
   VEC_REL Real data with tolerance relative to the vector magnitude of the array of numbers on the line
   ABS Real data with absolute tolerance.

The specification lines of "output.correct" and "output.now" must match up. That is, the N^th
specification line of "output.correct" and the N^th specification line in "output.now" must have the
same <tag_string> and the same number of values. Also the <tag_string> used for different lines
should be unique.

Example specification lines for "output.now":
"Orbit at s = 34.5" REL 1e-4 0.023 -0.0079 0.0060
"Element Names" STR "ABC cat" "XYZ dog"

In this example the "Orbit at s = 34.5" data has values 0.023, -0.0079, and 0.0060. The tolerance
is relative so to pass testing the corresponding values in "output.correct" must not differ by more
than 0.01%

5. CMake setup for compiling/linking the program:
   A) Add a `cmake.zzz` script in the regression_tests directory.
   Use, say, the existing `cmake.twiss_track_test` file as a template.
   In the cmake.zzz makefile substitute `zzz` for `twiss_track_test`.
   B) Edit the `CMakeLists.txt` in the regression_tests directory
   and add `cmake.zzz` to the EXE_SPECS list.

6. To generate the `output.correct` file:
   A) first run the program to generate the `output.now` file.  
   B) Copy `output.now` to `output.correct`.

7. You are now done. The scripts/run_tests.py script will parse the TESTS.LIST file in
   order to determine the list of programs to be run.

#### Guidelines

- Make the running time for programs short. Typical running times are less than
  a minute.

## Using Regression Programs for Custom Exploration

Sometimes it is convenient to be able to run a regression test with a modified
input file and for the test to output more than it would output normally. This
"custom exploration" mode has been built into a number of tests. For example,
the `mat6_calc_method_test` has a "custom exploration" mode that can be used to
compare the transfer matrix as computed with different tracking methods. In all
programs that have been configured for custom exploration, the custom
exploration is triggered if an input file is specified on the command line. If
custom exploration is triggered, a flag called print_extra is set to True to
signal custom output.
