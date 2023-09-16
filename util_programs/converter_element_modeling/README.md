This directory contains the two programs required to generate data for a Bmad converter element,
`converter_simulation` and `converter_fitter`.  To build them, see the documentation in the `doc`
folder, which details the required dependencies.

There is also a test program, `converter_test`, in the `test_prog` folder.  
The program uses data that has generated from geant4 and simulates particles going through a converter. 
The lattice used by this test program is named `lat.bmad` and an example lattice is in the `test_prog` directory.
To build the test program, copy the `test_prog` directory somewhere else on your system, and build like
any other ACC executable: ```cp -r test_prog /path/to/destination cd /path/to/destination/test_prog
mk ```


