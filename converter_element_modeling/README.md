This directory contains the two programs required to generate a Bmad converter element, `converter_simulation` and `converter_fitter`.
To build them, see the documentation in the `doc` folder, which details the required dependencies.
We have also included a test program, `converter_test`, in the `test_prog` folder.
This doesn't require any special dependencies to build, but it is useless without a converter element file.
To build the test program, copy the `test_prog` directory somewhere else on your system, and build like any other ACC executable:
```
cp -r test_prog /path/to/destination
cd /path/to/destination/test_prog
mk
```
