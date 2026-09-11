# Bmad-Ecosystem repository

Bmad toolkit (library) for the simulation of charged particles and X-rays in accelerators and storage rings. This is the primary repository for the various libraries and programs that comprise the Bmad ecosystem. For details, see the Bmad website at [https://www.classe.cornell.edu/bmad/](https://www.classe.cornell.edu/bmad/).

## Manuals

- [Bmad manual](https://www.classe.cornell.edu/bmad/manual.html)
- [Tao manual](https://www.classe.cornell.edu/bmad/tao.html)
- [Bmad & Tao tutorial](https://www.classe.cornell.edu/bmad/tao.html)
- [Long_term_tracking program manual](https://www.classe.cornell.edu/bmad/other_manuals.html)
- [Manuals for other Bmad-based programs](https://www.classe.cornell.edu/bmad/other_manuals.html)

## Bmad Installation

Bmad can be installed pre-compiled or from source. Detailed unstructions at <https://wiki.classe.cornell.edu/ACC/ACL/OffsiteDoc>.

### Pre-compiled from conda-forge

The simplest way to install Bmad is from [conda-forge](https://conda-forge.org). For the regular Bmad (OpenMP enabled), install the latest version using:

```zsh
conda install -c conda-forge bmad
```

For the MPI-enabled code, install the latest version using:

```zsh
conda install -c conda-forge bmad="*=mpi_openmpi*"
```

This will add all of the appropriate executables to the environment's PATH.

## Compile from Source

If you want to compile Bmad directly,
download a [Release](https://github.com/bmad-sim/bmad-ecosystem/releases)
(or click on link on right hand side of this page and download the **bmad_dist.tar.gz** file.
ignore the _source code_ files))
and follow the setup instructions at <https://wiki.classe.cornell.edu/ACC/ACL/OffsiteDoc>.

### Developer Setup (for people involved in Bmad development)

Developers should clone this repository, as well as the external packages repository:

```bash
git clone https://github.com/bmad-sim/bmad-ecosystem.git
git clone https://github.com/bmad-sim/bmad-external-packages.git
```

The external packages repository is simply a set of libraries needed by Bmad.

```bash
cd bmad-ecosystem
rm ../bmad-external-packages/README.md   # Do not copy this file
cp -r ../bmad-external-packages/* .
```

If this is the first time,
follow the setup instructions at <https://wiki.classe.cornell.edu/ACC/ACL/OffsiteDoc>.
Otherwise if the environment has been setup, to build do:

```bash
cd bmad-ecosystem
source util/dist_source_me
util/dist_build_production
```

### Conda-based development

In order to build Bmad without building each dependency one-by-one, you can use
conda to create an environment with all of the necessary build tools and
dependencies.

First, create a build environment:

```
conda env create -n bmad-build -f .github/bmad-build-env.yaml
conda activate bmad-build
```

This is the same environment used in GitHub Actions continuous integration.

Next, in `util/dist_prefs`:

1. Set `ACC_CONDA_BUILD` to `Y`
2. Set `ACC_CONDA_PATH` to `$CONDA_PREFIX`
3. Set `ACC_PLOT_PACKAGE` to `pgplot` (or `none` if desirable)
4. For PyTao usage, ensure `ACC_ENABLE_SHARED` is set to `Y` (if applicable)

Then:

```bash
source util/dist_source_me
util/dist_build_production
# or util/dist_build_debug
```

## Contributing to Bmad: Pull Requests

What is a Pull Request? A Pull Request (PR) is a mechanism for requesting that changes that you have made
to a copy of this repository (bmad-sim/bmad-ecosystem) are integrated (merged) into this repository.

The **main** branch of bmad-ecosystem is the central branch where all changes are merged into.

Pull Requests start with changes you have made to a branch that is not **main**. The PR is then a request for the changes you have made
to be merged with **main**.

Your copy of the bmad-ecosystem repository can be a
[fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks)
or simply a [clone](https://github.com/git-guides/git-clone).
Note: The procedure for
[creating a PR](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request)
when using a fork is somewhat different than when using a clone.

## Experimental: Direct CMake Build

We are in the process of replacing the old build system with a direct cmake build. Currently this should work to build tao and the bsim applications, as well as being able to write your own applications linked against the Bmad libraries. Only plplot plotting is supported at the moment. The build assumes that you have the dependent packages built separately or installed on your system. This: <https://github.com/bmad-sim/bmad-dependencies-lean> should build the packages that are not commonly installed/available on most systems, and assumes that fftw, gsl, hdf5, and lapack development libraries are already installed on your system.

A standard CMake build process involves a configuration step, followed by build and install steps. Furthermore, when you build with CMake, you build in a directory that is different from the one where the sources are. 

First, start with the configuration step. You should have the dependencies built already. Go into the top-level directory for bmad-ecosystem, and type
```
mkdir build
cd build
```
This creates the build directory and makes that the current directory. This directory could have been anywhere other than the source tree itself, so feel free to make a directory elsewhere. Now do the configuration step,
```
cmake -DCMAKE_INSTALL_PREFIX=$HOME/where/i/install/bmad -DCMAKE_PREFIX_PATH=$HOME/where/i/installed/the/packages ..
```
You can also add `-DCMAKE_BUILD_TYPE=Debug` to get a debug build, and `-DENABLE_OPENMP=ON` to enable OpenMP. The `..` at the end is really the path to the source tree, so if you want your build directory to be somewhere other than a subdirectory of the source tree, you need to replace `..` with the path to the source tree. Next, build and install:
```
cmake --build .
cmake --install .
```
To run tao, you simply need to add `$HOME/where/i/install/bmad/bin` to your path. If you are using pytao, you will need to add `$HOME/where/i/install/bmad/lib` to the `LD_LIBRARY_PATH` environment variable or otherwise tell pytao where the library is.

To build an application that uses the Bmad library, you simply build the application as you normally would with CMake, adding
```
find_package(Bmad REQUIRED)
target_link_libraries(my_program Bmad::bmad)
```
and adding `-DCMAKE_PREFIX_PATH=$HOME/where/i/install/bmad` to the CMake configure step. An example `CMakeLists.txt` for a simple executable would be:
```
cmake_minimum_required(VERSION 3.14)
project(my_project LANGUAGES Fortran)

find_package(Bmad REQUIRED)
add_executable(my_program my_program.f90)
target_link_libraries(my_program Bmad::bmad)
install(TARGETS my_program)
```
and you would configure with
```
cmake -DCMAKE_INSTALL_PREFIX=$HOME/where/my_program/goes -DCMAKE_PREFIX_PATH=$HOME/where/i/install/bmad ..
```
You can even avoid the install prefix for my_program if you like, skip the install step, and run the program out of the build tree.
