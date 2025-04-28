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

### Compile from source

If you want to compile Bmad directly,
download a [Release](https://github.com/bmad-sim/bmad-ecosystem/releases)
(or click on link on right hand side of this page and download the **bmad_dist.tar.gz** file.
ignore the _source code_ files))
and follow the setup instructions at <https://wiki.classe.cornell.edu/ACC/ACL/OffsiteDoc>.

## Developer Setup (for people involved in Bmad development)

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
