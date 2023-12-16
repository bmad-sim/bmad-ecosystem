# Bmad-Ecosystem repository
Bmad toolkit (library) for the simulation of charged particles and X-rays in accelerators and storage rings. This is the primary repository for the various libraries and programs that comprise the Bmad ecosystem. For details, see the Bmad website at [https://www.classe.cornell.edu/bmad/](https://www.classe.cornell.edu/bmad/).

## Manuals
+ [Bmad manual](https://www.classe.cornell.edu/bmad/manual.html)
+ [Tao manual](https://www.classe.cornell.edu/bmad/tao.html)
+ [Bmad & Tao tutorial](https://www.classe.cornell.edu/bmad/tao.html)
+ [Long_term_tracking program manual](https://www.classe.cornell.edu/bmad/other_manuals.html)
+ [Manuals for other Bmad-based programs](https://www.classe.cornell.edu/bmad/other_manuals.html)

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

If you want to compile Bmad directly, download a [Release](https://github.com/bmad-sim/bmad-ecosystem/releases) (or click on link on right hand side of this page and download the **bmad_dist.tar.gz** file (ignore the *source code* files) and follow the setup instructions at <https://wiki.classe.cornell.edu/ACC/ACL/OffsiteDoc>.


## Developer Setup

Developers should clone this repository, as well as its external dependencies:

```bash
git clone https://github.com/bmad-sim/bmad-ecosystem.git
git clone https://github.com/bmad-sim/bmad-external-deps.git
```

The external dependencies repository is simply a set of compressed files. A simple bash script is provided to extract these into the `bmad-ecosystem`:
```bash
cd bmad-ecosystem
bash util/extract_external_deps
```

To build everything:
```bash
cd bmad-ecosystem
source util/dist_source_me
util/dist_build_production
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
