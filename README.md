# bmad-ecosystem repository
Bmad toolkit (library) for the simulation of charged particles and X-rays in accelerators and storage rings.

**Bmad web site at: [https://www.classe.cornell.edu/bmad/](https://www.classe.cornell.edu/bmad/)**

## Bmad Setup

Two possibilities for setting up Bmad:
- Use **conda-forge** (a Python package manager). Instructions at <https://wiki.classe.cornell.edu/ACC/ACL/OffsiteDoc>. 

Or if you want to compile Bmad directly:
- Download a **Release** (click on link on right hand side of this page and download the **bmad_dist.tar.gz** file [Ignore the *source code* files.]) and follow the setup instructions at <https://wiki.classe.cornell.edu/ACC/ACL/OffsiteDoc>.


## Developer Setup

Developers should clone this repository, as well as its external dependencies:

```bash
git clone https://github.com/bmad-sim/bmad-ecosystem.git
git clone https://github.com/bmad-sim/bmad-external-deps.git
```

The external dependencies repository is simply a set of compressed files. A simple bash script is provided to extract these into the `bmad-ecoystem`:
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
