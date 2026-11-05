# ADR-0078: The Windows CI jobs share one cache of conda package archives

**Status:** Accepted
**Date:** 2026-11-05

## Context

GitHub keeps 10 GB of Actions caches per repository and evicts the least recently used beyond that.  On 2026-11-04 the repository held 8.9 GB in 23 caches, all on develop, and no Linux cache at all.  Each Windows job cached its conda package folder under its own key (`windows-conda-cpp-...` and `windows-conda-<python>-...`), about 2.6 GB apiece: six of them do not fit, and while they took turns evicting each other, they evicted every other platform's caches and every pull request's.  A PR restores only its own caches and develop's, so PRs compiled cold.

Three things made the folders large and many:

- **What the environment installs.**  `environment-win.yml` listed the `compilers` metapackage, whose Fortran compiler (flang) brings llvm-tools, lld and clang: 830 MB of a 1153 MB download.  The builds use only clang-cl from it.  conda's git (115 MB) duplicated the runner's, and numpy and scipy took MKL (109 MB) as their BLAS.
- **Unpacked copies.**  A conda package folder holds each package twice, the `.conda` archive and its unpacked folder, so 1.15 GB of downloads made a 2.6 GB cache.
- **One cache per Python version, keyed by the workflow file.**  The solves for 3.10 and 3.14 share 282 MB of 320 MB; only python, numpy, scipy, eigenpy, boost-python and mpi4py differ (~40 MB).  And `hashFiles(..., '.github/workflows/build_and_test.yml')` in the key made every workflow edit a new set of caches, while the old set waited to be evicted.

## Decision

- `environment-win.yml` asks for `clang >=20` (clang-cl) instead of `compilers`, drops git, and pins numpy's BLAS to OpenBLAS (`libblas=*=*openblas`).  The download is ~320 MB.  The Clang floor is the one Visual Studio 2026's STL demands (STL1000), so the builds no longer define `_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH`.
- Every Windows job uses one cache, `windows-conda-archives-<hash of environment-win.yml>-<yyyy-MM>`, of the package archives only (`conda_pkgs/*.conda`, `conda_pkgs/*.tar.bz2`); conda unpacks an archive it finds.  The C++ test job and the wheel jobs restore it and, on a miss, save it; whichever saves first wins, and the others download their few Python-specific packages.  The wheel tests restore it and save nothing.
- The key names the environment file and the month, not the workflow, so a workflow edit keeps the cache and a new month brings in conda-forge's updates (a hit saves nothing, so without the month new package versions would be downloaded on every run for good).

## Consequences

- Do not go back to a cache per Python version or per job of the whole package folder: at ~2.6 GB apiece they crowd the other platforms out of the 10 GB.
- Do not put the workflow file back into the conda key; put the environment file, and anything else that changes what conda installs, there instead.
- Do not reinstate `compilers` to get clang-cl; ask for `clang`.  If the STL raises its Clang floor, raise the pin.
- A package added to `environment-win.yml` changes the key, so the next run saves a new cache; the old one ages out.
- The cibuildwheel key of the Linux and macOS wheel jobs still hashes the workflow file; it caches downloads only and is small, so it was left alone.
