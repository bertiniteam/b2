# ADR-0078: The Windows CI jobs cache conda package archives, keyed by what was installed

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
- Without conda's lld, CMake took `lld-link`, `llvm-rc` and `llvm-mt` from the runner's own LLVM (`C:\Program Files\LLVM`), which nothing pins.  The Windows builds pass `CMAKE_LINKER`, `CMAKE_AR`, `CMAKE_RC_COMPILER` and `CMAKE_MT` explicitly: `link.exe` and `lib.exe` from the Visual Studio whose STL the build already uses (`VCToolsInstallDir`), `rc.exe` and `mt.exe` from its Windows SDK (`WindowsSdkVerBinPath`), as full paths, since Git's `link.exe` can come first on PATH.
- The first run without lld showed that b2 used Boost.JSON header-only but let Boost auto-link `boost_json.lib` by name on Windows; the old linker found one in conda's library folder.  `bertini2` now defines `BOOST_JSON_NO_LIB`.
- The cache holds the package archives only (`conda_pkgs/*.conda`, `conda_pkgs/*.tar.bz2`); conda unpacks an archive it finds.
- The key is what was installed: `windows-conda-archives-<python>-<SHA-256 of the sorted conda list --explicit --md5>`, with `cpp` for the C++ test job.  Each job restores the newest cache for its Python by prefix, installs (conda reuses what it finds and downloads what changed), records the installed packages, prunes the folder to their archives, and saves under the record's hash.  An unchanged environment finds its key taken and saves nothing; a conda-forge release, or an edit to `environment-win.yml`, makes a new key, and the superseded cache ages out after 7 days unused.  The wheel tests restore the cache of their Python and save nothing.
- One cache per Python version (~320 MB each, ~1.9 GB for the six) rather than one shared cache: the shared part (282 MB of 320) is stored six times, but a shared cache keyed by contents needs a job that solves all six before any Windows job starts, on the critical path of every run.
- The workflow file is not in the key, so editing the workflow keeps the cache.  A key that hit forever would be wrong the other way: a hit saves nothing, so every conda-forge release since the save would be downloaded on every run.

## Consequences

- Do not go back to caching the whole package folder: at ~2.6 GB per Python version it crowds the other platforms out of the 10 GB.
- Do not put the workflow file, the environment file or a date into the conda key; the record of what was installed already changes exactly when the contents do.
- Do not drop the pruning step: without it the restored archives of superseded versions are saved again, and the cache grows with every release.
- Do not reinstate `compilers` to get clang-cl; ask for `clang`.  If the STL raises its Clang floor, raise the pin.
- Do not drop the explicit MSVC tools for "CMake finds a linker anyway": it finds whatever LLVM the runner image ships.  Do not pass them as bare names: Git's `link.exe` is not a linker.
- Do not remove `BOOST_JSON_NO_LIB` while Boost.JSON is used header-only.
- In the wheel jobs the record comes after the build, because the build's `conda install python=... --update-deps` can change the environment.
- The cibuildwheel key of the Linux and macOS wheel jobs still hashes the workflow file; it caches downloads only and is small, so it was left alone.
