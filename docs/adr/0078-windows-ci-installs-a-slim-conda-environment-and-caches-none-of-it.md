# ADR-0078: The Windows CI jobs install a slim conda environment and cache none of it

**Status:** Accepted
**Date:** 2026-11-06

## Context

GitHub keeps 10 GB of Actions caches per repository and evicts the least recently used beyond that.  On 2026-11-04 the repository held 8.9 GB in 23 caches, all on develop, and no Linux cache at all.  Each Windows job cached its conda package folder under its own key (`windows-conda-cpp-...` and `windows-conda-<python>-...`), about 2.6 GB apiece: six of them do not fit, and while they took turns evicting each other, they evicted every other platform's caches and every pull request's.  A PR restores only its own caches and develop's, so PRs compiled cold.

Three things made the folders large and many:

- **What the environment installs.**  `environment-win.yml` listed the `compilers` metapackage, whose Fortran compiler (flang) brings llvm-tools, lld and clang: 830 MB of a 1153 MB download.  The builds use only clang-cl from it.  conda's git (115 MB) duplicated the runner's, and numpy and scipy took MKL (109 MB) as their BLAS.
- **Unpacked copies.**  A conda package folder holds each package twice, the `.conda` archive and its unpacked folder, so 1.15 GB of downloads made a 2.6 GB cache.
- **One cache per Python version, keyed by the workflow file.**  The solves for the Python versions share 282 MB of ~320 MB.  And `hashFiles(..., '.github/workflows/build_and_test.yml')` in the key made every workflow edit a new set of caches.

## What was measured

With the environment slimmed (below), PR #527 cached the package archives only, ~325 MB per job, keyed by the SHA-256 of `conda list --explicit --md5`, and logged for every job how many archives it restored, how many conda reused, and how many it downloaded.  Over three runs (2026-11-05):

| | jobs | conda install (setup-miniconda) |
|---|---|---|
| no cache to restore | 16 | 116-252 s; 130-179 s in the last run |
| cache restored (81 archives, 326.7 MB, in 5 s) | 2 | 211 s and 178 s |

In the job that restored the cache, conda installed the same 81 packages, and all 81 archives were rewritten after the restore: conda used none of them and fetched every package again.  The install time varies by more than two minutes from runner to runner, so runner speed, not the cache, decides it.  And the most a cache that conda did reuse could save is the download of ~325 MB from conda-forge's CDN, which costs about what restoring a cache of the same size does (5 s here).

## Decision

- **No conda cache on Windows.**  The three Windows jobs (C++ tests, wheels, wheel tests) install the environment from conda-forge on every run; `CONDA_PKGS_DIRS` is not set.  The ~1.9 GB the per-Python caches would hold stays free for the caches that pay: ccache on every platform, and the Linux and macOS caches.
- **The environment is slim.**  `environment-win.yml` asks for `clang >=20` (clang-cl) instead of `compilers`, drops git, and pins numpy's BLAS to OpenBLAS (`libblas=*=*openblas`).  The download is ~325 MB instead of 1153 MB.  The Clang floor is the one Visual Studio 2026's STL demands (STL1000), so the builds no longer define `_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH`.
- **The link tools are Visual Studio's.**  Without conda's lld, CMake took `lld-link`, `llvm-rc` and `llvm-mt` from the runner's own LLVM (`C:\Program Files\LLVM`), which nothing pins.  `tools/msvc_link_tools.ps1` enters the developer environment of the newest Visual Studio (found with `vswhere`; the one clang-cl finds by itself, whose STL it compiles against) for `LIB`, `LIBPATH` and `INCLUDE`, which `link.exe` needs and the conda activation does not set in PowerShell, puts PATH back so the build keeps conda's cmake, ninja, ccache and clang-cl, and returns `link.exe` and `lib.exe` from its MSVC toolset and `rc.exe` and `mt.exe` from its Windows SDK.  The builds pass them as `CMAKE_LINKER`, `CMAKE_AR`, `CMAKE_RC_COMPILER` and `CMAKE_MT`, by full path, since Git's `link.exe` can come first on PATH.
- **Boost.JSON does not auto-link.**  The first run without lld failed with `could not open 'boost_json.lib'`: b2 uses Boost.JSON header-only, but on Windows Boost's auto-linking asked for the library by name, and only the old linker happened to find one in conda's library folder.  `bertini2` defines `BOOST_JSON_NO_LIB`.
- **The Windows build steps log what they used and fail where they fail**: the main conda package versions, which cmake, ninja, clang-cl and ccache run, `clang-cl --version`, and the compiler and link tools in `CMakeCache.txt`.  A native command's exit code does not stop a PowerShell step, so the first failed link surfaced only as ten tests "Not Run"; the build steps now exit with the build's code.

## Consequences

- Do not add a conda cache back to the Windows jobs without measuring that it saves time: log, as PR #527 did, how many restored archives conda reuses and how long the install takes with and without it.  A cache that only mirrors a fast download costs cache space and saves nothing.
- If a cache is ever measured to pay, do not cache the whole package folder (~2.6 GB per Python version crowds the other platforms out of the 10 GB) or key it on the workflow file (every edit makes new caches).
- Do not reinstate `compilers` to get clang-cl; ask for `clang`.  If the STL raises its Clang floor, raise the pin.
- Do not drop the explicit MSVC tools for "CMake finds a linker anyway": it finds whatever LLVM the runner image ships.  Do not pass them as bare names: Git's `link.exe` is not a linker.  Do not leave Visual Studio's PATH in place after entering its developer shell: its CMake and Ninja would replace conda's.
- Do not remove `BOOST_JSON_NO_LIB` while Boost.JSON is used header-only.
- The cibuildwheel key of the Linux and macOS wheel jobs still hashes the workflow file; it caches downloads only and is small, and it was not measured here.
