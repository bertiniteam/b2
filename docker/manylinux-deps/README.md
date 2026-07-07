# `b2-manylinux-deps` — prebuilt Boost + eigenpy for the wheel matrix

The Linux wheel jobs used to recompile **Boost.Python + eigenpy in every run** — the
bulk of the ~hour CI. `libboost_python3X` is ABI-locked per CPython, so it can't just be
built once and shared inside a single run. This image bakes it (and eigenpy) **once per
CPython**, so wheel builds just link against prebuilt libraries.

## Why an image, not `actions/cache`

The GitHub **Actions cache is 10 GB/repo with LRU eviction**; Boost+eigenpy prefixes
across (OS × Python) overrun it and thrash. A **GHCR image** lives in a *separate*
package-storage quota (unbounded for public images) and the runner layer-caches the pull,
so it sidesteps the ceiling entirely.

## The ABI contract (load-bearing)

Prebuilt `libboost_python` must match, exactly: the manylinux glibc/libstdc++ ABI, the
CPython ABI, the Boost version, and the Eigen/eigenpy it links. A mismatch is not a build
error — it's a runtime `SIGABRT`/`SIGSEGV` (this repo has that history, ADR-0006). So:

- The image **tag encodes the toolchain key**:
  `boost<BOOST>-eigen<EIGEN>-eigenpy<EIGENPY>-mlx2_34`.
- The versions are also stamped as image **labels** (`org.bertini.*`).
- **Main CI pins the exact key** (never `:latest`) and **asserts** the labels match its
  own `BOOST_VERSION`/`EIGENPY_VERSION` before building. Bump a version → new tag → new
  image; a stale image can never be silently consumed.

## The complete build+test env

The image is the whole Linux env, factored and assembled once — not just the C++ deps:

- `/opt/deps/<cpXY-cpXY>/` — Boost (incl. `libboost_python3X`) + eigenpy for that CPython.
- `/usr/local/` — Eigen headers (CPython-independent).
- **OpenMPI** (`/usr/lib64/openmpi`, on `PATH`) — so `pip install mpi4py` builds in the
  wheel *test* venv, which **un-skips the MPI tests** (`test_mpi_zerodim.py`,
  `test_doc_example_scripts.py`). Toward the 0-skiptests goal.
- **ccache** and a **pinned patchelf 0.17.2.1** — moved off `CIBW_BEFORE_ALL_LINUX`.

Main CI (`build_and_test.yml`) then sets, per active Python tag:
`CMAKE_PREFIX_PATH=/opt/deps/<tag>` and `LD_LIBRARY_PATH=/opt/deps/<tag>/lib:...`;
`CIBW_BEFORE_ALL_LINUX` collapses to ~nothing and `CIBW_BEFORE_BUILD_LINUX` drops from
"compile Boost+eigenpy" to a no-op. Python test deps stay in `CIBW_TEST_REQUIRES_LINUX`
(cibuildwheel's test phase is an isolated venv, so it re-pip-installs them — cheap for
pure-Python deps; `mpi4py` now *builds* because OpenMPI is present).

## Building / publishing

Built and pushed by `.github/workflows/build-ci-image.yml` → `ghcr.io/bertiniteam/b2-manylinux-deps`.
Triggers: a push to `develop` touching `docker/manylinux-deps/**`, or manual dispatch.
It **must run in `bertiniteam/b2`** (its token can only push to that org's GHCR).

**One-time setup:** after the first push, make the package **public** in the org's
Packages settings so any CI (incl. fork PRs) can pull it without auth. (Owner action —
GHCR package visibility is a package setting, not something CI can flip.)

## Follow-ups (not yet done)

- **Size:** Boost headers (~200 MB) are duplicated per Python for recipe fidelity. A
  later pass can install the CPython-independent headers + non-python libs once and keep
  only `libboost_python` + eigenpy per Python, shrinking the image several-fold.
- **macOS / Windows:** no container there. The same idea via **prebuilt tarballs as
  release assets** (also off the 10 GB cache ceiling) is the fast-follow.
- **Version single-sourcing:** the pins here mirror `build_and_test.yml`'s env. Keep them
  in lockstep (a mismatch just misses the assert and rebuilds — safe, but slow).
