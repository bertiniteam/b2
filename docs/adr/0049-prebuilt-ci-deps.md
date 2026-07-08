# ADR-0049: Prebuilt CI dependencies — image (Linux), tarballs (macOS), conda (Windows)

**Status:** Accepted
**Date:** 2026-07-07

## Context

Building the wheels recompiled **Boost (incl. Boost.Python) + eigenpy from source in every
CI run, on every (platform × CPython) cell** — the dominant cost of the ~hour-long matrix.
This is not gratuitous: `libboost_python3X` is **ABI-locked to a single CPython version**, so
one build cannot be shared across interpreters within a run; each of ~15 cells rebuilt it.

The obvious fix — cache the built prefixes with `actions/cache` — was tried and **fails**:
the GitHub Actions cache is **10 GB per repo with LRU eviction**, and Boost+eigenpy prefixes
across the matrix overrun it and thrash (evict → rebuild → re-cache → evict). So caching is
the wrong tool; the deps must live in storage *outside* that ceiling.

A second, load-bearing constraint: prebuilt `libboost_python` must match, **exactly**, the
runtime ABI — glibc/libstdc++ (manylinux), CPython ABI, Boost version, and the Eigen/eigenpy
it links. A mismatch is not a build error but a runtime `SIGABRT`/`SIGSEGV` (this repo has
that history — see ADR-0006). Any prebuild scheme must make a mismatch **impossible to
consume silently**.

## Decision

Prebuild the dependencies **once per toolchain** and have CI *consume* them, per platform:

- **Linux → a GHCR image.** `ghcr.io/bertiniteam/b2-manylinux-deps`, built by
  `.github/workflows/build-ci-image.yml`: a **matrix** builds Boost.Python + eigenpy for each
  CPython in parallel (inside the manylinux container via `docker run`,
  `docker/manylinux-deps/build-python-deps.sh`), and a thin **assemble** job unpacks them into
  the image alongside the system layer (eigen/ccache/patchelf/**OpenMPI**). The wheel job sets
  `CIBW_MANYLINUX_X86_64_IMAGE` to the pinned tag; `CIBW_BEFORE_BUILD_LINUX` just symlinks the
  active Python's prefix (`/opt/deps/<cpXY-cpXY>` → `/opt/deps/current`). Image/package storage
  is a separate, effectively-unbounded quota (sidesteps the 10 GB ceiling); the runner
  layer-caches the pull.

- **macOS → GitHub Release-asset tarballs.** No container, so the analog is prebuilt
  `deps-macos14-arm64-<cpXY>-<key>.tar.gz` on the `ci-deps` prerelease, built by
  `build-macos-deps.yml`. `CIBW_BEFORE_BUILD_MACOS` downloads the matching bundle (with a
  **build-from-source fallback** if it's absent), extracting to `/tmp/deps-py` so its
  `hardcode-dll-paths` stay valid for delocate. Also off the 10 GB ceiling (release storage).

- **Windows → nothing new.** conda-forge already ships Boost/boost-python/eigenpy as binaries,
  and `environment-win.yml` already includes `msmpi`+`mpi4py`. The only from-source cost is
  bertini itself (a ccache experiment addresses recompiles).

**The ABI contract (the crux).** The image *tag* and the tarball *name* are the **toolchain
key**: `boost<B>-eigen<E>-eigenpy<EP>-mlx2_34`. The consumer pins that exact key, derived from
the versions. If a version is bumped but the prebuilt artifact for the new key doesn't exist
yet, the **pull 404s and the job fails loudly** — the implicit ABI assert. There is no path by
which a stale, wrong-ABI binary is silently consumed.

**Single source of versions.** `.github/ci-deps-versions.env` holds `BOOST_VERSION` /
`EIGEN_VERSION` / `EIGENPY_VERSION`, read by all three workflows. It is in the **path filters**
of the image and macOS producers, so a bump there both reruns the wheel CI and **retriggers the
prebuild** — no more editing versions in four places and forgetting to rebuild the deps.

**0 skiptests.** The image (and macOS via `brew open-mpi`) provide OpenMPI so `pip install
mpi4py` succeeds in the test venv; the MPI test modules then run at 1 rank instead of skipping,
matching Windows (which already ships mpi4py).

## Consequences

- Wheel CI stops recompiling Boost per run; the prebuild runs **rarely** — only when the
  Dockerfile/recipe or a pinned version changes (or on manual dispatch).
- The Linux build+test env is a single, transparent, versioned artifact — the natural
  foundation for a **conda-forge feedstock (with MPI variants)** and a **Homebrew tap** later.
- The image must be built **in `bertiniteam/b2`** (its token pushes to that org's GHCR) and the
  package made **public** once (so fork-PR CI can pull it).
- Gotcha fixed here: `test_wheels_linux_macos` *excludes* ubuntu, so an ubuntu-only os matrix
  (a non-develop PR base, or `-f os=linux`/`os=windows`) expanded to zero combinations and
  failed the run with nothing to point at; it is now guarded with
  `contains(os, 'macos-14')`.
- The image is intentionally larger than necessary (Boost headers duplicated per CPython for
  recipe fidelity); a later pass can share the CPython-independent headers to shrink it.

See `docker/manylinux-deps/README.md` for the operational detail.
