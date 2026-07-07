# syntax=docker/dockerfile:1
#
# Prebuilt dependencies for bertini2 Linux wheels: Boost (incl. a per-CPython
# Boost.Python) + eigenpy, on top of the manylinux image cibuildwheel already uses.
#
# WHY: the wheel matrix recompiles Boost.Python + eigenpy in every run (the ~hour
# sink), because libboost_python3X is ABI-locked per CPython.  Baking them into a
# GHCR image moves that cost out of every run and off the 10 GB Actions-cache ceiling
# (image storage is a separate, effectively-unbounded quota; the runner layer-caches
# the pull).  The build recipe below MIRRORS build_and_test.yml's
# CIBW_BEFORE_ALL/BEFORE_BUILD_LINUX exactly -- just relocated here.
#
# ABI CONTRACT (load-bearing): the image *tag* encodes the toolchain
# (boost/eigen/eigenpy + manylinux), and main CI ASSERTS a match before using it.  A
# version bump = a new tag = a new image; a stale image can never be silently consumed.
# This repo has a SIGABRT history from ABI drift (ADR-0006), so this is not optional.
#
# Per-CPython prefixes land at /opt/deps/<tag> (e.g. /opt/deps/cp310-cp310).  Main CI
# points CMAKE_PREFIX_PATH / LD_LIBRARY_PATH at the active Python's prefix.
#
# SIZE NOTE: Boost headers (~200 MB) are duplicated across the per-Python prefixes for
# correctness/simplicity (this mirrors the known-good recipe).  A follow-up can share
# the CPython-independent headers + non-python libs once and keep only libboost_python
# + eigenpy per Python -- see README.md.

ARG MANYLINUX_IMAGE=quay.io/pypa/manylinux_2_34_x86_64
FROM ${MANYLINUX_IMAGE}

# Pins -- MUST stay in lockstep with BOOST_VERSION / EIGENPY_VERSION in build_and_test.yml.
ARG BOOST_VERSION=1.90.0
ARG EIGEN_VERSION=3.4.0
ARG EIGENPY_VERSION=3.13.0
# CPython tags to build for (manylinux ships these under /opt/python/<tag>).
ARG PYTHON_TAGS="cp310-cp310 cp311-cp311 cp312-cp312 cp313-cp313 cp314-cp314"

LABEL org.opencontainers.image.source="https://github.com/bertiniteam/b2"
LABEL org.opencontainers.image.description="manylinux_2_34 + prebuilt Boost/Boost.Python + eigenpy for bertini2 wheels"
# Toolchain stamp -- main CI reads these labels to assert an ABI match.
LABEL org.bertini.boost_version="${BOOST_VERSION}"
LABEL org.bertini.eigen_version="${EIGEN_VERSION}"
LABEL org.bertini.eigenpy_version="${EIGENPY_VERSION}"

# System layer (CPython-independent): numeric -devel libs, ccache, and OpenMPI so
# `pip install mpi4py` can build in the wheel *test* venv -- which un-skips the MPI tests
# (test_mpi_zerodim.py, test_doc_example_scripts.py).  This is the "complete build+test
# env" idea: everything the Linux wheel run needs, factored here and assembled once.
RUN yum install -y epel-release \
 && yum install -y wget gmp-devel mpfr-devel libmpc-devel libtool ccache openmpi openmpi-devel \
 && yum clean all \
 && echo /usr/lib64/openmpi/lib > /etc/ld.so.conf.d/openmpi.conf && ldconfig

# OpenMPI's mpicc/mpirun are not on the default PATH on AlmaLinux; expose them so mpi4py
# builds and (for the dedicated MPI job / serial runs) mpirun is available.
ENV PATH=/usr/lib64/openmpi/bin:${PATH}

# Pin patchelf 0.17.2.1 (mirrors CIBW_BEFORE_ALL_LINUX): the image's stock patchelf
# mis-rewrites the bundled bertini2 CLI's ELF and segfaults it at load.
RUN python3 -m venv /tmp/pe \
 && /tmp/pe/bin/pip install -q "patchelf==0.17.2.1" \
 && install -m755 /tmp/pe/bin/patchelf /usr/local/bin/patchelf \
 && rm -rf /tmp/pe \
 && patchelf --version

# Eigen (header-only, CPython-independent) -> /usr/local, same as CIBW_BEFORE_ALL_LINUX.
RUN cd /tmp \
 && wget -q "https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz" \
 && tar xzf "eigen-${EIGEN_VERSION}.tar.gz" \
 && cmake -S "eigen-${EIGEN_VERSION}" -B /tmp/eigen-bld -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release \
 && cmake --install /tmp/eigen-bld \
 && rm -rf /tmp/eigen*

# Fetch Boost + eigenpy sources once.
RUN cd /tmp \
 && BU="boost_$(echo "${BOOST_VERSION}" | tr . _)" \
 && wget -q "https://archives.boost.io/release/${BOOST_VERSION}/source/${BU}.tar.bz2" \
 && wget -q "https://github.com/stack-of-tasks/eigenpy/releases/download/v${EIGENPY_VERSION}/eigenpy-${EIGENPY_VERSION}.tar.gz"

# Build Boost.Python + eigenpy per CPython into /opt/deps/<tag>.  A fresh Boost tree
# per Python avoids cross-variant b2 contamination.  Mirrors CIBW_BEFORE_BUILD_LINUX.
RUN set -eux; cd /tmp; \
    BU="boost_$(echo "${BOOST_VERSION}" | tr . _)"; \
    for tag in ${PYTHON_TAGS}; do \
      PY="/opt/python/${tag}/bin/python"; \
      [ -x "$PY" ] || { echo "missing interpreter: $PY"; exit 1; }; \
      PREFIX="/opt/deps/${tag}"; \
      "$PY" -m pip install -q numpy scipy; \
      PY_VER="$("$PY" -c 'import sys;print(f"{sys.version_info.major}.{sys.version_info.minor}")')"; \
      PY_INC="$("$PY" -c 'import sysconfig;print(sysconfig.get_path("include"))')"; \
      PY_LIB="$("$PY" -c 'import sysconfig;print(sysconfig.get_config_var("LIBDIR") or "")')"; \
      rm -rf "/tmp/${BU}"; tar xjf "/tmp/${BU}.tar.bz2" -C /tmp; \
      cd "/tmp/${BU}"; \
      ./bootstrap.sh --with-python="$PY" --prefix="$PREFIX"; \
      printf 'using python : %s : %s : %s : %s ;\n' "$PY_VER" "$PY" "$PY_INC" "$PY_LIB" > user-config.jam; \
      ./b2 install -j"$(nproc)" --user-config=user-config.jam python="$PY_VER" --without-mpi; \
      cd /tmp; \
      rm -rf "eigenpy-${EIGENPY_VERSION}"; tar xzf "eigenpy-${EIGENPY_VERSION}.tar.gz"; \
      cmake -S "eigenpy-${EIGENPY_VERSION}" -B "/tmp/eigenpy-bld-${tag}" \
        -DCMAKE_PREFIX_PATH="$PREFIX" -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON -DPython3_EXECUTABLE="$PY" \
        -DPython3_NumPy_INCLUDE_DIR="$("$PY" -c 'import numpy;print(numpy.get_include())')" \
        -DBUILD_TESTING=OFF -DCMAKE_INSTALL_DO_STRIP=ON; \
      cmake --build "/tmp/eigenpy-bld-${tag}" -j2 --target install; \
      find "$PREFIX" -name '*.so*' -type f -exec strip --strip-unneeded {} + 2>/dev/null || true; \
      rm -rf "/tmp/eigenpy-bld-${tag}"; \
    done; \
    rm -rf "/tmp/${BU}" "/tmp/${BU}.tar.bz2" /tmp/eigenpy-${EIGENPY_VERSION}*
