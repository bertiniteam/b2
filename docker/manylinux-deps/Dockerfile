# syntax=docker/dockerfile:1
#
# Assembly image: the manylinux system layer + prebuilt per-CPython Boost/eigenpy bundles COPYed in.
#
# The slow per-CPython Boost+eigenpy builds run in PARALLEL matrix jobs (build-ci-image.yml, via
# build-python-deps.sh) and are handed here as tarballs.  This Dockerfile therefore only installs
# the CPython-independent system layer and unpacks those tarballs, so it builds in minutes instead
# of the ~hour a serial per-Python loop took.
#
# NB: not standalone-buildable -- it expects deps/*.tar.gz in the build context, which the CI
# assemble job populates from the matrix artifacts.
#
# ABI CONTRACT: the image *tag* encodes the toolchain (boost/eigen/eigenpy/flint + manylinux),
# stamped also as labels; main CI pins the exact tag.  See README.md, ADR-0006, ADR-0049.

ARG MANYLINUX_IMAGE=quay.io/pypa/manylinux_2_34_x86_64
FROM ${MANYLINUX_IMAGE}

ARG BOOST_VERSION=1.90.0
ARG EIGEN_VERSION=3.4.0
ARG EIGENPY_VERSION=3.13.0
ARG GMP_VERSION=6.3.0
ARG MPFR_VERSION=4.2.2
ARG MPC_VERSION=1.3.1
ARG FLINT_VERSION=3.5.0

LABEL org.opencontainers.image.source="https://github.com/bertiniteam/b2"
LABEL org.opencontainers.image.description="manylinux_2_34 + prebuilt Boost/Boost.Python + eigenpy + GMP/MPFR/MPC/FLINT for bertini2 wheels"
LABEL org.bertini.boost_version="${BOOST_VERSION}"
LABEL org.bertini.eigen_version="${EIGEN_VERSION}"
LABEL org.bertini.eigenpy_version="${EIGENPY_VERSION}"
LABEL org.bertini.gmp_version="${GMP_VERSION}"
LABEL org.bertini.mpfr_version="${MPFR_VERSION}"
LABEL org.bertini.mpc_version="${MPC_VERSION}"
LABEL org.bertini.flint_version="${FLINT_VERSION}"

# System layer (CPython-independent): ccache, OpenMPI (so `pip install mpi4py` builds in the wheel
# test venv -> un-skips the MPI tests), and a pinned patchelf (the stock one corrupts the bundled
# bertini2 CLI's ELF and segfaults it at load).  /usr/local/lib joins the loader path so the
# number libraries built below are found by the build, by ldd, and so by auditwheel.
RUN yum install -y epel-release \
 && yum install -y wget xz libtool ccache openmpi openmpi-devel \
 && yum clean all \
 && echo /usr/lib64/openmpi/lib > /etc/ld.so.conf.d/openmpi.conf \
 && echo /usr/local/lib > /etc/ld.so.conf.d/local.conf \
 && ldconfig
ENV PATH=/usr/lib64/openmpi/bin:${PATH}
RUN python3 -m venv /tmp/pe \
 && /tmp/pe/bin/pip install -q "patchelf==0.17.2.1" \
 && install -m755 /tmp/pe/bin/patchelf /usr/local/bin/patchelf \
 && rm -rf /tmp/pe \
 && patchelf --version

# The number libraries, from source into /usr/local, one consistent set.  FLINT 3 (which carries
# Arb, the ball arithmetic behind b2's interval types) needs GMP >= 6.2.1, and the AlmaLinux 9 base
# ships 6.2.0, so GMP is built here -- fat on x86_64, so the CPU-specific code is chosen at run time
# rather than fixed to the build host -- and MPFR, MPC and FLINT are built against it.  Everything
# is stripped: an unstripped libflint.so adds about 25 MB to every wheel.
RUN cd /tmp \
 && wget -q "https://ftp.gnu.org/gnu/gmp/gmp-${GMP_VERSION}.tar.xz" \
 && tar xJf "gmp-${GMP_VERSION}.tar.xz" \
 && cd "gmp-${GMP_VERSION}" \
 && FAT="" && { [ "$(uname -m)" != x86_64 ] || FAT=--enable-fat; } \
 && ./configure --prefix=/usr/local --disable-static $FAT \
 && make -j"$(nproc)" \
 && make install-strip \
 && ldconfig \
 && rm -rf /tmp/gmp-*
RUN cd /tmp \
 && wget -q "https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.xz" \
 && tar xJf "mpfr-${MPFR_VERSION}.tar.xz" \
 && cd "mpfr-${MPFR_VERSION}" \
 && ./configure --prefix=/usr/local --disable-static --with-gmp=/usr/local \
 && make -j"$(nproc)" \
 && make install-strip \
 && ldconfig \
 && rm -rf /tmp/mpfr-*
RUN cd /tmp \
 && wget -q "https://ftp.gnu.org/gnu/mpc/mpc-${MPC_VERSION}.tar.gz" \
 && tar xzf "mpc-${MPC_VERSION}.tar.gz" \
 && cd "mpc-${MPC_VERSION}" \
 && ./configure --prefix=/usr/local --disable-static --with-gmp=/usr/local --with-mpfr=/usr/local \
 && make -j"$(nproc)" \
 && make install-strip \
 && ldconfig \
 && rm -rf /tmp/mpc-*
RUN cd /tmp \
 && wget -q "https://github.com/flintlib/flint/releases/download/v${FLINT_VERSION}/flint-${FLINT_VERSION}.tar.gz" \
 && tar xzf "flint-${FLINT_VERSION}.tar.gz" \
 && cd "flint-${FLINT_VERSION}" \
 && ./configure --prefix=/usr/local --disable-static --with-gmp=/usr/local --with-mpfr=/usr/local \
 && make -j"$(nproc)" \
 && make install \
 && strip --strip-unneeded /usr/local/lib/libflint.so.* \
 && ldconfig \
 && rm -rf /tmp/flint-* \
 && ls -l /usr/local/lib/libgmp.so* /usr/local/lib/libmpfr.so* /usr/local/lib/libmpc.so* /usr/local/lib/libflint.so*

# Eigen headers (CPython-independent) -> /usr/local, for the bertini compile in main CI.
RUN cd /tmp \
 && wget -q "https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz" \
 && tar xzf "eigen-${EIGEN_VERSION}.tar.gz" \
 && cmake -S "eigen-${EIGEN_VERSION}" -B /tmp/eigen-bld -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release \
 && cmake --install /tmp/eigen-bld \
 && rm -rf /tmp/eigen*

# Prebuilt per-CPython Boost + eigenpy bundles from the parallel build-deps matrix jobs.  Each
# tarball is rooted at opt/deps/<tag>, so `tar -C /` recreates /opt/deps/<tag>.
COPY deps/ /tmp/deps-tars/
RUN set -eux; for f in /tmp/deps-tars/*.tar.gz; do echo "unpacking $f"; tar xzf "$f" -C /; done \
 && rm -rf /tmp/deps-tars \
 && ls -1 /opt/deps
