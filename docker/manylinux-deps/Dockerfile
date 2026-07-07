# syntax=docker/dockerfile:1
#
# Assembly image: the manylinux system layer + prebuilt per-CPython Boost/eigenpy bundles COPYed in.
#
# The slow per-CPython Boost+eigenpy builds run in PARALLEL matrix jobs (build-ci-image.yml, via
# build-python-deps.sh) and are handed here as tarballs.  This Dockerfile therefore only installs
# the CPython-independent system layer and unpacks those tarballs, so it builds in ~a couple of
# minutes instead of the ~hour a serial per-Python loop took.
#
# NB: not standalone-buildable -- it expects deps/*.tar.gz in the build context, which the CI
# assemble job populates from the matrix artifacts.
#
# ABI CONTRACT: the image *tag* encodes the toolchain (boost/eigen/eigenpy + manylinux), stamped
# also as labels; main CI pins the exact tag and asserts the labels match.  See README.md, ADR-0006.

ARG MANYLINUX_IMAGE=quay.io/pypa/manylinux_2_34_x86_64
FROM ${MANYLINUX_IMAGE}

ARG BOOST_VERSION=1.90.0
ARG EIGEN_VERSION=3.4.0
ARG EIGENPY_VERSION=3.13.0

LABEL org.opencontainers.image.source="https://github.com/bertiniteam/b2"
LABEL org.opencontainers.image.description="manylinux_2_34 + prebuilt Boost/Boost.Python + eigenpy for bertini2 wheels"
LABEL org.bertini.boost_version="${BOOST_VERSION}"
LABEL org.bertini.eigen_version="${EIGEN_VERSION}"
LABEL org.bertini.eigenpy_version="${EIGENPY_VERSION}"

# System layer (CPython-independent): numeric -devel libs, ccache, OpenMPI (so `pip install mpi4py`
# builds in the wheel test venv -> un-skips the MPI tests), and a pinned patchelf (the stock one
# corrupts the bundled bertini2 CLI's ELF and segfaults it at load).
RUN yum install -y epel-release \
 && yum install -y wget gmp-devel mpfr-devel libmpc-devel libtool ccache openmpi openmpi-devel \
 && yum clean all \
 && echo /usr/lib64/openmpi/lib > /etc/ld.so.conf.d/openmpi.conf && ldconfig
ENV PATH=/usr/lib64/openmpi/bin:${PATH}
RUN python3 -m venv /tmp/pe \
 && /tmp/pe/bin/pip install -q "patchelf==0.17.2.1" \
 && install -m755 /tmp/pe/bin/patchelf /usr/local/bin/patchelf \
 && rm -rf /tmp/pe \
 && patchelf --version

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
