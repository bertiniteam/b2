#!/usr/bin/env bash
#
# Build Boost (incl. Boost.Python) + eigenpy for ONE CPython tag into /opt/deps/<tag> and tar it.
#
# Run inside the manylinux container, ONE tag per invocation, so build-ci-image.yml can fan the
# per-CPython builds across a matrix (parallel) instead of a serial Dockerfile loop.  The assembly
# Dockerfile then just COPYs the resulting tarballs in.  Recipe mirrors the wheel job's historical
# CIBW_BEFORE_BUILD_LINUX exactly.
#
# Usage (from the repo root, cwd mounted at /work):
#   docker run --rm -v "$PWD:/work" -w /work -e BOOST_VERSION -e EIGEN_VERSION -e EIGENPY_VERSION \
#     quay.io/pypa/manylinux_2_34_x86_64 bash docker/manylinux-deps/build-python-deps.sh cp312-cp312
#
# Emits ./deps-<tag>.tar.gz (paths rooted at opt/deps/<tag>, so the image extracts it with -C /).

set -euxo pipefail
TAG="${1:?usage: build-python-deps.sh <cpXY-cpXY>}"
: "${BOOST_VERSION:?}"; : "${EIGEN_VERSION:?}"; : "${EIGENPY_VERSION:?}"
OUT="$PWD"   # mounted workspace, before we cd away

yum install -y wget gmp-devel mpfr-devel libmpc-devel libtool >/dev/null
yum clean all

# Eigen -> /usr/local, needed only to BUILD eigenpy here (the assembly image ships its own copy for
# the bertini compile; the tarball carries only /opt/deps/<tag>).
cd /tmp
wget -q "https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz"
tar xzf "eigen-${EIGEN_VERSION}.tar.gz"
cmake -S "eigen-${EIGEN_VERSION}" -B /tmp/eigen-bld -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_BUILD_TYPE=Release
cmake --install /tmp/eigen-bld

PY="/opt/python/${TAG}/bin/python"
[ -x "$PY" ] || { echo "missing interpreter: $PY"; exit 1; }
"$PY" -m pip install -q numpy scipy
PREFIX="/opt/deps/${TAG}"
PY_VER="$("$PY" -c 'import sys;print(f"{sys.version_info.major}.{sys.version_info.minor}")')"
PY_INC="$("$PY" -c 'import sysconfig;print(sysconfig.get_path("include"))')"
PY_LIB="$("$PY" -c 'import sysconfig;print(sysconfig.get_config_var("LIBDIR") or "")')"

cd /tmp
BU="boost_$(echo "${BOOST_VERSION}" | tr . _)"
wget -q "https://archives.boost.io/release/${BOOST_VERSION}/source/${BU}.tar.bz2"
tar xjf "${BU}.tar.bz2"
cd "${BU}"
./bootstrap.sh --with-python="$PY" --prefix="$PREFIX"
printf 'using python : %s : %s : %s : %s ;\n' "$PY_VER" "$PY" "$PY_INC" "$PY_LIB" > user-config.jam
./b2 install -j"$(nproc)" --user-config=user-config.jam python="$PY_VER" --without-mpi

cd /tmp
wget -q "https://github.com/stack-of-tasks/eigenpy/releases/download/v${EIGENPY_VERSION}/eigenpy-${EIGENPY_VERSION}.tar.gz"
tar xzf "eigenpy-${EIGENPY_VERSION}.tar.gz"
cmake -S "eigenpy-${EIGENPY_VERSION}" -B /tmp/eigenpy-bld \
  -DCMAKE_PREFIX_PATH="$PREFIX" -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON -DPython3_EXECUTABLE="$PY" \
  -DPython3_NumPy_INCLUDE_DIR="$("$PY" -c 'import numpy;print(numpy.get_include())')" \
  -DBUILD_TESTING=OFF -DCMAKE_INSTALL_DO_STRIP=ON
cmake --build /tmp/eigenpy-bld -j2 --target install

find "$PREFIX" -name '*.so*' -type f -exec strip --strip-unneeded {} + 2>/dev/null || true

tar czf "${OUT}/deps-${TAG}.tar.gz" -C / "opt/deps/${TAG}"
chmod a+r "${OUT}/deps-${TAG}.tar.gz"
ls -lh "${OUT}/deps-${TAG}.tar.gz"
