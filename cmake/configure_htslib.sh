#!/usr/bin/env bash

set -euxo pipefail

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"
readonly ENABLE_CLOUD_IO="${3:-OFF}"
readonly INSTALL_PREFIX="${4:-}"

readonly ZLIBNG_BUILD_DIR=$(realpath "${ROOT_DIR}/../zlib-ng-build")
readonly LIBDEFLATE_INC_DIR=$(realpath "${ROOT_DIR}/../libdeflate-src")
readonly LIBDEFLATE_LIB_DIR=$(realpath "${ROOT_DIR}/../libdeflate-build")

readonly LANCET_OPT_FLAGS=$(grep 'LANCET_OPT_FLAGS' "${ROOT_DIR}/../../CMakeCache.txt" | cut -d= -f2-)

echo "HTSLIB CONFIGURE ROOT_DIR : ${ROOT_DIR}"
echo "HTSLIB CONFIGURE C_COMPILER : ${C_COMPILER}"
echo "HTSLIB CONFIGURE ENABLE_CLOUD_IO : ${ENABLE_CLOUD_IO}"
echo "HTSLIB CONFIGURE INSTALL_PREFIX : ${INSTALL_PREFIX}"

CLOUD_FLAGS="--disable-libcurl"
if [ "${ENABLE_CLOUD_IO}" = "ON" ]; then
    CLOUD_FLAGS="--enable-libcurl --enable-s3 --enable-gcs"
fi

# Build CPPFLAGS/LDFLAGS from vendored build dirs + inherited env.
# When CMAKE_INSTALL_PREFIX points to a conda host prefix (rattler-build),
# also add its include/lib so htslib's autotools can find system deps
# (liblzma, libcurl, openssl) that live there instead of /usr/local.
PREFIX_CPPFLAGS="${CPPFLAGS:-} -I${ZLIBNG_BUILD_DIR} -I${LIBDEFLATE_INC_DIR}"
PREFIX_LDFLAGS="${LDFLAGS:-} -L${ZLIBNG_BUILD_DIR} -L${LIBDEFLATE_LIB_DIR}"
if [ -n "${INSTALL_PREFIX}" ] && [ -d "${INSTALL_PREFIX}/include" ]; then
    PREFIX_CPPFLAGS="${PREFIX_CPPFLAGS} -I${INSTALL_PREFIX}/include"
    PREFIX_LDFLAGS="${PREFIX_LDFLAGS} -L${INSTALL_PREFIX}/lib"
fi

env CC="${C_COMPILER}" CFLAGS="${LANCET_OPT_FLAGS:-}" \
  "${ROOT_DIR}"/configure --prefix="${ROOT_DIR}" ${CLOUD_FLAGS} --disable-plugins --with-libdeflate \
  CPPFLAGS="${PREFIX_CPPFLAGS}" \
  LDFLAGS="${PREFIX_LDFLAGS}"
