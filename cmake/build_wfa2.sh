#!/usr/bin/env bash

set -euxo pipefail

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"
readonly CXX_COMPILER="${3}"

readonly LANCET_OPT_FLAGS=$(grep 'LANCET_OPT_FLAGS' "${ROOT_DIR}/../../CMakeCache.txt" | cut -d= -f2-)

echo "WFA2 BUILD ROOT_DIR     : ${ROOT_DIR}"
echo "WFA2 BUILD C_COMPILER   : ${C_COMPILER}"
echo "WFA2 BUILD CXX_COMPILER : ${CXX_COMPILER}"
echo "WFA2 BUILD OPT_FLAGS    : ${LANCET_OPT_FLAGS}"

# WFA2's Makefile line 48 uses target-specific append: 'all: CC_FLAGS+=-O3 -march=native'.
# In GNU Make, target-specific += appends AFTER command-line values, so -march=native
# would override our -march=x86-64-v3 (GCC uses the last -march= flag). Strip it.
sed -i 's/-march=native//' "${ROOT_DIR}/Makefile"

# Create required output directories (Makefile's 'setup' target does this,
# but 'lib_wfa' doesn't depend on 'setup').
mkdir -p "${ROOT_DIR}/build" "${ROOT_DIR}/build/cpp" "${ROOT_DIR}/lib"

# Build only static libraries (libwfa.a + libwfacpp.a), skip tools/examples.
# libwfacpp.a is self-contained (includes all C objects + C++ bindings).
env CC="${C_COMPILER}" CXX="${CXX_COMPILER}" \
  make -C "${ROOT_DIR}" lib_wfa \
    BUILD_TOOLS=0 BUILD_EXAMPLES=0 \
    CC_FLAGS="${LANCET_OPT_FLAGS}"
