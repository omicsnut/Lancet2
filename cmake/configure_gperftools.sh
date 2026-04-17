#!/usr/bin/env bash

set -euxo pipefail

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"
readonly CXX_COMPILER="${3}"

readonly LANCET_OPT_FLAGS=$(grep 'LANCET_OPT_FLAGS' "${ROOT_DIR}/../../CMakeCache.txt" | cut -d= -f2-)
readonly CFLAGS="${LANCET_OPT_FLAGS}"

echo "GPERFTOOLS CONFIGURE ROOT_DIR : ${ROOT_DIR}"
echo "GPERFTOOLS CONFIGURE C_COMPILER : ${C_COMPILER}"
echo "GPERFTOOLS CONFIGURE CXX_COMPILER : ${CXX_COMPILER}"

# --enable-frame-pointers: force gperftools to use the frame-pointer-based stack
# unwinder (stacktrace_generic_fp-inl.h) instead of libgcc's _Unwind_Backtrace.
#
# The libgcc unwinder reads DWARF .eh_frame tables and takes internal locks inside
# the SIGPROF signal handler, which crashes with SIGSEGV in x86_64_fallback_frame_state
# when the profiled thread is in a frame without DWARF info (e.g., statically-linked
# third-party libraries).  The frame-pointer unwinder walks the RBP chain directly —
# no DWARF parsing, no locks, fully async-signal-safe.
#
# This flag sets -DFORCED_FRAME_POINTERS at compile time, which places generic_fp
# first in gperftools' all_impls[] array.  This is the ONLY way to select it — an
# environment variable (TCMALLOC_STACKTRACE_METHOD) doesn't work because gperftools
# reads it during static initialization, before main() runs.
#
# See gperftools_crash_analysis.md for the full root-cause analysis.
env CC="${C_COMPILER}" CXX="${CXX_COMPILER}" CFLAGS="${CFLAGS}" "${ROOT_DIR}"/configure --prefix="${ROOT_DIR}" \
  --enable-static --disable-shared --disable-heap-profiler --disable-heap-checker --disable-debugalloc --enable-frame-pointers
