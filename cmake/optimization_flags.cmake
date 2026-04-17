# ═══════════════════════════════════════════════════════════════════════════════
# Optimization Flags — frame pointers, architecture targeting, static linking
#
# Configures compiler flags for all build modes plus Release-specific
# optimization flags for maximum throughput on supported platforms
# (Linux x86_64, Linux ARM64, macOS ARM64).
# ═══════════════════════════════════════════════════════════════════════════════

# ── Static Build ──────────────────────────────────────────────────────────────
# Prefer .a over .so/.dylib for find_package(), and pass -static to linker.
# macOS static override is handled in platform_checks.cmake.
if (${LANCET_BUILD_STATIC})
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
endif ()

# ── Frame Pointers (all build modes) ─────────────────────────────────────────
# Preserve frame pointers in every build type (Debug, Release, RelWithDebInfo).
# Costs 1 register (~1-2% on x86_64) but enables:
#   - gperftools' frame-pointer stack unwinder (generic_fp), the only crash-safe
#     unwinder for CPU profiling from a signal handler
#   - Reliable backtraces in the crash handler, perf, GDB, and Instruments
#
# Set on CMAKE_C/CXX_FLAGS (base flags) so it applies unconditionally.
# Also included in LANCET_OPT_FLAGS below to propagate to external dependencies
# (spoa, mimalloc, zlib-ng, htslib, etc.) via their configure scripts.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")

# ── Release Optimization Flags ────────────────────────────────────────────────
# Each flag and its relevance to Lancet2:
#
# -O3
#   Maximum optimization: loop unrolling, auto-vectorization, aggressive
#   inlining. Critical for alignment scoring and k-mer hashing hot paths.
#
# -pipe
#   Uses pipes instead of temporary files between compiler stages.
#   Speeds up compilation with no runtime effect.
#
# -DNDEBUG
#   Disables assert() and LANCET_ASSERT macros in Release builds. Required
#   because we override CMAKE_CXX_FLAGS_RELEASE entirely (CMake's default
#   Release flags include this, but our override replaces them).
#
# -fno-math-errno
#   Does not set errno after math library calls (pow, log, sqrt). Lancet2
#   never checks errno for math functions — this is pure overhead removal.
#
# -fno-trapping-math
#   Assumes floating-point operations never generate hardware traps.
#   Lancet2 does not install FP exception handlers.
#
# NOTE: -ffast-math was removed because its dangerous sub-flags caused
# intermittent segfaults via silent UB:
#   -ffinite-math-only  makes std::isnan()/std::isinf() return false,
#                       causing NaN to propagate silently as UB.
#   -funsafe-math-optimizations  allows reassociation that can introduce
#                       NaN in code that was otherwise safe.
# The safe sub-flags above give ~80% of the speedup with zero risk.
# Lancet2's FP-heavy code (polar_coords.h, scoring) is already hand-
# optimized to avoid expensive libm calls.
#
# -fno-omit-frame-pointer
#   Duplicated here from the base CMAKE_C/CXX_FLAGS above so that external
#   dependencies also get frame pointers when their configure scripts read
#   LANCET_OPT_FLAGS from CMakeCache.txt.  See the "Frame Pointers" section.
#
# Architecture flags (selected automatically per platform):
#   ARM64:  always -mcpu=native (ARM binaries are rarely cross-distributed)
#   x86_64: -march=x86-64-v3 (portable, AVX2 + BMI2, Haswell 2013+)
#           -march=native    (non-portable, when -DLANCET_NATIVE_BUILD=ON)

if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64" OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "aarch64")
	# ARM64: always native — ARM binaries target the build host's microarchitecture
	set(LANCET_OPT_FLAGS "-O3 -pipe -fno-math-errno -fno-trapping-math -fno-omit-frame-pointer -DNDEBUG -mcpu=native"
			CACHE STRING "Lancet optimization flags" FORCE)
elseif (LANCET_NATIVE_BUILD)
	# x86 native: uses all instruction sets available on this CPU (non-portable)
	message(STATUS "LANCET_NATIVE_BUILD=ON: using -march=native (non-portable binary)")
	set(LANCET_OPT_FLAGS "-O3 -pipe -fno-math-errno -fno-trapping-math -fno-omit-frame-pointer -DNDEBUG -march=native"
			CACHE STRING "Lancet optimization flags" FORCE)
else ()
	# x86 portable: AVX2 + BMI2 baseline (covers ~95% of modern server CPUs)
	set(LANCET_OPT_FLAGS "-O3 -pipe -fno-math-errno -fno-trapping-math -fno-omit-frame-pointer -DNDEBUG -march=x86-64-v3"
			CACHE STRING "Lancet optimization flags" FORCE)
endif ()

# Clang: no longer need -Wno-nan-infinity-disabled since -ffinite-math-only
# is no longer enabled (was caused by -ffast-math).

# GCC on Linux: -fno-semantic-interposition allows the compiler to assume that
# global symbols will not be interposed at runtime, enabling more aggressive
# devirtualization and inlining across translation units.
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
	set(LANCET_OPT_FLAGS "${LANCET_OPT_FLAGS} -fno-semantic-interposition")
endif()

# Linux: -fno-plt avoids the Procedure Linkage Table indirection for external
# function calls, reducing branch mispredictions. Not available on macOS
# (Mach-O uses a different dynamic linking mechanism).
if (CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
	set(LANCET_OPT_FLAGS "${LANCET_OPT_FLAGS} -fno-plt")
endif()

# ── Apply flags to Release configuration ─────────────────────────────────────
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} ${LANCET_OPT_FLAGS}"
		CACHE INTERNAL "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${LANCET_OPT_FLAGS}"
		CACHE INTERNAL "" FORCE)

# ── Apply flags to RelWithDebInfo configuration ──────────────────────────────
# Same optimization flags as Release, plus -g for DWARF debug symbols.
# Used exclusively with LANCET_PROFILE_MODE=ON for source-level profiling
# via pprof --list. The -g flag adds ~2× binary size but zero runtime cost.
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS} ${LANCET_OPT_FLAGS} -g"
		CACHE INTERNAL "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${LANCET_OPT_FLAGS} -g"
		CACHE INTERNAL "" FORCE)
