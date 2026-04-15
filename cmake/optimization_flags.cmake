# ═══════════════════════════════════════════════════════════════════════════════
# Optimization Flags — architecture targeting, compiler flags, static linking
#
# Configures Release-mode compiler flags for maximum throughput on supported
# platforms (Linux x86_64, Linux ARM64, macOS ARM64).
# ═══════════════════════════════════════════════════════════════════════════════

# ── Static Build ──────────────────────────────────────────────────────────────
# Prefer .a over .so/.dylib for find_package(), and pass -static to linker.
# macOS static override is handled in platform_checks.cmake.
if (${LANCET_BUILD_STATIC})
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
endif ()

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
# -ffast-math
#   Relaxes IEEE 754 floating-point compliance for speed. Implies:
#     -ffinite-math-only   — assumes no NaN/Inf (makes std::isnan() UB)
#     -fno-signed-zeros    — ignores -0.0 vs +0.0 distinction
#     -fno-trapping-math   — assumes FP ops never trap
#     -fno-math-errno      — does not set errno on math functions
#     -fassociative-math   — allows FP reassociation (e.g., sum reordering)
#     -freciprocal-math    — allows x/y → x*(1/y) transformation
#   Lancet2 design: uses std::optional<f64> instead of NaN sentinels
#   precisely because -ffinite-math-only is enabled.
#
# -fno-omit-frame-pointer
#   Preserves frame pointers in Release builds. Costs 1 register but enables
#   stack-based profiling with pprof, perf, and Instruments. Required for
#   production performance analysis.
#
# Architecture flags (selected automatically per platform):
#   ARM64:  always -mcpu=native (ARM binaries are rarely cross-distributed)
#   x86_64: -march=x86-64-v3 (portable, AVX2 + BMI2, Haswell 2013+)
#           -march=native    (non-portable, when -DLANCET_NATIVE_BUILD=ON)

if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64" OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "aarch64")
	# ARM64: always native — ARM binaries target the build host's microarchitecture
	set(LANCET_OPT_FLAGS "-O3 -pipe -ffast-math -DNDEBUG -mcpu=native"
			CACHE STRING "Lancet optimization flags" FORCE)
elseif (LANCET_NATIVE_BUILD)
	# x86 native: uses all instruction sets available on this CPU (non-portable)
	message(STATUS "LANCET_NATIVE_BUILD=ON: using -march=native (non-portable binary)")
	set(LANCET_OPT_FLAGS "-O3 -pipe -ffast-math -DNDEBUG -march=native"
			CACHE STRING "Lancet optimization flags" FORCE)
else ()
	# x86 portable: AVX2 + BMI2 baseline (covers ~95% of modern server CPUs)
	set(LANCET_OPT_FLAGS "-O3 -pipe -ffast-math -DNDEBUG -march=x86-64-v3"
			CACHE STRING "Lancet optimization flags" FORCE)
endif ()

# ── Compiler-Specific Flags ───────────────────────────────────────────────────
# Clang: suppress -Wnan-infinity-disabled warnings caused by -ffinite-math-only
#        in third-party headers that reference NaN/Inf constants.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
	set(LANCET_OPT_FLAGS "${LANCET_OPT_FLAGS} -Wno-nan-infinity-disabled")
endif()

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
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} ${LANCET_OPT_FLAGS} -fno-omit-frame-pointer"
		CACHE INTERNAL "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${LANCET_OPT_FLAGS} -fno-omit-frame-pointer"
		CACHE INTERNAL "" FORCE)
