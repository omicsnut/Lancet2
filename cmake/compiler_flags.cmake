# ═══════════════════════════════════════════════════════════════════════════════
# Compiler & Linker Flags — toolchain selection, base flags, per-configuration
#                            optimization, architecture targeting
#
# Configures the complete compiler and linker flag surface for all build
# configurations (Debug, Release, RelWithDebInfo) across supported platforms
# (Linux x86_64, Linux ARM64, macOS ARM64):
#
#   Static Linking      -static, .a library preference
#   Clang Runtime       libc++, compiler-rt, lld for a self-consistent toolchain
#   Frame Pointers      -fno-omit-frame-pointer for profiling and crash traces
#   Release Opts        -O3, safe math flags, platform-specific codegen tweaks
#   Architecture        -march=x86-64-v3 / -march=native / -mcpu=native
#   Flag Application    CMAKE_*_FLAGS_RELEASE, _RELWITHDEBINFO, _DEBUG assembly
#   Link-Time Opt       IPO/LTO for Release builds (if compiler supports it)
# ═══════════════════════════════════════════════════════════════════════════════

# ── Static Build ──────────────────────────────────────────────────────────────
# Prefer .a over .so/.dylib for find_package(), and pass -static to linker.
# macOS static override is handled in platform_checks.cmake.
if (${LANCET_BUILD_STATIC})
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
endif ()

# ── Clang Runtime Stack ───────────────────────────────────────────────────────
# Use clang's full runtime stack so the toolchain stays internally consistent:
# libc++ (instead of GCC's libstdc++), compiler-rt for builtins (instead of
# libgcc), libunwind for stack unwinding (instead of libgcc_s), and lld as the
# linker (instead of GNU ld). Vendored deps fetched via FetchContent inherit
# these flags automatically since they use the parent project's CMAKE_C/CXX_FLAGS.
#
# Why each piece:
#  -stdlib=libc++          Use clang's libc++ as the C++ standard library. Pinned
#                          in pixi.toml as the `libcxx` package.
#  -rtlib=compiler-rt      Use clang's runtime support library (128-bit math
#                          helpers, soft-float, etc.) instead of libgcc. This is
#                          the default on conda-forge clangxx_linux-64, but
#                          Debian's apt-installed clang (used in our Dockerfile)
#                          defaults to -rtlib=libgcc. Pin explicitly so both
#                          environments produce the same binary.
#  --unwindlib=libgcc      Use GCC's unwinder (libgcc_eh.a for static,
#                          libgcc_s.so for dynamic) for the C++ exception
#                          handling ABI (_Unwind_Resume, _Unwind_GetIP, etc.).
#                          Needed because -rtlib=compiler-rt replaces libgcc
#                          builtins but NOT the unwinder — without this,
#                          -static builds fail with undefined _Unwind_* symbols.
#                          Unrelated to gperftools' stack-walking libunwind
#                          (which we avoid in favor of frame-pointer unwinding).
#  -lc++abi                Explicitly link the C++ ABI library for static builds.
#                          Clang's driver with -stdlib=libc++ adds -lc++ but not
#                          -lc++abi because the shared libc++.so has ABI symbols
#                          merged in. The static libc++.a does not — it references
#                          __cxa_throw, RTTI vtables, etc. from libc++abi.a.
#  -lpthread               Transitive dependency of libc++abi.a, which uses
#                          pthreads for thread-safe exception handling.
#  -fuse-ld=lld            LLVM's linker. More reliable than GNU ld for clang
#                          static builds on conda-forge, faster link times,
#                          better diagnostics. The reference linker for Chromium,
#                          Firefox, and LLVM's own static builds.
#
# These apply to all clang builds (Release, Debug, RelWithDebInfo, sanitizer)
# regardless of LANCET_BUILD_STATIC, so the dynamic Docker build and the
# static single-binary distribution use the same C++ runtime.
#
# -stdlib=libc++ goes on both compile and link lines (compile: header search
#  paths; link: selects libc++.so/libc++.a). The rest are linker-only.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
	string(APPEND CMAKE_EXE_LINKER_FLAGS
		" -stdlib=libc++ -rtlib=compiler-rt --unwindlib=libgcc -lc++abi -lpthread -fuse-ld=lld")
endif ()

# ── Frame Pointers (all build modes) ─────────────────────────────────────────
# Preserve frame pointers in every build type (Debug, Release, RelWithDebInfo).
# Costs 1 register (~1-2% on x86_64) but enables:
#   - gperftools' frame-pointer stack unwinder (generic_fp), the only crash-safe
#     unwinder for CPU profiling from a signal handler
#   - Reliable backtraces in the crash handler, perf, GDB, and Instruments
#
# -momit-leaf-frame-pointer reclaims the RBP register in leaf functions (functions
# that don't call other functions), which are the tightest inner loops where
# register pressure matters most (DP matrices, k-mer hashing).  Profilers can
# still attribute time to leaf functions via the instruction pointer, and the
# generic_fp (safe variant) stops cleanly at the garbage RBP via bounds checking.
# gperftools itself uses this same flag combination since v2.11.
#
# Set on CMAKE_C/CXX_FLAGS (base flags) so it applies unconditionally.
# Also included in LANCET_OPT_FLAGS below to propagate to external dependencies
# (spoa, mimalloc, zlib-ng, htslib, etc.) via their configure scripts.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-omit-frame-pointer -momit-leaf-frame-pointer")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -momit-leaf-frame-pointer")

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
# -fno-omit-frame-pointer -momit-leaf-frame-pointer
#   Duplicated here from the base CMAKE_C/CXX_FLAGS above so that external
#   dependencies also get frame pointers when their configure scripts read
#   LANCET_OPT_FLAGS from CMakeCache.txt.  See the "Frame Pointers" section.
#
# Architecture flags (selected automatically per platform):
#   ARM64:  always -mcpu=native (ARM binaries are rarely cross-distributed)
#   x86_64: -march=x86-64-v3 (portable, AVX2 + BMI2, Haswell 2013+)
#           -march=native    (non-portable, when -DLANCET_NATIVE_BUILD=ON)

string(CONCAT LANCET_BASE_OPT "-O3 -pipe -fno-math-errno -fno-trapping-math"
	" -fno-omit-frame-pointer -momit-leaf-frame-pointer -DNDEBUG")

if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64" OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "aarch64")
	# ARM64: always native — ARM binaries target the build host's microarchitecture
	set(LANCET_OPT_FLAGS "${LANCET_BASE_OPT} -mcpu=native" CACHE STRING "Lancet optimization flags" FORCE)
elseif (LANCET_NATIVE_BUILD)
	# x86 native: uses all instruction sets available on this CPU (non-portable)
	message(STATUS "LANCET_NATIVE_BUILD=ON: using -march=native (non-portable binary)")
	set(LANCET_OPT_FLAGS "${LANCET_BASE_OPT} -march=native" CACHE STRING "Lancet optimization flags" FORCE)
else ()
	# x86 portable: AVX2 + BMI2 baseline (covers ~95% of modern server CPUs)
	set(LANCET_OPT_FLAGS "${LANCET_BASE_OPT} -march=x86-64-v3" CACHE STRING "Lancet optimization flags" FORCE)
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
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} ${LANCET_OPT_FLAGS}" CACHE INTERNAL "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${LANCET_OPT_FLAGS}" CACHE INTERNAL "" FORCE)

# ── Apply flags to RelWithDebInfo configuration ──────────────────────────────
# Same optimization flags as Release, plus -g for DWARF debug symbols.
# Used exclusively with LANCET_PROFILE_MODE=ON for source-level profiling
# via pprof --list. The -g flag adds ~2× binary size but zero runtime cost.
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS} ${LANCET_OPT_FLAGS} -g" CACHE INTERNAL "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${LANCET_OPT_FLAGS} -g" CACHE INTERNAL "" FORCE)

# ── Debug Arch ────────────────────────────────────────────────────────────────
# Apply -march / -mcpu to Debug builds. The Release / RelWithDebInfo blocks
# above get their arch flag from LANCET_OPT_FLAGS, but Debug builds otherwise
# inherit only CMake's default Debug flags, which include no arch flag at all.
# Without this, a Debug-built Lancet2 might run on hardware where the Release
# binary can't (or vice versa on ARM), masking codegen-divergence bugs that
# only surface under the platform's full ISA.
#
# Mirrors the LANCET_OPT_FLAGS branch structure above:
#   ARM64:  always -mcpu=native (matches Release)
#   x86_64: -march=x86-64-v3 by default, -march=native if LANCET_NATIVE_BUILD
#
# Sanitizer tasks (configure-asan, configure-tsan, configure-ubsan) override
# CMAKE_CXX_FLAGS wholesale in pixi.toml — they pass -march=x86-64-v3 directly
# in their CXX_FLAGS line. configure-msan deliberately stays at the SSE
# baseline because of historical AVX2-MSan instrumentation gaps.
if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm64" OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "aarch64")
	set(LANCET_DEBUG_ARCH_FLAG "-mcpu=native")
elseif (LANCET_NATIVE_BUILD)
	set(LANCET_DEBUG_ARCH_FLAG "-march=native")
else ()
	set(LANCET_DEBUG_ARCH_FLAG "-march=x86-64-v3")
endif ()
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${LANCET_DEBUG_ARCH_FLAG}" CACHE INTERNAL "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${LANCET_DEBUG_ARCH_FLAG}" CACHE INTERNAL "" FORCE)

# ── Link-Time Optimization ────────────────────────────────────────────────────
# Optional IPO. Do not use IPO if it's not supported by compiler.
include(CheckIPOSupported)
check_ipo_supported(RESULT COMPILER_SUPPORTS_IPO LANGUAGES C CXX)
if (COMPILER_SUPPORTS_IPO AND ${CMAKE_BUILD_TYPE} MATCHES Release)
	message(STATUS "Enabling Link Time Optimization because compiler supports it")
	set(ENABLE_LTO ON)
endif ()
