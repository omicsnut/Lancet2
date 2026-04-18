# ═══════════════════════════════════════════════════════════════════════════════
# Platform & Compiler Validation
#
# Enforces build requirements before any targets are defined:
#   - OS:       Linux or macOS (no Windows, no WSL)
#   - Arch:     x86_64 / amd64 or arm64 / aarch64 (no 32-bit)
#   - Compiler: GCC ≥ 12.0 or Clang ≥ 14.0 (C++20 concepts + ranges)
#   - Static:   macOS forces dynamic linking (Apple linker limitation)
#   - Default:  Release build type when none is specified
# ═══════════════════════════════════════════════════════════════════════════════

# ── OS check ──────────────────────────────────────────────────────────────────
if (NOT (CMAKE_HOST_SYSTEM_NAME MATCHES "Linux" OR CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin"))
	message(WARNING "Running Lancet build on a ${CMAKE_HOST_SYSTEM_NAME} machine")
	message(FATAL_ERROR "Lancet can only be built on systems with Linux or macOS Operating System")
endif ()

# ── Architecture check ────────────────────────────────────────────────────────
if (NOT (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "amd64" OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64" OR CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "aarch64"))
	message(WARNING "Running Lancet build on a ${CMAKE_HOST_SYSTEM_PROCESSOR} machine")
	message(FATAL_ERROR "Lancet can only be built on a 64-bit x86 or ARM64 machine")
endif ()

# ── Compiler vendor check ────────────────────────────────────────────────────
if (NOT ${CMAKE_C_COMPILER_ID} MATCHES "GNU|Clang|AppleClang" OR NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "GNU|Clang|AppleClang")
	message(WARNING "Running Lancet build with ${CMAKE_C_COMPILER} and ${CMAKE_CXX_COMPILER}")
	message(FATAL_ERROR "Lancet can only be built from source with GCC or Clang compilers")
endif ()

# ── Compiler version check ───────────────────────────────────────────────────
# GCC 12+ required for C++20 concepts, ranges, and std::format support.
if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU" AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "12.0")
	message(WARNING "Running Lancet build with GCC version ${CMAKE_CXX_COMPILER_VERSION}")
	message(FATAL_ERROR "Lancet requires GCC version to be at least 12.0 or greater")
endif ()

# Clang 14+ required for C++20 concepts and ranges constraints.
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang|AppleClang" AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "14.0")
	message(WARNING "Running Lancet build with Clang version ${CMAKE_CXX_COMPILER_VERSION}")
	message(FATAL_ERROR "Lancet requires Clang version to be at least 14.0 or greater")
endif ()

# ── macOS static build override ──────────────────────────────────────────────
# Apple's linker does not support fully static executables.
if (CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin" AND LANCET_BUILD_STATIC)
	message(WARNING "macOS does not support fully static builds. Forcing LANCET_BUILD_STATIC to OFF.")
	set(LANCET_BUILD_STATIC OFF CACHE STRING "Build a statically linked Lancet executable" FORCE)
endif ()

# ── Default build type ───────────────────────────────────────────────────────
if (NOT LANCET_PROFILE_MODE AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting CMAKE_BUILD_TYPE as Release since none was specified.")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

# ── Profiling build type override ────────────────────────────────────────────
# LANCET_PROFILE_MODE requires RelWithDebInfo: retains -O2 optimized codegen
# while preserving DWARF debug info for pprof's --list source-line annotation.
# Release strips debug info, making sub-function profiling impossible.
if (LANCET_PROFILE_MODE AND NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
	message(STATUS "LANCET_PROFILE_MODE=ON requires RelWithDebInfo for source-level profiling.")
	message(STATUS "Overriding CMAKE_BUILD_TYPE from '${CMAKE_BUILD_TYPE}' to 'RelWithDebInfo'.")
	set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
endif ()
