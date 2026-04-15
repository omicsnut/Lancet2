# ═══════════════════════════════════════════════════════════════════════════════
# Build Toolchain — CMake policies, ccache, Link-Time Optimization
#
# Configures build infrastructure that must be set before any targets:
#   - CMake policy defaults for FetchContent subprojects
#   - ccache integration (compile + link caching when available)
#   - Link-Time Optimization (Release builds only, if compiler supports it)
# ═══════════════════════════════════════════════════════════════════════════════

# ── Policy Defaults for FetchContent Subprojects ──────────────────────────────
# With cmake_minimum_required(VERSION 3.25), all these policies are already NEW
# in Lancet2's own scope. These CMAKE_POLICY_DEFAULT_* variables force NEW
# behavior on FetchContent-downloaded subprojects that may declare a lower
# cmake_minimum_required (e.g., abseil, spoa, zlib-ng).
#
# CMP0012: if() recognizes numbers and boolean constants
# CMP0048: project() supports VERSION keyword
# CMP0063: Honor visibility properties for all target types
# CMP0069: Enforce INTERPROCEDURAL_OPTIMIZATION (LTO)
# CMP0074: find_package uses <PackageName>_ROOT variables
# CMP0077: option() honors normal variables set before it
# CMP0135: FetchContent sets file timestamps to extraction time
set(CMAKE_POLICY_DEFAULT_CMP0012 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0063 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0074 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0135 NEW)

# ── ccache ────────────────────────────────────────────────────────────────────
# Use ccache if found to cache previously built object files
find_program(CCACHE_EXE ccache)
if (CCACHE_EXE)
	message(STATUS "Found ccache in PATH. Using ccache to speed up recompilation")
	set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ${CCACHE_EXE})
	set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ${CCACHE_EXE})
endif ()

# ── Link-Time Optimization ────────────────────────────────────────────────────
# Optional IPO. Do not use IPO if it's not supported by compiler.
include(CheckIPOSupported)
check_ipo_supported(RESULT COMPILER_SUPPORTS_IPO LANGUAGES C CXX)
if (COMPILER_SUPPORTS_IPO AND ${CMAKE_BUILD_TYPE} MATCHES Release)
	message(STATUS "Enabling Link Time Optimization because compiler supports it")
	set(ENABLE_LTO ON)
endif ()
