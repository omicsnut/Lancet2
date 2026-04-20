# ═══════════════════════════════════════════════════════════════════════════════
# Third-Party Dependencies — fetched and built from source
#
# All dependencies are pinned to exact versions for reproducible builds.
# FetchContent downloads at configure time; ExternalProject builds at build time.
#
# Core libraries:
#   mimalloc      — high-performance allocator (replaces system malloc)
#   abseil-cpp    — Abseil C++ utilities (containers, strings, hashing)
#   spdlog        — structured logging (bundles fmtlib)
#   CLI11         — command-line argument parsing
#   concurrentqueue — lock-free multi-producer/consumer queue
#
# Compression (required by HTSlib):
#   libdeflate    — fast DEFLATE/gzip compression
#   zlib-ng       — zlib-compatible compression (SIMD-optimized)
#
# Cloud I/O (optional, dynamic-linked; enabled via -DLANCET_ENABLE_CLOUD_IO=ON):
#   libcurl       — S3/GCS streaming (system, found via find_package)
#   OpenSSL       — TLS for cloud transport (system, found via find_package)
#
# Bioinformatics:
#   htslib        — BAM/CRAM/VCF I/O (ExternalProject, builds libhts.a)
#   minimap2      — read-to-haplotype alignment (ExternalProject, builds libminimap2.a)
#   WFA2-lib      — wavefront gap-affine alignment (ExternalProject, builds libwfacpp.a)
#   spoa          — SIMD Partial Order Alignment (graph-based MSA)
#
# Testing / Benchmarking / Profiling:
#   Catch2        — unit test framework (amalgamated, tests only)
#   benchmark     — Google Benchmark (benchmarks only)
#   gperftools    — CPU profiler (profile mode only)
# ═══════════════════════════════════════════════════════════════════════════════
include(ExternalProject)
include(FetchContent)
include(ProcessorCount)
ProcessorCount(NumCores)
find_program(MAKE_EXE NAMES gmake nmake make REQUIRED)

# Suppress CMake developer warnings from third-party FetchContent projects.
# Many dependencies (abseil, zlib-ng, spoa) emit CMP* policy warnings that
# clutter Lancet2's configure output and are not actionable from our side.
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON)

set(MI_SECURE OFF)
set(MI_PADDING OFF)
set(MI_BUILD_STATIC ON)
set(MI_BUILD_SHARED OFF)
set(MI_BUILD_OBJECT OFF)
set(MI_BUILD_TESTS OFF)
set(MI_OVERRIDE ON)
if (CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
	set(MI_OVERRIDE OFF)
endif ()
FetchContent_Declare(mimalloc GIT_REPOSITORY https://github.com/microsoft/mimalloc.git GIT_TAG v3.3.0 SYSTEM)
FetchContent_MakeAvailable(mimalloc)

FetchContent_Declare(abseil GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git GIT_TAG b9536c9 SYSTEM)
FetchContent_GetProperties(abseil)
if (NOT abseil_POPULATED)
	set(BUILD_TESTING OFF)
	set(ABSL_PROPAGATE_CXX_STD ON)
	set(ABSL_USE_SYSTEM_INCLUDES ON)
	FetchContent_Populate(abseil)
	add_subdirectory(${abseil_SOURCE_DIR} ${abseil_BINARY_DIR} SYSTEM)
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${abseil_SOURCE_DIR}/absl/copts)
	include(${abseil_SOURCE_DIR}/absl/copts/AbseilConfigureCopts.cmake)
endif ()

FetchContent_Declare(spdlog GIT_REPOSITORY https://github.com/gabime/spdlog.git GIT_TAG v1.17.0 SYSTEM)
FetchContent_MakeAvailable(spdlog)

FetchContent_Declare(cli11 GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git GIT_TAG v2.6.2 SYSTEM)
FetchContent_MakeAvailable(cli11)

FetchContent_Declare(concurrentqueue GIT_REPOSITORY https://github.com/cameron314/concurrentqueue.git GIT_TAG v1.0.5 SYSTEM)
FetchContent_GetProperties(concurrentqueue)
if (NOT concurrentqueue_POPULATED)
	FetchContent_Populate(concurrentqueue)
	add_library(concurrentqueue INTERFACE)
	target_include_directories(concurrentqueue SYSTEM INTERFACE "${concurrentqueue_SOURCE_DIR}")
endif ()

set(LIBDEFLATE_BUILD_STATIC_LIB ON)
set(LIBDEFLATE_BUILD_SHARED_LIB OFF)
set(LIBDEFLATE_BUILD_GZIP OFF)
set(LIBDEFLATE_BUILD_TESTS OFF)
set(LIBDEFLATE_USE_SHARED_LIB OFF)
FetchContent_Declare(libdeflate GIT_REPOSITORY https://github.com/ebiggers/libdeflate.git GIT_TAG v1.25 SYSTEM)
FetchContent_MakeAvailable(libdeflate)

set(ZLIB_COMPAT ON)
set(ZLIB_ENABLE_TESTS OFF)
set(ZLIBNG_ENABLE_TESTS OFF)
set(WITH_GTEST OFF)
set(WITH_BENCHMARKS OFF)
set(BUILD_SHARED_LIBS OFF)
set(WITH_NEW_STRATEGIES ON)
set(WITH_OPTIM ON)
set(WITH_NATIVE_INSTRUCTIONS OFF)
FetchContent_Declare(zlib-ng GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git GIT_TAG 2.3.3 SYSTEM)
FetchContent_MakeAvailable(zlib-ng)

set(HTSLIB_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/htslib")
set(LIB_HTS "${HTSLIB_ROOT_DIR}/libhts.a")
set(HTSLIB_CONFIG_PARAMS ${HTSLIB_ROOT_DIR} ${CMAKE_C_COMPILER} ${LANCET_ENABLE_CLOUD_IO})
ExternalProject_Add(htslib
		URL https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2
		URL_MD5 aee2c757fd8c88b9b5b61e8a1eae99de PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
		SOURCE_DIR ${HTSLIB_ROOT_DIR} BUILD_IN_SOURCE 1 INSTALL_COMMAND ""
		BUILD_COMMAND ${MAKE_EXE} -j${NumCores} lib-static BUILD_BYPRODUCTS ${LIB_HTS}
		CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/configure_htslib.sh ${HTSLIB_CONFIG_PARAMS}
		LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON USES_TERMINAL_DOWNLOAD OFF
		USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)
add_dependencies(htslib zlibstatic libdeflate_static)

set(MM2_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/minimap2")
set(LIB_MM2 "${MM2_ROOT_DIR}/libminimap2.a")
set(MM2_BUILD_PARAMS ${MM2_ROOT_DIR} ${CMAKE_C_COMPILER} ${CMAKE_HOST_SYSTEM_PROCESSOR})
ExternalProject_Add(minimap2
		URL https://github.com/lh3/minimap2/releases/download/v2.30/minimap2-2.30.tar.bz2
		URL_MD5 e016e3578bf6c763cefe08e7f22f440c PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps"
		SOURCE_DIR ${MM2_ROOT_DIR} BUILD_IN_SOURCE 1 CONFIGURE_COMMAND "" INSTALL_COMMAND ""
		BUILD_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/build_minimap2.sh ${MM2_BUILD_PARAMS}
		BUILD_BYPRODUCTS ${LIB_MM2} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)
add_dependencies(minimap2 zlibstatic)

set(spoa_optimize_for_native OFF)
FetchContent_Declare(spoa GIT_REPOSITORY https://github.com/rvaser/spoa GIT_TAG 4.1.5 SYSTEM)
FetchContent_MakeAvailable(spoa)

# ExternalProject (not FetchContent) because WFA2's CMakeLists.txt and Makefile
# both inject -march=native, which would pollute Lancet2's portable arch flags.
set(WFA2_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/WFA2-lib")
set(LIB_WFA2 "${WFA2_ROOT_DIR}/lib/libwfa.a")
set(LIB_WFA2CPP "${WFA2_ROOT_DIR}/lib/libwfacpp.a")
set(WFA2_BUILD_PARAMS ${WFA2_ROOT_DIR} ${CMAKE_C_COMPILER} ${CMAKE_CXX_COMPILER})
ExternalProject_Add(wfa2
		GIT_REPOSITORY https://github.com/smarco/WFA2-lib.git GIT_TAG v2.3.6
		PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${WFA2_ROOT_DIR}
		BUILD_IN_SOURCE 1 CONFIGURE_COMMAND "" INSTALL_COMMAND ""
		BUILD_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/build_wfa2.sh ${WFA2_BUILD_PARAMS}
		BUILD_BYPRODUCTS ${LIB_WFA2} ${LIB_WFA2CPP}
		LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)

if (LANCET_PROFILE_MODE)
	set(GPERFTOOLS_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/gperftools")
	set(GPERFTOOLS_INC_DIR "${GPERFTOOLS_ROOT_DIR}/include")
	set(LIB_PROFILER "${GPERFTOOLS_ROOT_DIR}/lib/libprofiler.a")
	set(GPERFTOOLS_CONFIG_PARAMS ${GPERFTOOLS_ROOT_DIR} ${CMAKE_C_COMPILER} ${CMAKE_CXX_COMPILER})
	ExternalProject_Add(gperftools
		URL https://github.com/gperftools/gperftools/releases/download/gperftools-2.18.1/gperftools-2.18.1.tar.gz
		URL_MD5 129c01f6f5297f0482b33d431b5ec555
		PREFIX "${CMAKE_CURRENT_BINARY_DIR}/_deps" SOURCE_DIR ${GPERFTOOLS_ROOT_DIR} BUILD_IN_SOURCE 1
		INSTALL_COMMAND ${MAKE_EXE} install BUILD_COMMAND ${MAKE_EXE} -j${NumCores}
		CONFIGURE_COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/configure_gperftools.sh ${GPERFTOOLS_CONFIG_PARAMS}
		BUILD_BYPRODUCTS ${LIB_PROFILER} LOG_DOWNLOAD ON LOG_CONFIGURE ON LOG_BUILD ON LOG_INSTALL ON
		USES_TERMINAL_DOWNLOAD OFF USES_TERMINAL_BUILD OFF USES_TERMINAL_INSTALL OFF)
endif ()

if (LANCET_TESTS)
	file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_deps/Catch2")
	set(CATCH_ROOT "${CMAKE_CURRENT_BINARY_DIR}/_deps/Catch2")
	set(CATCH_URL "https://github.com/catchorg/Catch2/releases/download/v3.14.0")
	set(CATCH_MD5c "c7c5431b9bef8a27bfad0a8a45581ba7")
	set(CATCH_MD5h "a0e371008e8fd95bab1ff2fc8ff2058c")
	file(DOWNLOAD "${CATCH_URL}/catch_amalgamated.cpp" "${CATCH_ROOT}/catch_amalgamated.cpp" EXPECTED_MD5 ${CATCH_MD5c})
	file(DOWNLOAD "${CATCH_URL}/catch_amalgamated.hpp" "${CATCH_ROOT}/catch_amalgamated.hpp" EXPECTED_MD5 ${CATCH_MD5h})
	add_library(Catch2 STATIC "${CATCH_ROOT}/catch_amalgamated.cpp" "${CATCH_ROOT}/catch_amalgamated.hpp")
	target_include_directories(Catch2 SYSTEM PUBLIC "${CATCH_ROOT}")
endif ()

if (LANCET_BENCHMARKS)
	set(BENCHMARK_ENABLE_TESTING OFF)
	set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
	set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF)
	set(BENCHMARK_ENABLE_INSTALL OFF)
	set(BENCHMARK_ENABLE_INSTALL OFF)
	set(BENCHMARK_INSTALL_DOCS OFF)
	set(BENCHMARK_ENABLE_LTO OFF)
	FetchContent_Declare(benchmark GIT_REPOSITORY https://github.com/google/benchmark.git GIT_TAG v1.9.5 SYSTEM)
	FetchContent_MakeAvailable(benchmark)
endif ()
