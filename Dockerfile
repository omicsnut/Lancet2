# ─── Builder stage ──────────────────────────────────────────────────────────
# Use pixi as the build-time toolchain manager so the compiler, libc++,
# libcurl, openssl, and every other library version exactly match what
# pixi-local and conda builds resolve to. This eliminates version drift
# between the three deployment paths (pixi, conda, docker).
#
# The pixi env is created inside the Docker builder layer at
# /Lancet2/.pixi/. We copy pixi.toml + pixi.lock first (separate Docker
# layer) so the slow `pixi install` step is cached when only source
# files change.
FROM debian:stable-slim AS builder
LABEL maintainer="Rajeeva Musunuri <rmusunuri@nygenome.org>"

# Minimal apt deps to bootstrap pixi. Everything else (clang, lld,
# libc++, libcurl, openssl, bzip2, lzma, cmake, ninja) comes from
# the pixi env, pinned to versions in pixi.toml.
RUN DEBIAN_FRONTEND="noninteractive" apt-get update && \
    apt-get install --yes --no-install-recommends \
        ca-certificates curl git xz-utils && \
    curl -fsSL https://pixi.sh/install.sh | bash && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
ENV PATH="/root/.pixi/bin:${PATH}"

WORKDIR /Lancet2

# Copy pixi.toml + lockfile FIRST (before source) so Docker layer
# caching keeps the slow `pixi install` step unchanged when only
# source files are modified. --frozen ensures the lockfile is used
# verbatim; CI fails fast if pixi.lock is stale relative to pixi.toml.
COPY pixi.toml pixi.lock ./
RUN pixi config set --local run-post-link-scripts insecure && \
    pixi install --frozen

# Now copy the rest of the source.
COPY . .

# Build with the same CMake flags as conda (STATIC=OFF + CLOUD_IO=ON),
# under pixi's env. Compiler (clang 22.1.4), C++ stdlib (libc++ 22.1.4),
# libcurl, openssl, etc. are all pinned by pixi.toml.
# CPU baseline -march=x86-64-v3 (Haswell 2013+: AVX2, BMI2, FMA, LZCNT)
# is the default in cmake/optimization_flags.cmake when LANCET_NATIVE_BUILD
# is OFF; we rely on that default here rather than passing an explicit
# arch flag. Pre-Haswell hardware will fail with an Illegal Instruction
# (SIGILL) crash on the first AVX2 instruction, typically during dynamic
# library init before main() runs. The hardware requirement is documented
# in the conda recipe description, the Dockerfile is intended for the
# same modern-server deployment context.
RUN pixi run cmake -GNinja \
        -B cmake-build-docker \
        -DCMAKE_C_COMPILER=clang \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DCMAKE_BUILD_TYPE=Release \
        -DLANCET_BUILD_STATIC=OFF \
        -DLANCET_ENABLE_CLOUD_IO=ON && \
    pixi run cmake --build cmake-build-docker -v

# Stage the binary and the .so files it depends on from the pixi env.
# Lancet2 with STATIC=OFF + CLOUD_IO=ON dynamically links against
# libc++.so.1, libc++abi.so.1, libcurl.so.4, libssl.so, libcrypto.so,
# libz.so.1, libbz2.so.1, liblzma.so.5, libunwind.so.1, and a few
# libcurl transitive deps (nghttp2, etc.). All come from the pixi env.
# We use ldd inside `pixi run` to find them, then cp -L to follow
# symlinks so the .so files (not just the version symlinks) end up
# in the staging dir.
RUN mkdir -p /staging/bin /staging/lib && \
    cp cmake-build-docker/Lancet2 /staging/bin/Lancet2 && \
    pixi run ldd cmake-build-docker/Lancet2 | \
        awk '/=>.*\.pixi\// {print $3}' | \
        sort -u | \
        xargs -I{} cp --no-clobber -L {} /staging/lib/ && \
    echo "── Bundled runtime libraries from pixi env ──" && \
    ls -la /staging/lib/

# ─── Runtime stage ──────────────────────────────────────────────────────────
# Debian stable-slim base. The bundled libs from the pixi env (libc++,
# libcurl, openssl, etc.) live at /opt/lancet2-libs and are reached via
# LD_LIBRARY_PATH. The host system's libc/libpthread/libm/libdl/libgcc_s
# resolve from Debian stable-slim — these are the basic C runtime that
# every Linux binary needs and they're stable across glibc versions for
# the symbols Lancet2 uses. (Builder and runtime stages are both Debian
# stable so their glibc matches by construction.)
#
# Image size: this approach trades a slightly bigger image (~90-100 MB
# vs ~80 MB previously) for exact toolchain match with pixi-local and
# conda builds. The bundled libs total ~10-15 MB.
FROM debian:stable-slim
RUN DEBIAN_FRONTEND="noninteractive" apt-get update && \
    apt-get install --yes --no-install-recommends \
        ca-certificates bash && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY --from=builder /staging/bin/Lancet2 /usr/bin/Lancet2
COPY --from=builder /staging/lib /opt/lancet2-libs

# CMake's default RPATH for build trees points at the build-time lib
# location (/Lancet2/.pixi/envs/default/lib). That path doesn't exist
# in the runtime image, so the dynamic loader falls through to
# LD_LIBRARY_PATH. CMake uses DT_RUNPATH (not DT_RPATH), so
# LD_LIBRARY_PATH takes precedence over the stale build-time RPATH.
ENV LD_LIBRARY_PATH=/opt/lancet2-libs:$LD_LIBRARY_PATH
